function [IFR_PSTH_Summary, all_IFR] = run_IFR_PSTH_byDepthMove(All_SpikesPerMove_Tbl, AO_spike_fs, varargin)

% Compute instantaneous firing rate (IFR; via getIFR.m) and PSTH for each MoveType × STN depth
% (and for each SpikeField C# if present).
% Uses stitched inputs from makeZetaInputs_fromAOStartStopTimes (true per-trial durations).
%
% OUTPUTS
%   IFR_PSTH_Summary : table with one row per (MoveType × Depth)
%   all_IFR          : cell array of structs containing arrays & metadata per row

p = inputParser;
% Data/schema/options (mirror ZETA wrapper defaults)
p.addParameter('UseMaxDur_s', [], @(x) isempty(x) || isscalar(x));
p.addParameter('PadITI_s', 0.005, @(x) isscalar(x) && x>=0);                % 5 ms gap inserted between trials on the stitched timeline (prevents overlap)
p.addParameter('MinDur_s', 0,        @(x) isscalar(x) && x>=0);
p.addParameter('MaxDur_s', inf,      @(x) isscalar(x) && x>0);

p.addParameter('DepthIDs', {'t','c','b'});
p.addParameter('MoveTypeOrder', {'HAND OC','HAND PS','ARM EF','REST'});

% spike fields: allow either SpikeField (single) or SpikeFields (multi/auto)
p.addParameter('SpikeField', 'C1', @(s) ischar(s) || isstring(s));
p.addParameter('SpikeFields', [], @(x) iscellstr(x) || isempty(x));

% start/stop bounds for spike field segments (based on align_SpikesPerMove_TTL)
p.addParameter('StartField', 'TTL_spk_idx_Start', @(s) ischar(s) || isstring(s)); % TTL_spk_idx_Start(move_i) - offset_spike_samples (50 ms), onset (t=0) minus 50ms
p.addParameter('EndField',  'TTL_spk_idx_End',   @(s) ischar(s) || isstring(s));  % from start of current rep to start of next rep (minus 10 frames), a full move rep

% optional lead-in/out (sec) padding:
p.addParameter('PreWindow_s',  0.050, @(x) isscalar(x) && x>=0);            % 50 ms pre-onset --> 150 ms
p.addParameter('PostWindow_s', 0.000, @(x) isscalar(x) && x>=0);            % 0  ms post-offset

% IFR-specific params (match getIFR signature)
p.addParameter('IFR_SmoothSd', 2, @(x) isscalar(x) && x>=0);   % Gaussian SD (#bins in MSD space)
p.addParameter('IFR_MinScale', [], @(x) isempty(x) || isscalar(x));
p.addParameter('IFR_Base',     1.5, @(x) isscalar(x) && x>1);
p.addParameter('IFR_UseParallel', [], @(x) isempty(x) || islogical(x));

% PSTH/plot/save
p.addParameter('BinSize_ms', 10, @(x) isscalar(x) && x>0);
p.addParameter('DoPlot', true, @islogical);
p.addParameter('SaveDir', '', @(s) ischar(s) || isstring(s));   % if non-empty, saves PNG + MAT
p.addParameter('CaseDate', '', @(s) ischar(s) || isstring(s));  % used in filenames
p.parse(varargin{:});
U = p.Results;


%% JNE Color scheme

JNE_colors = (...  % full palette
    [118,42,131;   % dark purple      (dorsal STN)
    175,141,195;   % lavender         (central STN)
    231,212,232;   % *light purple    (ventral STN)
    128,128,128;   % grey (REST / neutral)
    217,240,211;   % light green      (Hand OC)
    127,191,123;   % green            (Hand PS)
    27,120,55]) ... % dark green      (Arm EF)
    ./ 256; % rgb scaling

% Depth colors (purple shades)
purpleShades = ([ ...
    118,42,131;      % dorsal STN (t)
    175,141,195;     % central STN (c)
    231-15, 212-15, 232-15] ... % ventral STN (b)
    ./ 255); % /255 = standard

% Movement-context colors (greens)
greenShades = ([ ...
    128,128,128;     % REST  (grey)
    217-20,240-20,211-20;     % HAND OC  (light green)
    127,191,123;     % HAND PS  (green)
    27,120,55] ...   % ARM EF   (dark green)
    ./ 255);  % /255 = standard


% Maps for easy lookup
depthColorMap = containers.Map( ...
    {'t','c','b'}, ...
    {purpleShades(1,:), purpleShades(2,:), purpleShades(3,:)});

moveColorMap  = containers.Map( ...
    {'REST','HAND OC','HAND PS','ARM EF'}, ...
    {greenShades(1,:), greenShades(2,:), greenShades(3,:), greenShades(4,:)});


%% Initialize Spike fields and MoveTypes

% Normalize spike fields: prefer SpikeFields; else use SpikeField
spikeFields = U.SpikeFields;
if isempty(spikeFields)
    if ~isempty(U.SpikeField)
        spikeFields = {char(U.SpikeField)};
    else
        cand = All_SpikesPerMove_Tbl.Properties.VariableNames;
        isCnum = ~cellfun('isempty', regexp(cand,'^C\d+$','once'));
        spikeFields = cand(isCnum);
        keep = false(size(spikeFields));
        for k = 1:numel(spikeFields)
            col = All_SpikesPerMove_Tbl.(spikeFields{k});
            keep(k) = iscell(col) && any(~cellfun(@isempty,col));
        end
        spikeFields = spikeFields(keep);
    end
end

move_types = intersect(U.MoveTypeOrder, unique(All_SpikesPerMove_Tbl.MoveType),'stable');


%% Outer loop over spike fields (C#s)

rows = {};
all_IFR = {};

for SpkF = spikeFields
    curSF = SpkF{1};

    % Loop over MoveType × STN depth
    for m = 1:numel(move_types)
        for d = 1:numel(U.DepthIDs)
            mv = move_types{m};
            dz = U.DepthIDs{d};

            % subset rows
            move_tbl = All_SpikesPerMove_Tbl( ...
                strcmp(All_SpikesPerMove_Tbl.MoveType, mv) & ...
                contains(All_SpikesPerMove_Tbl.move_trial_ID, dz), :);

            if isempty(move_tbl), continue; end

            % --- Behavioral repetition count (unit/MUA independent) ---
            Start_idx = char(U.StartField);
            End_idx   = char(U.EndField);

            hasStartStop = ~isnan(move_tbl.(Start_idx)) & ~isnan(move_tbl.(End_idx)) & ...
                (move_tbl.(End_idx) > move_tbl.(Start_idx));

            repIdx_behavioral = find(hasStartStop);  % indices into move_tbl

            % Store behavioral rep count alongside neural-valid trial count
            nReps_behavioral = sum(hasStartStop);    % true rep count for y-axis


            %% ===== Helper Functions =========================
            % ------ Build stitched spike/event inputs --------
            % Rep-centric: do NOT drop empty spike trials (for consistent raster rep counts across SU vs MUA)
            if strcmp(mv, 'REST')
                % REST-specific behavior (new):
                [spkT, evTimes, useMaxDur, trialDur, kept_tbl] = makeZetaInputs_fromAOStartTimes_REST( ...
                    move_tbl, AO_spike_fs, ...
                    'UseMaxDur_s', U.UseMaxDur_s, ...
                    'PadITI_s',    U.PadITI_s, ...
                    'SpikeField',  curSF, ...
                    'StartField',  U.StartField, ...
                    'EndField',    U.EndField, ...
                    'StartOffset_s', 0.250, ...
                    'MinDur_s',      0.750, ...
                    'MaxDur_s',      1.000, ...
                    'PreWindow_s',   0.000, ...
                    'PostWindow_s',  U.PostWindow_s, ...
                    'DropEmptySpikeTrials', false);
            else
                % Active movement behavior (unchanged):
                [spkT, evTimes, useMaxDur, trialDur, kept_tbl] = makeZetaInputs_fromAOStartStopTimes( ...
                    move_tbl, AO_spike_fs, ...
                    'UseMaxDur_s',  U.UseMaxDur_s, ...
                    'PadITI_s',     U.PadITI_s, ...
                    'SpikeField',   curSF, ...
                    'StartField',   U.StartField, ...
                    'EndField',     U.EndField, ...
                    'PreWindow_s',  U.PreWindow_s, ...
                    'PostWindow_s', U.PostWindow_s, ...
                    'DropEmptySpikeTrials', false);  % false: rep-centric, true: unit-centric
            end


            if isempty(spkT) || isempty(evTimes), continue; end
            nTrials = size(evTimes,1); % neural trials contributing to IFR/ZETA

            % Map neural trials back onto behavioral rep indices
            repIdx_neural = repIdx_behavioral(1:nTrials);


            %% --- compute IFR (getIFR works on stitched timeline like zetatest) ---
            [vecTime, vecRate, sIFR] = getIFR( ...
                spkT, evTimes, useMaxDur, ...
                U.IFR_SmoothSd, U.IFR_MinScale, U.IFR_Base, U.IFR_UseParallel);


            %% Plotting:

            % --- PSTH across trials (aligned to each event onset) ---
            bin_w   = U.BinSize_ms / 1000;              % seconds
            tmin    = -U.PreWindow_s;                   % true onset - 50 ms
            tmax    =  useMaxDur;                       % cap PSTH to the IFR window
            edges   = tmin:bin_w:tmax;
            centers = edges(1:end-1) + bin_w/2;

            % build raster (relative times) and histogram
            ras_t = [];
            ras_tr = [];
            counts = zeros(1, numel(edges)-1);

            for i = 1:nTrials
                on  = evTimes(i,1);
                hi  = on + tmax;
                lo  = on + tmin;
                mask = spkT >= lo & spkT < hi;
                rel  = spkT(mask) - on;

                if ~isempty(rel)
                    rel = rel(:);  % forces column

                    ras_t   = [ras_t; rel];

                    % Map neural trials back onto behavioral rep indices
                    % ras_tr  = [ras_tr; i*ones(numel(rel),1)]; % unit-centric
                    ras_tr = [ras_tr; repIdx_neural(i) .* ones(size(rel))]; % rep-centric
                    counts  = counts + histcounts(rel, edges);
                end
            end

            % PSTH bins from -PreWindow_s to useMaxDur (s), normalized to Hz
            psth_Hz = (counts ./ nTrials) ./ bin_w;      % spikes/trial/bin  -> Hz


            %% --- Compute mean IFR and max IFR in peri-movement window ---
            % Use [-PreWindow, UseMaxDur] (same as PSTH & ZETA window)
            IFR_mean_Hz = NaN;
            IFR_max_Hz  = NaN;
            if ~isempty(vecTime) && ~isempty(vecRate)
                t_vec = vecTime(:);
                r_vec = vecRate(:);
                mask_mean = (t_vec >= -U.PreWindow_s) & (t_vec <= useMaxDur);
                if any(mask_mean)
                    IFR_mean_Hz = mean(r_vec(mask_mean), 'omitnan'); % mean IFR in a defined epoch
                    IFR_max_Hz = max(r_vec(mask_mean), [], 'omitnan'); % max IFR in a defined epoch
                end
            end


            %% Plotting update
            % shift every x-axis quantity by -U.PreWindow_s and draw xline(0).

            % Make x=0 correspond to the true onset (+ U.PreWindow_s, 0+50ms)
            trueOnset_s = U.PreWindow_s;  % 0.05

            % Shift x-axis so 0 = true onset
            ras_t0    = ras_t - trueOnset_s;
            centers0  = centers - trueOnset_s;
            tmin0     = tmin - trueOnset_s;
            tmax0     = tmax - trueOnset_s;


            % --- Simple figures (raster + IFR/PSTH overlay) ---
            hFig = [];
            if U.DoPlot
                hFig = figure('Color','w','Position',[100 100 900 600]);
                tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

                % Depth label
                depthMap = containers.Map({'t','c','b'},{'dorsal STN','central STN','ventral STN'});
                if isKey(depthMap,dz); depthLbl = depthMap(dz); else, depthLbl = dz; end

                % Colors for this depth
                if isKey(depthColorMap, dz)
                    depthCol = depthColorMap(dz);
                else
                    depthCol = [0 0 0];
                end
                % Colors for this movement
                if isKey(moveColorMap, mv)
                    moveCol = moveColorMap(mv);
                else
                    moveCol = [0.3 0.3 0.3];
                end

                %% Raster
                ax1 = nexttile;
                if ~isempty(ras_t)
                    % raster dots colored by depth
                    scatter(ras_t0, ras_tr, 6, depthCol, 'filled'); hold on;
                end

                % time-zero line in neutral dark gray
                xline(0,'-','Color',[0.2 0.2 0.2]); % 0 = trueOnset_s
                grid on;

                % y-axis: integer ticks only
                ylim([0.5 nTrials+0.5]);
                yticks(1:nTrials);
                ylabel('Trial Number');

                % x-axis
                xlim([-U.PreWindow_s, useMaxDur]); % [tmin tmax]
                xlabel('Time from onset (s)');

                % Title format: Raster | MoveType - Depth (n=## reps)
                % title(sprintf('Raster | %s - %s - %s (n=%d reps)', mv, depthLbl, curSF, nTrials));
                sfLabel = strrep(curSF,'_','\_');  % escape '_' so it's not a subscript
                title(sprintf('Raster | %s - %s - %s (n=%d reps)', mv, depthLbl, sfLabel, nTrials));


                %% IFR + PSTH overlay (two y-axes)
                ax2 = nexttile;

                % Left y-axis: IFR
                yyaxis left
                if ~isempty(vecTime)
                    % vecTime is stitched absolute; convert to peristimulus using modulo onsets
                    % The IFR from getIFR is already averaged across trials on the stitched grid,
                    % but plotted against vecTime (stitched). For display, we re-index around onsets:
                    % Simpler: plot IFR against vecTime relative to first onset block (vecTime covers 0..useMaxDur for each block).
                    % The average IFR returned is arranged on sIFR.vecTime which repeats per trial block.
                    % To show as peri-stimulus, collapse by modulo useMaxDur (safe when PadITI_s > 0).
                    blockDur = useMaxDur + U.PadITI_s;

                    % peri-event time where 0 = padded start
                    relIFR_t_pad0 = mod(vecTime, blockDur);

                    % shift so that time runs from -PreWindow to useMaxDur
                    relIFR_t = mod(vecTime + U.PreWindow_s, blockDur) - U.PreWindow_s;

                    % shift so 0 = true onset
                    relIFR_t0 = relIFR_t_pad0 - U.PreWindow_s; 

                    keepIFR  = relIFR_t0 >= -U.PreWindow_s & relIFR_t0 <= useMaxDur; % [tmin, tmax]

                    plot(relIFR_t0(keepIFR), vecRate(keepIFR), 'LineWidth', 2, ...
                        'Color', depthCol);
                    hold on;
                end
                ylabel('IFR (Hz)');
                yl = ylim; ylim([0 max(yl(2), eps)]); % Force IFR y-axis to start at 0

                % Right y-axis: PSTH
                yyaxis right
                plot(centers0, psth_Hz, '--', 'LineWidth', 2, ...
                    'Color', moveCol);
                ylabel(sprintf('PSTH (Hz, %d ms bins)', U.BinSize_ms));
                y2 = ylim; ylim([0 max(y2(2), eps)]);  % Force PSTH y-axis to start at 0

                % Shared x settings
                xlim([-U.PreWindow_s, useMaxDur]); % [tmin tmax]
                % time-zero line in neutral dark gray
                xline(0,'-','Color',[0.2 0.2 0.2]); % 0 = trueOnset_s
                grid on;
                xlabel('Time from onset (s)');
                title('IFR (solid) + PSTH (dashed)');
            end

            %% --- Collect outputs ---
            rows(end+1,:) = { ...
                curSF, mv, dz, nTrials, nReps_behavioral, useMaxDur, ...
                U.BinSize_ms, mean(trialDur,'omitnan'), std(trialDur,'omitnan'), ...
                mean(vecRate,'omitnan'), max(vecRate,[],'omitnan'), ...
                IFR_mean_Hz, IFR_max_Hz, ... % explicitly windowed mean, max
                centers, centers0, psth_Hz, ...
                vecTime, vecRate};

            sOut = struct;
            sOut.SpikeField = curSF;
            sOut.MoveType = mv;
            sOut.Depth = dz;
            sOut.nTrials = nTrials;
            sOut.nReps_behavioral = nReps_behavioral;
            sOut.useMaxDur_s = useMaxDur;
            sOut.PreWindow_s = U.PreWindow_s;
            sOut.PostWindow_s = U.PostWindow_s;
            sOut.binSize_ms  = U.BinSize_ms;
            sOut.edges = edges; 
            sOut.centers = centers; sOut.centers_t0 = centers0;
            sOut.psth_Hz = psth_Hz;
            sOut.raster_t = ras_t; sOut.raster_t0 = ras_t0;
            sOut.raster_trial = ras_tr; 
            sOut.vecTime = vecTime; sOut.vecRate = vecRate; sOut.sIFR = sIFR;
            sOut.IFR_mean_Hz  = IFR_mean_Hz;
            sOut.IFR_max_Hz   = IFR_max_Hz;

            all_IFR{end+1,1} = sOut;

            % --- Save?
            if ~isempty(U.SaveDir)
                if ~exist(U.SaveDir,'dir'); mkdir(U.SaveDir); end
                tag = sprintf('%s_%s_%s', mv, dz, curSF);
                if ~isempty(U.CaseDate)
                    base = sprintf('%s_IFR-PSTH_%s', U.CaseDate, tag);
                else
                    base = sprintf('IFR-PSTH_%s', tag);
                end

                if ~isempty(hFig) && isvalid(hFig)
                    print(hFig, fullfile(U.SaveDir, [base,'.png']), '-dpng','-r300');
                    close(hFig);
                end
                save(fullfile(U.SaveDir, [base,'.mat']), '-struct','sOut','-v7.3');
            end

        end % depth
    end % move
end % spike field


%% ---- finalize table ----
varNames = { ...
    'SpikeField','MoveType','Depth', ...
    'nTrials', 'nReps_behavioral', 'UseMaxDur_s', ...
    'BinSize_ms','MeanTrialDur_s','StdTrialDur_s', ...
    'MeanIFR_Hz', 'MaxIFR_Hz', ...
    'IFR_mean_Hz', 'IFR_max_Hz',...
    'PSTH_TimeCenters_s','PSTH_TimeCenters_t0_s','PSTH_Hz', ...
    'IFR_Time_s','IFR_Hz'};

if isempty(rows)
    IFR_PSTH_Summary = cell2table(cell(0,numel(varNames)), 'VariableNames', varNames);
else
    IFR_PSTH_Summary = cell2table(rows, 'VariableNames', varNames);
end


%% Compute baseline IFR per unit × depth (REST trials only)
% and baseline-normalized IFR mean and IFR max

nRows = height(IFR_PSTH_Summary);
IFR_baseline_Hz = nan(nRows,1); % mean IFR of REST trials per unit x depth
IFR_mean_baselineNorm = nan(nRows,1); % baseline-normalized mean IFR
IFR_max_baselineNorm  = nan(nRows,1); % baseline-normalized max IFR

tbl = IFR_PSTH_Summary;

% Ensure grouping variables are categorical
if ~iscategorical(tbl.SpikeField)
    tbl.SpikeField = categorical(tbl.SpikeField);
end
if ~iscategorical(tbl.Depth)
    tbl.Depth = categorical(tbl.Depth);
end
if ~iscategorical(tbl.MoveType)
    tbl.MoveType = categorical(tbl.MoveType);
end

% epsilon floor
epsBase = 1e-6;

% Group by unit × depth
[G_unit_d, udLevels] = findgroups(tbl.SpikeField, tbl.Depth);

for g = 1:numel(udLevels)
    idx  = (G_unit_d == g);
    rest_rows = tbl(idx,:);

    % REST rows within this unit × depth
    rest_mask = rest_rows.MoveType == "REST";

    if ~any(rest_mask)
        % No REST baseline for this depth for this unit
        continue
    end

    % Mean REST IFR across ALL REST rows for this unit × depth
    rest_mean_vals   = rest_rows.IFR_mean_Hz(rest_mask);
    baseline_val = mean(rest_mean_vals, 'omitnan');

    % Guard against NaN, zero, or tiny baseline
    if isnan(baseline_val) || baseline_val == 0 || baseline_val < epsBase
        IFR_baseline_Hz(idx)  = NaN;
        IFR_mean_baselineNorm(idx) = NaN;
        IFR_max_baselineNorm(idx) = NaN;
    else
        % Map local index back to global indices in IFR_PSTH_Summary
        IFR_baseline_Hz(idx) = baseline_val;

        % Baseline-normalized: (IFR - baseline)/baseline
        IFR_mean_baselineNorm(idx) = (tbl.IFR_mean_Hz(idx) - baseline_val) ./ baseline_val;
        IFR_max_baselineNorm(idx) = (tbl.IFR_max_Hz(idx) - baseline_val) ./ baseline_val;
    end
end

% Attach to summary table
IFR_PSTH_Summary.IFR_baseline_Hz      = IFR_baseline_Hz;
IFR_PSTH_Summary.IFR_mean_baselineNorm     = IFR_mean_baselineNorm;
IFR_PSTH_Summary.IFR_max_baselineNorm      = IFR_max_baselineNorm;

% Mirror into all_IFR structs
for i = 1:nRows
    all_IFR{i}.IFR_baseline_Hz          = IFR_baseline_Hz(i);
    all_IFR{i}.IFR_mean_baselineNorm    = IFR_mean_baselineNorm(i);
    all_IFR{i}.IFR_max_baselineNorm     = IFR_max_baselineNorm(i);
end


%% Normalize mean IFR per unit x depth (within-unit z-score of IFR_mean_Hz per depth across unit's conditions)
% normalize within each unit per depth across all active movement types (exclude REST)

% exclude rest
isActive = tbl.MoveType ~= "REST";

IFR_mean_Znorm = nan(height(tbl),1);

% groups by unit x depth
for g = 1:numel(udLevels)
    idxGroup = (G_unit_d == g);

    % compute mu/sigma using ACTIVE rows only within group (unit × depth)
    idxFit_Active = idxGroup & isActive;
    vals = tbl.IFR_mean_Hz(idxFit_Active);
    if sum(~isnan(vals)) >= 2
        mu  = mean(vals, 'omitnan');
        sig = std(vals,  'omitnan');
        % apply z-score to all rows within group (including REST); REST becomes "relative to active distribution"
        if sig > 0
            IFR_mean_Znorm(idxGroup) = (tbl.IFR_mean_Hz(idxGroup) - mu) ./ sig;
        else
            IFR_mean_Znorm(idxGroup) = 0;
        end
    else
        % if you only have 0–1 active conditions for this unit×depth, define 0
        IFR_mean_Znorm(idxGroup) = 0;
    end
end

% Attach to summary table
IFR_PSTH_Summary.IFR_mean_Znorm = IFR_mean_Znorm;

% Mirror into all_IFR structs
for i = 1:height(tbl)
    all_IFR{i}.IFR_mean_Znorm = IFR_mean_Znorm(i);
end

end