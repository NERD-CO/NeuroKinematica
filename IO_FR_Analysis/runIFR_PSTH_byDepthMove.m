function [IFR_PSTH_Summary, all_IFR] = runIFR_PSTH_byDepthMove(All_SpikesPerMove_Tbl, AO_spike_fs, varargin)

% Compute instantaneous firing rate (IFR; via getIFR.m) and PSTH for each MoveType × STN depth.
% Uses stitched inputs from makeZetaInputs_fromAOStartStopTimes (true per-trial durations).
%
% OUTPUTS
%   IFR_PSTH_Summary : table with one row per (MoveType × Depth)
%   all_IFR          : cell array of structs containing arrays & metadata per row

p = inputParser;
% Data/schema/options (mirror ZETA wrapper defaults)
p.addParameter('UseMaxDur_s', [], @(x) isempty(x) || isscalar(x));
p.addParameter('PadITI_s', 0.005, @(x) isscalar(x) && x>=0);
p.addParameter('DepthIDs', {'t','c','b'});
p.addParameter('MoveTypeOrder', {'HAND OC','HAND PS','ARM EF','REST'});
% spike fields: allow either SpikeField (single) or SpikeFields (multi/auto)
p.addParameter('SpikeField', 'C1', @(s) ischar(s) || isstring(s));
p.addParameter('SpikeFields', [], @(x) iscellstr(x) || isempty(x));
% start/stop bounds for spike field segments:
p.addParameter('StartField', 'TTL_spk_idx_Start', @(s) ischar(s) || isstring(s));
p.addParameter('StopField',  'TTL_spk_idx_End',   @(s) ischar(s) || isstring(s)); %%%
% optional lead-in/out (sec) padding:
p.addParameter('PreWindow_s',  0.050, @(x) isscalar(x) && x>=0);
p.addParameter('PostWindow_s', 0.000, @(x) isscalar(x) && x>=0);

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


% Normalize spike fields: prefer SpikeFields; else use SpikeField; else auto-detect C\d+ with any spikes
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

rows = {};
all_IFR = {};

% Outer loop over spike fields (C#s)
for SpkF = spikeFields
    curSF = SpkF{1};

    % Loop over MoveType × STN depth
    for m = 1:numel(move_types)
        for d = 1:numel(U.DepthIDs)
            mv = move_types{m};
            dz = U.DepthIDs{d};

            move_tbl = All_SpikesPerMove_Tbl( ...
                strcmp(All_SpikesPerMove_Tbl.MoveType, mv) & ...
                contains(All_SpikesPerMove_Tbl.move_trial_ID, dz), :);

            if isempty(move_tbl), continue; end

            % Helper function:
            % --- Stitched inputs (true [start stop], optional pre/post kept in spikes) ---
            [spkT, evTimes, useMaxDur, trialDur, kept_tbl] = makeZetaInputs_fromAOStartStopTimes( ...
                move_tbl, AO_spike_fs, ...
                'UseMaxDur_s', U.UseMaxDur_s, ...
                'PadITI_s',    U.PadITI_s, ...
                'SpikeField',  curSF, ...
                'StartField',  U.StartField, ...
                'StopField',   U.StopField, ...
                'PreWindow_s', U.PreWindow_s, ...
                'PostWindow_s',U.PostWindow_s);

            if isempty(spkT) || isempty(evTimes), continue; end
            nTrials = size(evTimes,1);

            % --- IFR (getIFR works on stitched timeline like zetatest) ---
            [vecTime, vecRate, sIFR] = getIFR( ...
                spkT, evTimes, useMaxDur, ...
                U.IFR_SmoothSd, U.IFR_MinScale, U.IFR_Base, U.IFR_UseParallel);

            % --- PSTH across trials (aligned to each event onset) ---
            bin_w   = U.BinSize_ms / 1000;              % seconds
            tmin    = -U.PreWindow_s;
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
                    ras_t   = [ras_t; rel];
                    ras_tr  = [ras_tr; i*ones(numel(rel),1)];
                    counts  = counts + histcounts(rel, edges);
                end
            end
        end


        % PSTH bins from -PreWindow_s to useMaxDur (s), normalized to Hz
        psth_Hz = (counts ./ nTrials) ./ bin_w;      % spikes/trial/bin  -> Hz

        % --- Simple figures (raster + IFR/PSTH overlay) ---
        hFig = [];
        if U.DoPlot
            hFig = figure('Color','w','Position',[100 100 900 600]);
            tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

            % Depth label
            depthMap = containers.Map({'t','c','b'},{'dorsal STN','central STN','ventral STN'});
            if isKey(depthMap,dz); depthLbl = depthMap(dz); else, depthLbl = dz; end

            % Raster
            ax1 = nexttile;
            if ~isempty(ras_t)
                scatter(ras_t, ras_tr, 6, 'k', 'filled'); hold on;
            end
            xline(0,'r-'); grid on;

            % y-axis: integer ticks only
            ylim([0.5 nTrials+0.5]);
            yticks(1:nTrials);
            ylabel('Trial Rep');

            % x-axis
            xlim([tmin tmax]);
            xlabel('Time from onset (s)');

            % Title format: Raster | MoveType - Depth (n=## reps)
            title(sprintf('Raster | %s - %s - %s (n=%d reps)', mv, depthLbl, curSF, nTrials));
            

            % IFR + PSTH overlay (two y-axes)
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
                relIFR_t = mod(vecTime, blockDur);
                keepIFR  = relIFR_t >= 0 & relIFR_t <= tmax; % [0,useMaxDur]
                plot(relIFR_t(keepIFR), vecRate(keepIFR), 'LineWidth', 2); hold on;
            end
            ylabel('IFR (Hz)');
            yl = ylim; ylim([0 max(yl(2), eps)]); % Force IFR y-axis to start at 0

            % Right y-axis: PSTH
            yyaxis right
            plot(centers, psth_Hz, '--', 'LineWidth', 2);
            ylabel(sprintf('PSTH (Hz, %d ms bins)', U.BinSize_ms));
            y2 = ylim; ylim([0 max(y2(2), eps)]);  % Force PSTH y-axis to start at 0

            % Shared x settings
            xlim([tmin tmax]); grid on;
            xlabel('Time from onset (s)');
            title('IFR (solid) + PSTH (dashed)');
        end

        % --- Collect outputs ---
        rows(end+1,:) = { ...
            curSF, mv, dz, nTrials, useMaxDur, ...
            U.BinSize_ms, mean(trialDur,'omitnan'), std(trialDur,'omitnan'), ...
            max(vecRate,[],'omitnan'), ...
            centers, psth_Hz, ...
            vecTime, vecRate};

        sOut = struct;
        sOut.SpikeField = curSF;
        sOut.MoveType = mv; sOut.Depth = dz; sOut.nTrials = nTrials;
        sOut.useMaxDur_s = useMaxDur;
        sOut.PreWindow_s = U.PreWindow_s; sOut.PostWindow_s = U.PostWindow_s;
        sOut.binSize_ms  = U.BinSize_ms;
        sOut.edges = edges; sOut.centers = centers; sOut.psth_Hz = psth_Hz;
        sOut.raster_t = ras_t; sOut.raster_trial = ras_tr;
        sOut.vecTime = vecTime; sOut.vecRate = vecRate; sOut.sIFR = sIFR;

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
    end
end

varNames = { ...
    'SpikeField','MoveType','Depth','nTrials','UseMaxDur_s', ...
    'BinSize_ms','MeanDur_s','StdDur_s', ...
    'MaxIFR_Hz','PSTH_TimeCenters_s','PSTH_Hz', ...
    'IFR_Time_s','IFR_Hz'};

if isempty(rows)
    IFR_PSTH_Summary = cell2table(cell(0,numel(varNames)), 'VariableNames', varNames);
else
    IFR_PSTH_Summary = cell2table(rows, 'VariableNames', varNames);
end
end

