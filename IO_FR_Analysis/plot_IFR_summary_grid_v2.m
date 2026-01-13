function hFig = plot_IFR_summary_grid_v2(IFR_PSTH_Summary, all_IFR, spikeFieldKey, depthOrder, moveOrder, varargin)
% plot_IFR_summary_grid
% Rows = depths (each depth has 2 rows: raster then IFR/PSTH)
% Cols = MoveTypes
% Saves tiled summary figure optionally.

% -----------------------
% Defaults / ordering
% -----------------------
spikeFieldKey = char(string(spikeFieldKey));

if nargin < 4 || isempty(depthOrder)
    depthOrder = {'t','c','b'};
end
if nargin < 5 || isempty(moveOrder)
    moveOrder  = {'HAND OC','HAND PS','ARM EF'};
end

% -----------------------
% Optional inputs (assign varargin)
% -----------------------
p = inputParser;
p.addParameter('SaveDir', '', @(s) ischar(s) || isstring(s));
p.addParameter('CaseDate','', @(s) ischar(s) || isstring(s));

p.addParameter('PadITI_s',    0.005, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('PreWindow_s', 0.050, @(x) isnumeric(x) && isscalar(x) && x>=0);

p.addParameter('XLim', [-0.05 1.60], @(x) isnumeric(x) && numel(x)==2); % interpret U.XLim(1) as -PreWindow_s

p.addParameter('AddGroupSummary', true, @islogical);
p.addParameter('MoveXMax', containers.Map({'HAND OC','HAND PS','ARM EF'},{0.6, 0.8, 1.4}));
p.addParameter('SmoothWin', 11, @(x) isnumeric(x) && isscalar(x) && x>=1);

p.addParameter('IFR_YLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
p.addParameter('IFR_YTicks', 0:10:50, @(x) isnumeric(x));
p.addParameter('Summary_YLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));

p.addParameter('AddLegends', true, @islogical)
p.addParameter('GroupMaster', table(), @(x) istable(x));

p.addParameter('GroupField_Spike', "SpikeField", @(x) isstring(x) || ischar(x));
p.addParameter('GroupField_Time',  "IFR_Time_s", @(x) isstring(x) || ischar(x));
p.addParameter('GroupField_Hz',    "IFR_Hz",     @(x) isstring(x) || ischar(x));

% group-summary timebase
p.addParameter('Group_DT', 0.002, @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0)); % e.g., 0.01 (10ms)
p.addParameter('Group_MaxN', 4000, @(x) isnumeric(x) && isscalar(x) && x>=100);

p.parse(varargin{:});
U = p.Results;

% Make the group summary robust if GroupMaster is empty
if U.AddGroupSummary && (isempty(U.GroupMaster) || height(U.GroupMaster)==0)
    warning('AddGroupSummary=true but GroupMaster not provided/empty. Disabling group summary row.');
    U.AddGroupSummary = false;
end


% -----------------------
% JNE color scheme
% -----------------------

% Depth colors (purple shades)
purpleShades = ([ ...
    118,42,131;      % dorsal STN (t)
    175,141,195;     % central STN (c)
    231-15,212-15,232-15] ... % ventral STN (b)
    ./ 255); % /255 = standard

% Movement-context colors (greens)
greenShades = ([ ...
    217-20,240-20,211-20;     % HAND OC  (light green)
    127,191,123;              % HAND PS  (green)
    27,120,55] ...            % ARM EF   (dark green)
    ./ 255);  % /255 = standard

depthColor = containers.Map({'t','c','b'}, {purpleShades(1,:), purpleShades(2,:), purpleShades(3,:)});
moveColor  = containers.Map({'HAND OC','HAND PS','ARM EF'}, {greenShades(1,:), greenShades(2,:), greenShades(3,:)});

depthName = containers.Map({'t','c','b'}, {'dorsal','central','ventral'});

nDepths = numel(depthOrder);
nMoveTs = numel(moveOrder);

% -----------------------
% Build key -> index map from all_IFR
% -----------------------
key2idx = containers.Map;
for i = 1:numel(all_IFR)
    s = all_IFR{i};
    sf = char(string(s.SpikeField));
    mv = char(string(s.MoveType));
    dz = char(string(s.Depth));
    key = sprintf('%s|%s|%s', sf, mv, dz);
    if ~isKey(key2idx, key)
        key2idx(key) = i;
    end
end

% -----------------------
% Figure layout
% -----------------------
hFig = figure('Color','w','Position',[50 50 1500 900]);

nExtraRows = 0;
if U.AddGroupSummary
    nExtraRows = 1; % one extra IFR-only summary row
end
tlo  = tiledlayout(2*nDepths + nExtraRows, nMoveTs, 'TileSpacing','compact','Padding','compact');

% Reserve space on the right for legend
tlo.Units = 'normalized';
tlo.Position = [0.05 0.06 0.80 0.90];  % leave ~20% width for legend


%% Main loop:

for di = 1:nDepths
    dz = char(string(depthOrder{di}));
    cDepth = [0 0 0];
    if isKey(depthColor, dz), cDepth = depthColor(dz); end

    dzLbl = dz;
    if isKey(depthName, dz), dzLbl = depthName(dz); end

    for mi = 1:nMoveTs
        mv = char(string(moveOrder{mi}));

        % Per-move visualization window
        xMax = U.XLim(2);
        if isKey(U.MoveXMax, mv)
            xMax = U.MoveXMax(mv);
        end
        thisXLim = [U.XLim(1) xMax];

        cMove = [0 0 0];
        if isKey(moveColor, mv), cMove = moveColor(mv); end

        key = sprintf('%s|%s|%s', spikeFieldKey, mv, dz);

        tileRaster = (2*(di-1))*nMoveTs + mi;       % rows 1,3,5...
        tileCurve  = (2*(di-1)+1)*nMoveTs + mi;     % rows 2,4,6...

        %% ----- Raster tile -----
        axR = nexttile(tlo, tileRaster);
        cla(axR); hold(axR,'on'); grid(axR,'on');

        if isKey(key2idx, key)
            sI = all_IFR{key2idx(key)};

            % Prefer precomputed raster_t0 over raster_t
            if isfield(sI,'raster_t0') && ~isempty(sI.raster_t0)
                tR = sI.raster_t0;
            elseif isfield(sI,'raster_t') && ~isempty(sI.raster_t)
                if isfield(sI,'PreWindow_s') && ~isempty(sI.PreWindow_s)
                    tR = sI.raster_t - U.PreWindow_s;
                else
                    tR = sI.raster_t; % fallback to original
                end
            end
            % Plot raster (points)
            if ~isempty(tR)
                scatter(axR, tR, sI.raster_trial, 6, cMove, 'filled', ... % update from sI.raster to sI.raster_t0
                    'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor','none');
            end
            xline(axR, 0, '-', 'Color',[0.3 0.3 0.3]); % true onset at 0

            % FIXED x-axis
            xlim(axR, thisXLim);

            % y-axis = repetition count (total trials/events)
            if isfield(sI,'nReps_behavioral') && ~isempty(sI.nReps_behavioral)
                nRep = double(sI.nReps_behavioral);   % rep-centric y-axis
            else
                nRep = double(sI.nTrials);            % fallback
            end

            ylim(axR, [0.5 nRep+0.5]);

            % y-tick spacing
            if nRep >= 25 % or 30
                yTickStep = 5;
            elseif nRep >= 10 % or 12
                yTickStep = 2;
            else
                yTickStep = 1;
            end
            yticks(axR, 1:yTickStep:nRep);
            % ylabel(axR, 'Repetition');

            % y-axis label should be Movement Number (or Trial Number)
            yLbl = 'Trial Number';

            if mi == 1
                ylabel(axR, yLbl);
            else
                ylabel(axR, '');
            end

        else
            axis(axR,'off');
        end

        if di == 1
            title(axR, mv, 'Interpreter','none');
        end
        if mi == 1
            % depth label once per depth block
            text(axR, thisXLim(1)-0.02*(thisXLim(2)-thisXLim(1)), 0.5, dzLbl, ...
                'Rotation', 90, 'HorizontalAlignment','right', 'VerticalAlignment','bottom', ...
                'FontWeight','bold');
        end


        %% ----- IFR + PSTH tile -----
        axC = nexttile(tlo, tileCurve);
        cla(axC); hold(axC,'on'); grid(axC,'on');

        if isKey(key2idx, key)
            sI = all_IFR{key2idx(key)};

            % % PSTH dashed (move color)
            % if isfield(sI,'centers') && ~isempty(sI.centers) && isfield(sI,'psth_Hz') && ~isempty(sI.psth_Hz)
            %     plot(axC, sI.centers, sI.psth_Hz, '--', 'LineWidth', 1.6, 'Color', cMove);
            % end

            % % PSTH: prefer centers_t0 if present
            % if isfield(sI,'centers_t0') && ~isempty(sI.centers_t0) && isfield(sI,'psth_Hz') && ~isempty(sI.psth_Hz)
            %     tC = sI.centers_t0;
            % elseif isfield(sI,'centers') && ~isempty(sI.centers) && isfield(sI,'psth_Hz') && ~isempty(sI.psth_Hz)
            %     if isfield(sI,'PreWindow_s') && ~isempty(sI.PreWindow_s)
            %         tC = sI.centers - sI.PreWindow_s;
            %     else
            %         tC = sI.centers;
            %     end
            % else
            %     tC = [];
            % end
            % % Plot PSTH (dashed)
            % if ~isempty(tC) && isfield(sI,'psth_Hz') && ~isempty(sI.psth_Hz)
            %     plot(axC, tC, sI.psth_Hz, '--', 'LineWidth', 1.6, 'Color', cMove);
            % end


            % IFR solid (depth color)
            if ~isempty(sI.vecTime) && ~isempty(sI.vecRate)
                t_ifr = sI.vecTime(:);
                r_ifr = sI.vecRate(:);

                % Smooth (movmean is stable and fast)
                r_sm = movmean(r_ifr, U.SmoothWin, 'omitnan');

                % % vecTime is stitched, per-category plots used modulo.
                % % For grid, keep it simple: plot what falls in fixed peri window.
                % mask  = (t_ifr >= thisXLim(1)) & (t_ifr <= thisXLim(2));
                % plot(axC, t_ifr(mask), r_sm(mask), 'LineWidth', 1.8, 'Color', cDepth);

                % Determine PreWindow + window length
                preW = 0;
                if isfield(sI,'PreWindow_s') && ~isempty(sI.PreWindow_s)
                    preW = double(sI.PreWindow_s);
                end

                % convert stitched vecTime -> peri-event t0
                % Use the actual category window (preferred) + PadITI
                useMaxDur = [];
                if isfield(sI,'useMaxDur_s') && ~isempty(sI.useMaxDur_s)
                    useMaxDur = double(sI.useMaxDur_s);
                else
                    useMaxDur = U.XLim(2); % fallback
                end

                blockDur = useMaxDur + U.PadITI_s;

                % peri-event time where 0 = padded start of each block
                rel_pad0 = mod(t_ifr, blockDur);

                % shift so 0 = true onset (i.e., after PreWindow_s)
                t0 = rel_pad0 - preW;

                % only plot within this tile's window
                keep = (t0 >= thisXLim(1)) & (t0 <= thisXLim(2)) & isfinite(r_sm);
                plot(axC, t0(keep), r_sm(keep), 'LineWidth', 1.8, 'Color', cDepth);
            end

            xline(axC, 0, '-', 'Color',[0.3 0.3 0.3]);
            xlim(axC, thisXLim);

            % Add y-label "Hz"
            if mi == 1
                ylabel(axC, 'Hz');
            else
                ylabel(axC, '');
            end

            % yticks(axC, []); % work on yticks next
            % y ticks + optional ylim
            yticks(axC, U.IFR_YTicks);
            if ~isempty(U.IFR_YLim)
                ylim(axC, U.IFR_YLim);
            end


        else
            axis(axC,'off');
        end

        if di == nDepths
            xlabel(axC,'Time from onset (s)');
        end
    end
end


%% Bottom-row: across-subject mean IFR per MoveType (3 traces: t/c/b)

if U.AddGroupSummary
    for mi = 1:nMoveTs
        mv = char(string(moveOrder{mi}));

        % Per-move visualization window
        xMax = U.XLim(2); % movetype-specific useMaxDur
        if isKey(U.MoveXMax, mv), xMax = U.MoveXMax(mv); end
        thisXLim = [U.XLim(1) xMax];

        tileSummary = (2*nDepths)*nMoveTs + mi; % first tile after the main grid
        axS = nexttile(tlo, tileSummary);
        cla(axS); hold(axS,'on'); grid(axS,'on');

        [tC, M] = computeGroupIFR_from_MasterZETA( ...
            U.GroupMaster, spikeFieldKey, mv, depthOrder, thisXLim, U.SmoothWin, ...
            string(U.GroupField_Spike), string(U.GroupField_Time), string(U.GroupField_Hz), ...
            U.Group_DT, U.Group_MaxN, ...
            U.PreWindow_s, U.PadITI_s, xMax);

        % Plot in depth colors
        for di2 = 1:nDepths
            dz2 = char(string(depthOrder{di2}));
            cDepth = [0.4 0.4 0.4];
            if isKey(depthColor, dz2), cDepth = depthColor(dz2); end
            plot(axS, tC, M(di2,:), 'LineWidth', 2.2, 'Color', cDepth);
        end

        xline(axS, 0, '-', 'Color',[0.3 0.3 0.3]);
        xlim(axS, thisXLim);

        if mi == 1
            ylabel(axS, 'Hz');
        end
        xlabel(axS, 'Time from onset (s)');

        yticks(axS, U.IFR_YTicks);
        if ~isempty(U.Summary_YLim)
            ylim(axS, U.Summary_YLim);
        elseif ~isempty(U.IFR_YLim)
            ylim(axS, U.IFR_YLim);
        end

        title(axS, sprintf('%s | Group mean', mv), 'Interpreter','none');
    end
end

title(tlo, sprintf('SpikeField %s: Raster + IFR by Depth × MoveType', spikeFieldKey), ...
    'Interpreter','none');


%% Add legend in top-right margin outside the tiles

if U.AddLegends
    axL = axes(hFig, 'Units','normalized', 'Position',[0.87 0.80 0.12 0.18], 'Visible','off');
    hold(axL,'on');

    % Depth legend (lines)
    hD = gobjects(1, nDepths);
    depthLabels = strings(1,nDepths);
    for di2 = 1:nDepths
        dz2 = char(string(depthOrder{di2}));
        cD = depthColor(dz2);
        hD(di2) = plot(axL, nan, nan, '-', 'LineWidth', 2.5, 'Color', cD);

        if isKey(depthName, dz2)
            depthLabels(di2) = string(depthName(dz2));
        else
            depthLabels(di2) = string(dz2);
        end
    end

    % MoveType legend (dots)
    hM = gobjects(1, nMoveTs);
    moveLabels = strings(1,nMoveTs);
    for mi2 = 1:nMoveTs
        mv2 = char(string(moveOrder{mi2}));
        cM = moveColor(mv2);
        hM(mi2) = scatter(axL, nan, nan, 30, cM, 'filled');
        moveLabels(mi2) = string(mv2);
    end

    lgd = legend(axL, [hD(:); hM(:)], [depthLabels(:); moveLabels(:)], ...
        'Location','northwest');
    lgd.Box = 'off';
end


%% Save alongside per-category IFR plots

if ~isempty(U.SaveDir)
    if ~exist(U.SaveDir,'dir'), mkdir(U.SaveDir); end

    if ~isempty(U.CaseDate)
        base = sprintf('%s_IFRGrid_%s', char(string(U.CaseDate)), spikeFieldKey);
    else
        base = sprintf('IFRGrid_%s', spikeFieldKey);
    end
    base = strrep(base,' ','_');

    print(hFig, fullfile(U.SaveDir, [base,'.png']), '-dpng','-r300');
    savefig(hFig, fullfile(U.SaveDir, [base,'.fig']));
end


%% Local helper: Bottom-row across-subject mean IFR (3 traces per MoveType)

    function [tC, meanByDepth] = computeGroupIFR_from_MasterZETA( ...
            MasterT, spikeFieldKey, mv, depthOrder, thisXLim, smoothWin, ...
            colSpike, colTime, colHz, ...
            forceDT, maxN, ...
            preW_default, padITI_s, useMaxDur_default)

        % Expects MasterT has columns:
        % SpikeField, MoveType, Depth, IFR_Time_s, IFR_Hz
        % where IFR_Time_s and IFR_Hz are vectors (often stored as cells).

        xMin = thisXLim(1);
        xMax = thisXLim(2);

        if isempty(MasterT) || ~istable(MasterT)
            error('GroupMaster must be a non-empty table (MasterZETA or MasterZETA_MUA).');
        end

        colSpike = string(colSpike);
        colTime  = string(colTime);
        colHz    = string(colHz);

        need = [string("MoveType"), string("Depth"), colSpike, colTime, colHz];
        missing = setdiff(need, string(MasterT.Properties.VariableNames));
        if ~isempty(missing)
            error('GroupMaster missing required columns: %s', strjoin(cellstr(missing), ', '));
        end

        % Masks for this spikeField + MoveType subset (used for dt inference + row selection)
        maskSF = strcmpi(string(strtrim(MasterT.(colSpike))), string(spikeFieldKey));
        maskMV = strcmpi(string(strtrim(MasterT.MoveType)),   string(mv));


        % ---- Build a common timebase tC ----
        xMin = thisXLim(1);
        xMax = thisXLim(2);

        % If user forces dt, use it
        dt = [];
        if ~isempty(forceDT)
            dt = double(forceDT);
        else
            % Try infer dt from the first valid trace
            idx0   = find(maskSF & maskMV, 1, 'first');

            if ~isempty(idx0)
                [tv, ~] = unpackTrace(MasterT.(colTime)(idx0), MasterT.(colHz)(idx0));

                % Only infer dt if time looks sane
                tv = tv(isfinite(tv));
                if numel(tv) >= 10
                    d = diff(tv);
                    d = d(isfinite(d) & d>0);
                    if ~isempty(d)
                        dtCand = median(d);

                        % Guardrails: reject absurd dt (too tiny -> huge arrays)
                        if isfinite(dtCand) && dtCand >= 1e-4 && dtCand <= 0.5
                            dt = dtCand;
                        end
                    end
                end
            end
        end

        % Fallback dt if inference failed
        if isempty(dt) || ~isfinite(dt) || dt <= 0
            dt = 0.01; % 10 ms default
        end

        % Cap total points to avoid memory blowups
        span = xMax - xMin;
        N_est = floor(span/dt) + 1;

        if N_est > maxN
            dt = span / (maxN-1);
        end

        tC = linspace(xMin, xMax, floor(span/dt)+1);
        tC = tC(:)';  % row vector

        nT = numel(tC);

        meanByDepth = nan(numel(depthOrder), nT);

        % Compute mean trace per depth
        for di2 = 1:numel(depthOrder)
            dz2 = char(string(depthOrder{di2}));

            maskDZ = strcmpi(string(strtrim(MasterT.Depth)), string(dz2));
            rows = find(maskSF & maskMV & maskDZ);

            Y = nan(numel(rows), nT);
            keepRow = false(numel(rows),1);

            for rr = 1:numel(rows)
                r = rows(rr);

                [tv, yv] = unpackTrace(MasterT.(colTime)(r), MasterT.(colHz)(r));
                if numel(tv) < 6 || numel(yv) < 6, continue; end

                % Smooth per-trace (movmean is stable)
                yv = movmean(yv, smoothWin, 'omitnan');

                % % Restrict to visualization window
                % keep = (tv >= xMin) & (tv <= xMax) & isfinite(tv) & isfinite(yv);
                % tv2 = tv(keep);
                % yv2 = yv(keep);
                %
                % if numel(tv2) < 6, continue; end
                %
                % % Interpolate to common grid
                % yi = interp1(tv2, yv2, tC, 'linear', NaN);

                % --- Convert stitched time -> peri-event time (0 = true onset) ---
                preW = double(preW_default);
                useMaxDur = double(useMaxDur_default);

                % Optional: if GroupMaster contains per-row PreWindow_s / UseMaxDur_s, prefer those
                if ismember("PreWindow_s", string(MasterT.Properties.VariableNames))
                    try
                        preW = double(MasterT.PreWindow_s(r));
                    catch
                    end
                end
                if ismember("UseMaxDur_s", string(MasterT.Properties.VariableNames))
                    try
                        useMaxDur = double(MasterT.UseMaxDur_s(r));
                    catch
                    end
                end

                blockDur = useMaxDur + double(padITI_s);

                % peri-event time where 0 = padded start of each stitched block
                rel_pad0 = mod(tv, blockDur);

                % shift so that 0 = true onset (after PreWindow)
                t0 = rel_pad0 - preW;

                % Restrict to visualization window in t0 coordinates
                keep = (t0 >= xMin) & (t0 <= xMax) & isfinite(t0) & isfinite(yv);
                t02 = t0(keep);
                yv0 = yv(keep);

                if numel(t02) < 6, continue; end

                % --- Bin onto common grid instead of interp1 (robust + smooth) ---
                dtC = median(diff(tC));                 % bin width from common grid
                edges = [tC - dtC/2, tC(end) + dtC/2];  % nT+1 edges

                % % Convert times to bin indices in [1..nT]
                % binIdx = round((t0 - tC(1)) / dtC) + 1;
                % inRange = (binIdx >= 1) & (binIdx <= nT);
                % binIdx = binIdx(inRange);
                % y0 = yv0(inRange);
                binIdx = discretize(t02, edges);        % size == numel(t02)
                valid  = ~isnan(binIdx) & isfinite(yv0);

                if numel(binIdx) < 6
                    continue
                end

                % % nnz: Number of nonzero matrix elements
                % if nnz(valid) < 6
                %     continue
                % end

                % mean in each bin (averages duplicates)
                yb = accumarray(binIdx(valid), yv0(valid), [nT 1], @mean, NaN);
                yi = yb.';  % row vector 1 x nT


                % % --- Interp1 requires unique sample points (t02). After mod(), duplicates are expected.
                % % Quantize to dt grid to merge near-duplicates
                % t02q = round(t02 / dt) * dt;
                %
                % % Aggregate duplicates by averaging yv2 at identical t02 values.
                % [t02s, ord] = sort(t02q(:));
                % yv2s = yv2(:);
                % yv2s = yv2s(ord);
                %
                % [tU, ~, ic] = unique(t02s, 'stable');
                % yU = accumarray(ic, yv2s, [], @mean);
                %
                % % Guard: need enough unique points
                % if numel(tU) < 6
                %     continue;
                % end
                %
                % % Interpolate to common grid (tC already in shifted coordinates)
                % yi = interp1(tU, yU, tC, 'linear', NaN);
                

                Y(rr,:) = yi;
                keepRow(rr) = true;
            end

            if any(keepRow)
                meanByDepth(di2,:) = mean(Y(keepRow,:), 1, 'omitnan');
                meanByDepth(di2,:) = movmean(meanByDepth(di2,:), smoothWin, 'omitnan');
            end
        end

    end


    function [tv, yv] = unpackTrace(tCellOrVec, yCellOrVec)
        % Handles cells, numeric arrays, and strings (rare but defensive)

        tv = []; yv = [];

        % Time
        if iscell(tCellOrVec), tCellOrVec = tCellOrVec{1}; end
        if isstring(tCellOrVec) || ischar(tCellOrVec), tCellOrVec = str2num(char(tCellOrVec)); end
        if isnumeric(tCellOrVec), tv = double(tCellOrVec(:)); end

        % IFR
        if iscell(yCellOrVec), yCellOrVec = yCellOrVec{1}; end
        if isstring(yCellOrVec) || ischar(yCellOrVec), yCellOrVec = str2num(char(yCellOrVec)); end
        if isnumeric(yCellOrVec), yv = double(yCellOrVec(:)); end

        % Basic cleanup
        if isempty(tv) || isempty(yv), tv = []; yv = []; return; end
        n = min(numel(tv), numel(yv));
        tv = tv(1:n);
        yv = yv(1:n);

        % Ensure increasing time
        [tv, ord] = sort(tv);
        yv = yv(ord);

    end

end



