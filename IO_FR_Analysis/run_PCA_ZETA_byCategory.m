function run_PCA_ZETA_byCategory(MasterZETA, varargin)

% run_PCA_ZETA_byCategory
%
% Do PCA on ZETA temporal deviation vectors (ZETA_vecD) for each
% MoveType × STN depth, after resampling onto a common time axis given by
% PSTH_TimeCenters_s.

% 1) Navigate to MasterZeta file location
% Aggr_ZETA_dir = [FR_Kin_Dir, filesep, 'Aggregate Zeta Plots'];
% cd(Aggr_ZETA_dir)

% 2) Load master ZETA file (produced by aggregate_ZETA_and_plot): MasterZETA_AllData.mat
% MasterZETA = readtable('MasterZeta_AllSubjects.csv');
% MasterZETA.Properties.VariableNames;

% 3) Extract ZETA values from MasterZeta file for each MoveType × STN depth
% across subject spike clusters
% Build a matrix: X=[timeBins × units]

% 4) Run PCA on ZETA temporal deviation vectors (Zeta_vecD) per category
% % Category: each MoveType per STN depth (E.g., Hand OC x dorsal STN)
% x-axis (t, time): PSTH_TimeCenters_s (or Zeta_vecT), time vector (10 ms bins), uniform across subjects
% y-axis (d(t), features): ZETA_vecD, temporal deviation signal d(t) that ZETA is computed from
% % One vector (signals per unit) per MoveType × depth, defined over time support ZETA_vecT
%
% Then plot PC1 per MoveType for each Depth (3 rows: t/c/b).

p = inputParser;
p.addParameter('SavePath','', @(x)ischar(x)||isstring(x));
p.addParameter('FR_Kin_Dir','', @(x)ischar(x)||isstring(x)); % optional helper
% PC1 peak per MoveType
p.addParameter('AddPC1PeakLines', true, @islogical);
p.addParameter('PeakWindow_s', [0.0 1.6], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('PeakMode', 'abs', @(s) any(strcmpi(string(s), ["max","abs"]))); % max(PC1) or max(abs(PC1))
p.addParameter('PeakLineWidth', 1.5, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('PeakLineStyle', '--', @(s) ischar(s) || isstring(s));
% smoothing-based peak finding
p.addParameter('PeakSmooth_N', 9, @(x) isnumeric(x) && isscalar(x) && x>=1); % odd recommended, PeakSmooth_N=9 is ~90 ms if dt=0.01. 
p.addParameter('PeakSmooth_Method', 'movmean', @(s) any(strcmpi(string(s), ["movmean","movmedian"])));

p.parse(varargin{:});
U = p.Results;


% Default SavePath if not provided
if isempty(U.SavePath)
    if ~isempty(U.FR_Kin_Dir)
        U.SavePath = fullfile(char(U.FR_Kin_Dir), 'Aggregate Zeta Plots', 'PCA Plots');
    else
        U.SavePath = fullfile(pwd, 'PCA Plots');
    end
end
if ~exist(U.SavePath,'dir'), mkdir(U.SavePath); end


fprintf('[OK] Loaded MasterZETA_AllData: %d total entries\n', height(MasterZETA));

%% --- Sanity check required columns ---

% Adjust these names if your table uses Zeta_vecD/Zeta_vecT instead:
reqVars = {'MoveType','Depth','ZETA_vecD','ZETA_vecT', ...
    'PSTH_TimeCenters_s','SpikeField','PrettyLabel'};
missing = setdiff(reqVars, MasterZETA.Properties.VariableNames);
if ~isempty(missing)
    error('MasterZETA is missing required variables: %s', strjoin(missing,', '));
end

%% Restrict to rows that have both ZETA vecD and PSTH time centers

hasZ = ~cellfun(@isempty, MasterZETA.ZETA_vecD);
hasT = ~cellfun(@isempty, MasterZETA.PSTH_TimeCenters_s);
M    = MasterZETA(hasZ & hasT, :);

if isempty(M)
    error('No rows with both ZETA_vecD and PSTH_TimeCenters_s found.');
end

fprintf('[OK] Unique MoveTypes: %s\n', strjoin(unique(M.MoveType,'stable'), ', '));

%% MoveType and STN Depth mapping

moveTypes  = unique(M.MoveType,'stable');
depths     = {'t','c','b'};  % dorsal, central, ventral
depthNames = containers.Map({'t','c','b'}, ...
    {'dorsal STN','central STN','ventral STN'});

% % Colors by MoveType
% mtColor = containers.Map( ...
%     {'HAND OC','HAND PS','ARM EF','REST'}, ...
%     {[0.95 0.60 0.10], ... % Hand OC = orange
%      [0.20 0.65 0.30], ... % Hand PS = green
%      [0.15 0.45 0.85], ... % Arm EF  = blue
%      [0.60 0.60 0.60]});   % Rest    = gray

% -----------------------
% JNE color scheme (canonical)
% -----------------------
purpleShades = ([ ...
    118,42,131;      % dorsal STN (t)
    175,141,195;     % central STN (c)
    231-15, 212-15, 232-15] ... % ventral STN (b)
    ./ 255);

greenShades = ([ ...
    128,128,128;     % REST
    217-15,240-15,211-15;     % HAND OC
    127,191,123;     % HAND PS
    27,120,55] ...   % ARM EF
    ./ 255);

depthColorMap = containers.Map({'t','c','b'}, {purpleShades(1,:), purpleShades(2,:), purpleShades(3,:)});
moveColorMap  = containers.Map({'REST','HAND OC','HAND PS','ARM EF'}, {greenShades(1,:), greenShades(2,:), greenShades(3,:), greenShades(4,:)});
fallbackCol   = [0.5 0.5 0.5];


%% Storage for PC1 waveforms per MoveType × Depth

PC1 = struct();
for d = 1:numel(depths)
    for m = 1:numel(moveTypes)
        key = makeKey(depths{d}, moveTypes{m});
        PC1.(key).t       = [];
        PC1.(key).pc1     = [];
        PC1.(key).nUnits  = 0;
        PC1.(key).latent  = [];
        PC1.(key).scores  = [];
        PC1.(key).labels  = [];
    end
end

%% ----- MAIN LOOP: per MoveType × Depth -----

for d = 1:numel(depths)
    dz = depths{d};

    for m = 1:numel(moveTypes)
        mv = moveTypes{m};

        % Subset rows in this category
        catRows = M(M.MoveType==mv & M.Depth==dz, :);
        if height(catRows) < 2
            % Not enough units to do PCA
            continue;
        end

        % Reference time axis from PSTH_TimeCenters_s (first valid row)
        time_ref = catRows.PSTH_TimeCenters_s{1};
        if isempty(time_ref) || ~isnumeric(time_ref)
            warning('Empty or non-numeric PSTH_TimeCenters_s for %s × %s; skipping.', mv, dz);
            continue;
        end
        time_ref = time_ref(:);  % column vector

        nUnits = height(catRows);
        X       = nan(numel(time_ref), nUnits);  % time × units
        keepUnit = false(1, nUnits);

        for u = 1:nUnits
            dVec = catRows.ZETA_vecD{u};
            tVec   = catRows.ZETA_vecT{u};

            if isempty(dVec) || isempty(tVec), continue; end
            dVec = double(dVec(:));
            tVec   = double(tVec(:));

            if numel(dVec) < 3 || numel(tVec) ~= numel(dVec)
                continue;
            end

            % Interpolate ZETA deviation onto common PSTH time base
            try
                X(:,u) = interp1(tVec, dVec, time_ref, 'linear', 'extrap');
                keepUnit(u) = true;
            catch ME
                warning('interp1 failed for unit %d in %s × %s: %s', ...
                    u, mv, dz, ME.message);
            end
        end

        % Keep only columns with actual data
        X = X(:, keepUnit);
        if size(X,2) < 2
            % need at least 2 units to meaningfully do PCA
            continue;
        end

        % Remove timepoints that are NaN across all units
        goodTime = any(isfinite(X), 2);
        X    = X(goodTime, :);
        t_use = time_ref(goodTime);

        % Replace remaining NaNs (per unit) with unit-wise mean over non-NaN times
        for j = 1:size(X,2)
            nanMask = ~isfinite(X(:,j));
            if any(nanMask)
                mu_j = mean(X(~nanMask,j), 'omitnan');
                X(nanMask,j) = mu_j;
            end
        end

        % ----- PCA -----
        % Observations = units, variables = timepoints → X'
        [coeff, score, latent] = pca(X', 'Centered', true);

        % coeff: [nTime × nPC], columns are PCs across time
        pc1 = coeff(:,1);  % time course of PC1

        % Flip sign so that the mean of PC1 in the post-onset window is positive
        postMask = t_use >= 0 & t_use <= (max(t_use) * 0.5);
        if any(postMask)
            if mean(pc1(postMask)) < 0
                pc1 = -pc1;
                score(:,1) = -score(:,1);
            end
        end

        key = makeKey(dz, mv);
        PC1.(key).t       = t_use;
        PC1.(key).pc1     = pc1;
        PC1.(key).nUnits  = size(X,2);
        PC1.(key).latent  = latent;
        PC1.(key).scores  = score;
        PC1.(key).labels  = catRows.PrettyLabel(keepUnit);

        fprintf('[PCA] %s × %s: %d units, %d time points\n', ...
            mv, dz, size(X,2), numel(t_use));
    end
end

%% ----- PLOT: PC1 per depth (rows) and MoveType (colors) -----

% Consistent MoveType plotting order
mtOrder = {'HAND OC','HAND PS','ARM EF','REST'};
mtOrder = intersect(mtOrder, moveTypes, 'stable');

% Consistent MoveType plotting order
mtOrder = {'HAND OC','HAND PS','ARM EF','REST'};
mtOrder = intersect(mtOrder, moveTypes, 'stable');

% ---- Compute global y-limits for PC1 across all depths × MoveTypes ----
allPC1 = [];
for d = 1:numel(depths)
    dz = depths{d};
    for m = 1:numel(mtOrder)
        mv  = mtOrder{m};
        key = makeKey(dz, mv);
        if isfield(PC1, key) && ~isempty(PC1.(key).pc1)
            allPC1 = [allPC1; PC1.(key).pc1(:)];
        end
    end
end

if ~isempty(allPC1)
    yPad  = 0.015;              % offset
    yLims = [min(allPC1)-yPad, max(allPC1)+yPad];
else
    yLims = [];
end
% -----------------------------------------------------------------------

hFig = figure('Color','w','Position',[100 100 950 800]);
tlo = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

for d = 1:numel(depths)
    dz = depths{d};
    ax = nexttile; hold(ax,'on'); grid(ax,'on');

    % Store traces plotted in this axis so we can compute peaks afterward
    plotted = struct('mv', {}, 't', {}, 'pc1', {}, 'c', {});

    for m = 1:numel(mtOrder)
        mv  = mtOrder{m};
        key = makeKey(dz, mv);

        if ~isfield(PC1, key), continue; end
        if isempty(PC1.(key).t), continue; end

        t_use = PC1.(key).t;
        pc1   = PC1.(key).pc1;

        if isKey(moveColorMap, mv)
            c = moveColorMap(mv);
        else
            c = fallbackCol;
        end

        plot(ax, t_use, pc1, 'LineWidth', 2, 'Color', c, 'DisplayName', mv);

        plotted(end+1) = struct('mv', mv, 't', t_use, 'pc1', pc1, 'c', c); %#ok<AGROW>
    end

    % Onset line
    xline(ax, 0, 'k--', 'LineWidth', 1);

    % ---- PC1 peak lines per MoveType ----
    if U.AddPC1PeakLines
        for k = 1:numel(plotted)
            mv  = plotted(k).mv;
            t   = plotted(k).t;
            y   = plotted(k).pc1;
            col = plotted(k).c;

            % optional: skip REST to reduce clutter
            if strcmpi(mv,'REST'), continue; end

            % call PC1 peak-finding helper
            % tPk = local_getPC1PeakTime(t, y, U.PeakWindow_s, U.PeakMode);
            tPk = local_getPC1PeakTime(t, y, U.PeakWindow_s, U.PeakMode, ...
                           U.PeakSmooth_N, U.PeakSmooth_Method);


            if ~isnan(tPk)
                xline(ax, tPk, U.PeakLineStyle, ...
                    'Color', col, 'LineWidth', U.PeakLineWidth, ...
                    'HandleVisibility','off'); % keep legend clean
            end
        end
    end
    % --------------------------------------

    xlabel(ax, 'Time from movement onset (s)');
    ylabel(ax, 'PC1 (a.u.)');

    if isKey(depthNames, dz)
        title(ax, sprintf('PC1 of ZETA deviation | %s', depthNames(dz)));
    else
        title(ax, sprintf('PC1 of ZETA deviation | depth %s', dz));
    end

    % xlim(ax, [-0.2 1.6]);  % match SU/MUA PC1 plots

    if ~isempty(yLims)
        ylim(ax, yLims);
    end
end

% Global legend
legend(tlo.Children(end), mtOrder, 'Location','northeastoutside');          % eastoutside; % change location to top right corner (in line with title)
title(tlo, 'PC1 of ZETA temporal deviation per MoveType × STN depth');


%% Save figure

outPNG = fullfile(U.SavePath, 'ZETA_PC1_SU_byDepth.png');  % in SU function
exportgraphics(hFig, outPNG, 'Resolution',300);

outFIG = fullfile(U.SavePath, 'ZETA_PC1_SU_byDepth.fig');
savefig(hFig, outFIG);


%% ----- PLOT (No REST): PC1 per depth (rows) and MoveType (colors) -----

mtOrder_NR = {'HAND OC','HAND PS','ARM EF'};
mtOrder_NR = intersect(mtOrder_NR, moveTypes, 'stable');

% y-lims for No REST
allPC1_NR = [];
for d = 1:numel(depths)
    dz = depths{d};
    for m = 1:numel(mtOrder_NR)
        mv  = mtOrder_NR{m};
        key = makeKey(dz, mv);
        if isfield(PC1, key) && ~isempty(PC1.(key).pc1)
            allPC1_NR = [allPC1_NR; PC1.(key).pc1(:)];
        end
    end
end

if ~isempty(allPC1_NR)
    yPadNR  = 0.015;
    yLimsNR = [min(allPC1_NR)-yPadNR, max(allPC1_NR)+yPadNR];
else
    yLimsNR = [];
end

hFigNR = figure('Color','w','Position',[100 100 950 800]);
tloNR  = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

for d = 1:numel(depths)
    dz = depths{d};
    ax = nexttile; hold(ax,'on'); grid(ax,'on');

    plotted = struct('mv', {}, 't', {}, 'pc1', {}, 'c', {});

    for m = 1:numel(mtOrder_NR)
        mv  = mtOrder_NR{m};
        key = makeKey(dz, mv);

        if ~isfield(PC1, key) || isempty(PC1.(key).t), continue; end

        t_use = PC1.(key).t;
        pc1   = PC1.(key).pc1;

        if isKey(moveColorMap, mv), c = moveColorMap(mv); else, c = fallbackCol; end

        plot(ax, t_use, pc1, 'LineWidth', 2, 'Color', c, 'DisplayName', mv);
        plotted(end+1) = struct('mv', mv, 't', t_use, 'pc1', pc1, 'c', c); %#ok<AGROW>
    end

    xline(ax, 0, 'k--', 'LineWidth', 1);

    if U.AddPC1PeakLines
        for k = 1:numel(plotted)
            tPk = local_getPC1PeakTime(plotted(k).t, plotted(k).pc1, ...
                U.PeakWindow_s, U.PeakMode, U.PeakSmooth_N, U.PeakSmooth_Method);

            if ~isnan(tPk)
                xline(ax, tPk, U.PeakLineStyle, 'Color', plotted(k).c, ...
                    'LineWidth', U.PeakLineWidth, 'HandleVisibility','off');
            end
        end
    end

    xlabel(ax, 'Time from movement onset (s)');
    ylabel(ax, 'PC1 (a.u.)');
    if isKey(depthNames, dz)
        title(ax, sprintf('PC1 of ZETA deviation | %s', depthNames(dz)));
    else
        title(ax, sprintf('PC1 of ZETA deviation | depth %s', dz));
    end

    if ~isempty(yLimsNR), ylim(ax, yLimsNR); end
end

legend(tloNR.Children(end), mtOrder_NR, 'Location','northeastoutside');
title(tloNR, 'PC1 of ZETA temporal deviation per MoveType × STN depth');

%% Save figure (No REST)

outPNG_NR = fullfile(U.SavePath, 'ZETA_PC1_SU_byDepth_NoREST.png');
exportgraphics(hFigNR, outPNG_NR, 'Resolution',300);

outFIG_NR = fullfile(U.SavePath, 'ZETA_PC1_SU_byDepth_NoREST.fig');
savefig(hFigNR, outFIG_NR);


end


%% --------- Helper: build a valid struct field name ---------

function key = makeKey(dz, mv)
% dz: depth string, e.g. 't'
% mv: MoveType string, e.g. 'HAND OC'
mv = char(mv);
mv_clean = regexprep(mv, '\W', '');  % remove non-alphanumeric chars
key = sprintf('%s_%s', dz, mv_clean); % e.g. 't_HANDOC'
end


%% Helper: PC1 peak time index per MoveType

function tPk = local_getPC1PeakTime(t, y, win_s, modeStr, smoothN, smoothMethod)
    tPk = NaN;
    if isempty(t) || isempty(y), return; end

    t = double(t(:));
    y = double(y(:));

    w0 = min(win_s); w1 = max(win_s);
    mask = isfinite(t) & isfinite(y) & (t >= w0) & (t <= w1);
    if nnz(mask) < 3, return; end

    tt = t(mask);
    yy = y(mask);

    % --- smooth only for peak finding ---
    smoothN = max(1, round(smoothN));
    if strcmpi(smoothMethod,'movmedian')
        yy_s = movmedian(yy, smoothN, 'omitnan');
    else
        yy_s = movmean(yy, smoothN, 'omitnan');
    end

    if strcmpi(modeStr,'abs')
        [~, iPk] = max(abs(yy_s));
    else
        [~, iPk] = max(yy_s);
    end

    tPk = tt(iPk);
end

