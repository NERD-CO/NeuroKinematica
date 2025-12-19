function run_PCA_ZETA_MUA_byCategory(MasterZETA_MUA, varargin)

% run_PCA_ZETA_MUA_byCategory
%
% PCA on MUA ZETA temporal deviation vectors (vecD_MUA) for each
% MoveType × STN depth, using new_vecTime_MUA as the time axis.
%
% Each row (subject/hemisphere, MoveType, Depth) = one "unit".
% Plots PC1 per MoveType for each Depth (3 rows: t/c/b).

p = inputParser;
p.addParameter('SavePath','', @(x)ischar(x)||isstring(x));
p.addParameter('FR_Kin_Dir','', @(x)ischar(x)||isstring(x)); % optional helper
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


fprintf('[OK] MasterZETA_MUA: %d total entries\n', height(MasterZETA_MUA));

% Required columns
reqVars = {'MoveType','Depth','vecD_MUA','new_vecTime_MUA','PrettyLabel'};
missing = setdiff(reqVars, MasterZETA_MUA.Properties.VariableNames);
if ~isempty(missing)
    error('MasterZETA_MUA is missing required variables: %s', strjoin(missing,', '));
end

% Restrict to rows that actually have ZETA MUA vectors
hasD = ~cellfun(@isempty, MasterZETA_MUA.vecD_MUA);
hasT = ~cellfun(@isempty, MasterZETA_MUA.new_vecTime_MUA);
M = MasterZETA_MUA(hasD & hasT, :);
if isempty(M)
    error('No rows with non-empty vecD_MUA and new_vecTime_MUA.');
end

% Required columns
reqVars = {'MoveType','Depth','vecD_MUA','new_vecTime_MUA','PrettyLabel'};
missing = setdiff(reqVars, MasterZETA_MUA.Properties.VariableNames);
if ~isempty(missing)
    error('MasterZETA_MUA is missing required variables: %s', strjoin(missing,', '));
end

% Keep only rows with actual ZETA TS vectors
hasD = ~cellfun(@isempty, MasterZETA_MUA.vecD_MUA);
hasT = ~cellfun(@isempty, MasterZETA_MUA.new_vecTime_MUA);
M = MasterZETA_MUA(hasD & hasT, :);
if isempty(M)
    error('No rows with non-empty vecD_MUA and new_vecTime_MUA.');
end

%% MoveType and STN Depth mapping

moveTypes  = unique(M.MoveType,'stable');
depths     = {'t','c','b'};
depthNames = containers.Map({'t','c','b'}, ...
    {'dorsal STN','central STN','ventral STN'});

% Colors by MoveType (same as SU)
% mtColor = containers.Map( ...
%     {'HAND OC','HAND PS','ARM EF','REST'}, ...
%     {[0.95 0.60 0.10], ... % Hand OC = orange
%     [0.20 0.65 0.30], ... % Hand PS = green
%     [0.15 0.45 0.85], ... % Arm EF  = blue
%     [0.60 0.60 0.60]});   % Rest    = gray

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
    217,240,211;     % HAND OC
    127,191,123;     % HAND PS
    27,120,55] ...   % ARM EF
    ./ 255);

depthColorMap = containers.Map({'t','c','b'}, {purpleShades(1,:), purpleShades(2,:), purpleShades(3,:)});
moveColorMap  = containers.Map({'REST','HAND OC','HAND PS','ARM EF'}, {greenShades(1,:), greenShades(2,:), greenShades(3,:), greenShades(4,:)});
fallbackCol   = [0.5 0.5 0.5];


%% Storage for PC1 per category

PC1_MUA = struct();
for d = 1:numel(depths)
    for m = 1:numel(moveTypes)
        rawKey = sprintf('%s_%s', depths{d}, moveTypes{m});
        key = matlab.lang.makeValidName(rawKey);

        PC1_MUA.(key).t      = [];
        PC1_MUA.(key).pc1    = [];
        PC1_MUA.(key).nUnits = 0;
        PC1_MUA.(key).latent = [];
        PC1_MUA.(key).scores = [];
        PC1_MUA.(key).labels = [];
    end
end


%% -------- MAIN LOOP: per MoveType × Depth --------
for d = 1:numel(depths)
    dz = depths{d};
    for m = 1:numel(moveTypes)
        mv = moveTypes{m};

        catRows = M(M.MoveType==mv & M.Depth==dz, :);
        if height(catRows) < 2
            continue;   % need >=2 units for PCA
        end

        % Reference time axis from first row (should be 0→UseMaxDur_s)
        time_ref = catRows.new_vecTime_MUA{1};
        if isempty(time_ref) || ~isnumeric(time_ref)
            warning('Empty/non-numeric new_vecTime_MUA for %s × %s; skipping.', mv, dz);
            continue;
        end
        time_ref = time_ref(:);

        nUnits = height(catRows);
        X = nan(numel(time_ref), nUnits);
        keepUnit = false(1, nUnits);

        for u = 1:nUnits
            dVec = catRows.vecD_MUA{u};
            n_tVec = catRows.new_vecTime_MUA{u};

            if isempty(dVec) || isempty(n_tVec), continue; end
            dVec = double(dVec(:));
            n_tVec = double(n_tVec(:));

            if numel(dVec) < 3 || numel(n_tVec) ~= numel(dVec)
                continue;
            end

            try
                X(:,u) = interp1(n_tVec, dVec, time_ref, 'linear', 'extrap');
                keepUnit(u) = true;
            catch ME
                warning('interp1 failed for MUA unit %d in %s × %s: %s', ...
                    u, mv, dz, ME.message);
            end
        end

        X = X(:, keepUnit);
        if size(X,2) < 2
            continue;
        end

        % Remove timepoints with all NaNs
        goodTime = any(isfinite(X), 2);
        X = X(goodTime, :);
        t_use = time_ref(goodTime);

        % Fill remaining NaNs with column means
        for j = 1:size(X,2)
            nanMask = ~isfinite(X(:,j));
            if any(nanMask)
                mu_j = mean(X(~nanMask,j), 'omitnan');
                X(nanMask,j) = mu_j;
            end
        end

        % PCA: observations = units, variables = timepoints
        [coeff, score, latent] = pca(X', 'Centered', true);

        pc1 = coeff(:,1);

        % Flip sign so that post-onset mean is positive
        % postMask = t_use >= 0 & t_use <= min(0.5*max(t_use), 0.8*max(t_use));
        postMask = t_use >= 0 & t_use <= (max(t_use) * 0.5);
        if any(postMask)
            if mean(pc1(postMask)) < 0
                pc1 = -pc1;
                score(:,1) = -score(:,1);
            end
        end

        rawKey = sprintf('%s_%s', dz, mv);
        key = matlab.lang.makeValidName(rawKey);

        PC1_MUA.(key).t      = t_use;
        PC1_MUA.(key).pc1    = pc1;
        PC1_MUA.(key).nUnits = size(X,2);
        PC1_MUA.(key).latent = latent;
        PC1_MUA.(key).scores = score;
        PC1_MUA.(key).labels = catRows.PrettyLabel(keepUnit);

        fprintf('[PCA MUA] %s × %s: %d units, %d time points\n', ...
            mv, dz, size(X,2), numel(t_use));
    end
end


%% -------- PLOT: PC1 per depth (rows) & MoveType (colors) --------

mtOrder = {'HAND OC','HAND PS','ARM EF','REST'};
mtOrder = intersect(mtOrder, moveTypes, 'stable');

% --- Compute global y-limits for MUA PC1 across all depths × MoveTypes ---
allPC1_MUA = [];
for d = 1:numel(depths)
    dz = depths{d};
    for m = 1:numel(mtOrder)
        mv = mtOrder{m};
        rawKey = sprintf('%s_%s', dz, mv);
        key    = matlab.lang.makeValidName(rawKey);

        if isfield(PC1_MUA, key) && ~isempty(PC1_MUA.(key).pc1)
            allPC1_MUA = [allPC1_MUA; PC1_MUA.(key).pc1(:)];
        end
    end
end

if ~isempty(allPC1_MUA)
    yPad   = 0.015;
    yLimsM = [min(allPC1_MUA)-yPad, max(allPC1_MUA)+yPad];
else
    yLimsM = [];
end
% -------------------------------------------------------------------------

hFig = figure('Color','w','Position',[100 100 950 800]);
tlo = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

for d = 1:numel(depths)
    dz = depths{d};
    ax = nexttile; hold(ax,'on'); grid(ax,'on');

    for m = 1:numel(mtOrder)
        mv = mtOrder{m};

        rawKey = sprintf('%s_%s', dz, mv);
        key    = matlab.lang.makeValidName(rawKey);

        if ~isfield(PC1_MUA, key), continue; end
        if isempty(PC1_MUA.(key).t), continue; end

        t_use = PC1_MUA.(key).t;
        pc1   = PC1_MUA.(key).pc1;

        if isKey(moveColorMap, mv)
            c = moveColorMap(mv);
        else
            c = fallbackCol;
        end

        plot(ax, t_use, pc1, 'LineWidth', 2, 'Color', c, 'DisplayName', mv);
    end

    xlabel(ax, 'Time from movement onset (s)');
    ylabel(ax, 'MUA ZETA PC1 (a.u.)');
    if isKey(depthNames, dz)
        title(ax, sprintf('MUA PC1 of ZETA deviation | %s', depthNames(dz)));
    else
        title(ax, sprintf('MUA PC1 of ZETA deviation | depth %s', dz));
    end

    % % color subplot titles by depth
    % if isKey(depthNames, dz)
    %     tt = title(ax, sprintf('MUA PC1 of ZETA deviation | %s', depthNames(dz)));
    % else
    %     tt = title(ax, sprintf('MUA PC1 of ZETA deviation | depth %s', dz));
    % end
    % if isKey(depthColorMap, dz)
    %     set(tt,'Color', depthColorMap(dz));
    % end


    xline(ax, 0, 'k--', 'LineWidth', 1);
    xlim(ax, [-0.2 1.6]);  % match SU/MUA PC1 plots

    % ----- apply uniform y-limits across all MUA tiles -----
    if ~isempty(yLimsM)
        ylim(ax, yLimsM);
    end
    % -----------------------------------------------
end

% One legend for all
lgd = legend(tlo.Children(end), mtOrder, 'Location','northeastoutside');    % eastoutside; % change location to top right corner (in line with title)
title(tlo, 'PC1 of MUA ZETA temporal deviation per MoveType × STN depth');


%% Save figure

outPNG = fullfile(U.SavePath, 'ZETA_PC1_MUA_byDepth.png');
exportgraphics(hFig, outPNG, 'Resolution',300);

outFIG = fullfile(U.SavePath, 'ZETA_PC1_MUA_byDepth.fig');
savefig(hFig, outFIG);


end