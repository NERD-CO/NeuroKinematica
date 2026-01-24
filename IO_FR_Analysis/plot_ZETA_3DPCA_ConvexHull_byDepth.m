function out = plot_ZETA_3DPCA_ConvexHull_byDepth(MasterZETA, varargin)
% plot_ZETA_3DPCA_ConvexHull_byDepth
%
% GOALS:
%   1) For each STN depth (t/c/b), compute PCA on time-aligned ZETA deviation
%      vectors across units (units = observations; timepoints = variables).
%   2) Plot 3D PC1–PC3 scatter per active MoveType (HAND OC/PS/ARM EF),
%      and overlay convex-hull projections onto the XY/XZ/YZ "walls".
%   3) Compute and plot silhouette values (PC space) per MoveType per depth,
%      and return a summary table of silhouette statistics.
%
% REQUIRED columns in MasterZETA (names as you used):
%   - MoveType (categorical or string)
%   - Depth (categorical or string; values include 't','c','b')
%   - ZETA_vecD (cell array of numeric vectors)
%   - ZETA_vecT (cell array of numeric vectors)
%   - PSTH_TimeCenters_s (cell array of numeric vectors)  OR use ZETA_vecT as reference
%
% OUTPUT:
%   out struct fields:
%     .PerDepth(d).Depth
%     .PerDepth(d).PCAScores  (table with PC1,PC2,PC3, MoveType)
%     .PerDepth(d).Coeff      (PCA coeff)
%     .PerDepth(d).Explained  (variance explained %)
%     .PerDepth(d).Silhouette (table of silhouette per unit)
%     .SilSummaryTbl          (summary across depth x MoveType)
%
% Example:
%   out = plot_ZETA_3DPCA_ConvexHull_byDepth(MasterZETA, ...
%       'SaveDir', fullfile(pwd,'PCA'), 'DoSave', true);

% -----------------------
% Parse inputs
% -----------------------
p = inputParser;

p.addParameter('DepthOrder', {'t','c','b'}, @(x) iscell(x) && ~isempty(x));
p.addParameter('ActiveMoveOrder', {'HAND OC','HAND PS','ARM EF'}, @(x) iscell(x) && ~isempty(x));

p.addParameter('SaveDir', '', @(x)ischar(x)||isstring(x));
p.addParameter('DoSave', false, @islogical);
p.addParameter('FigPrefix', 'ZETA_PCA3D', @(x)ischar(x)||isstring(x));

% Time alignment
p.addParameter('TimeSource', 'PSTH_TimeCenters_s', @(s) any(strcmpi(string(s), ["PSTH_TimeCenters_s","ZETA_vecT"])));
p.addParameter('InterpMethod', 'linear', @(s)ischar(s)||isstring(s));

% PCA sign convention (optional)
p.addParameter('FlipPC1PostOnset', true, @islogical);
p.addParameter('PostOnsetWindow_s', [0 0.8], @(x)isnumeric(x)&&numel(x)==2);

% Hull projection appearance
p.addParameter('HullFaceAlpha', 0.10, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('HullEdgeAlpha', 0.80, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('HullLineWidth', 1.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% Silhouette
p.addParameter('DoSilhouettePlot', true, @islogical);
p.addParameter('SilhouetteMethod', 'euclidean', @(s)ischar(s)||isstring(s));

% Plot geometry
p.addParameter('View', [-40 20], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('MarkerSize', 30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MarkerAlpha', 0.75, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);

p.parse(varargin{:});
U = p.Results;

if U.DoSave && isempty(U.SaveDir)
    U.SaveDir = fullfile(pwd, 'PCA_Plots');
end
if U.DoSave && ~exist(U.SaveDir,'dir')
    mkdir(U.SaveDir);
end

% -----------------------
% Sanity checks / required columns
% -----------------------
reqVars = {'MoveType','Depth','ZETA_vecD','ZETA_vecT'};
missing = setdiff(reqVars, MasterZETA.Properties.VariableNames);
if ~isempty(missing)
    error('MasterZETA missing required variables: %s', strjoin(missing, ', '));
end
if strcmpi(U.TimeSource,'PSTH_TimeCenters_s')
    if ~ismember('PSTH_TimeCenters_s', MasterZETA.Properties.VariableNames)
        warning('TimeSource=PSTH_TimeCenters_s requested but column missing. Falling back to ZETA_vecT.');
        U.TimeSource = 'ZETA_vecT';
    end
end

% Normalize types
MT = MasterZETA;
MT.MoveType = string(MT.MoveType);
MT.Depth    = string(MT.Depth);

% Restrict to active MoveTypes only
activeMoves = string(U.ActiveMoveOrder);
MT = MT(ismember(MT.MoveType, activeMoves), :);

% Require non-empty ZETA
hasZ = ~cellfun(@isempty, MT.ZETA_vecD) & ~cellfun(@isempty, MT.ZETA_vecT);
MT = MT(hasZ, :);

if isempty(MT)
    error('No usable rows after filtering to active MoveTypes with non-empty ZETA vectors.');
end

% -----------------------
% JNE color scheme
% -----------------------
greenShades = ([ ...
    217-20,240-15,211-20;     % HAND OC
    127-20,191-15,123-20;     % HAND PS
    27-20, 120-15, 55-20] ... % ARM EF
    ./ 255);

moveColorMap = containers.Map( ...
    {'HAND OC','HAND PS','ARM EF'}, ...
    {greenShades(1,:), greenShades(2,:), greenShades(3,:)} );

depthNames = containers.Map({'t','c','b'}, {'dorsal STN','central STN','ventral STN'});


%% ==============================
% Main function: per-depth loop
% -------------------------------

depthOrder = string(U.DepthOrder);

out = struct;
out.PerDepth = struct([]);
silSummaryRows = {};

out.Stats = struct();
out.Stats.KW    = struct('tbl', table(), 'mc', table());
out.Stats.ANOVA = struct('tbl', table(), 'mc', table());

% per-depth loop
for d = 1:numel(depthOrder)
    dz = char(depthOrder(d));
    Md = MT(MT.Depth == string(dz), :);
    if height(Md) < 3
        % too few units to do PC1-3 sensibly
        continue;
    end

    % ---- Build time-aligned matrix X (time x units) ----
    [X, tRef, mvPerUnit] = build_timeAlignedMatrix(Md, U.TimeSource, U.InterpMethod);

    if size(X,2) < 3
        continue;
    end

    %% ---- PCA ----
    % PCA expects observations as rows -> units x time
    [coeff, score, latent, ~, explained] = pca(X', 'Centered', true);

    % Ensure we have 3 PCs for plotting
    if size(score,2) < 3
        score(:,3) = 0;
    end

    % Optional: flip PC1 so post-onset mean deflection is positive
    if U.FlipPC1PostOnset
        postMask = (tRef >= U.PostOnsetWindow_s(1)) & (tRef <= U.PostOnsetWindow_s(2));
        if any(postMask)
            pc1 = coeff(:,1);
            if mean(pc1(postMask), 'omitnan') < 0
                coeff(:,1) = -coeff(:,1);
                score(:,1) = -score(:,1);
            end
        end
    end

    PC1 = score(:,1);
    PC2 = score(:,2);
    PC3 = score(:,3);

    % Store PCA scores table
    Tscore = table(PC1,PC2,PC3, mvPerUnit, 'VariableNames',{'PC1','PC2','PC3','MoveType'});


    %% ---- Figure: 3D PCA + hull projections (+ silhouette tile optionally) ----
    if U.DoSilhouettePlot
        hFig = figure('Color','w','Units','pixels','Position',[120 80 1200 650]);
        tlo = tiledlayout(hFig, 1, 2, 'TileSpacing','compact','Padding','compact');
        ax3 = nexttile(tlo, 1); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
        axS = nexttile(tlo, 2); hold(axS,'on'); grid(axS,'on'); box(axS,'on');
    else
        hFig = figure('Color','w','Units','pixels','Position',[120 80 900 650]);
        ax3 = axes('Parent', hFig); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
        axS = [];
    end

    %% --- 3D PCA scatter by MoveType ---

    mvOrder = U.ActiveMoveOrder;
    hLeg = gobjects(numel(mvOrder),1); % Collect handles for legend

    % Determine "wall" positions based on overall PC ranges
    [wallInfo, limInfo] = computeWalls(PC1,PC2,PC3);

    for mIdx = 1:numel(mvOrder)
        mv = mvOrder{mIdx};
        idx = (mvPerUnit == string(mv));
        if nnz(idx) < 3
            continue;
        end
        c = moveColorMap(mv);

        % % 3D scatter of PC scores
        % hLeg(mIdx) = scatter3(ax3, PC1(idx), PC2(idx), PC3(idx), ...
        %     U.MarkerSize, c, 'filled', ...
        %     'MarkerFaceAlpha', U.MarkerAlpha, 'MarkerEdgeColor','none', ...
        %     'DisplayName', mv);
        % hold on
        % % trajectory lines
        % plot3(ax3, PC1(idx), PC2(idx), PC3(idx), 'Color',c, 'LineWidth', 3)

        % -----------------------------------------------------------------
        %% 3D scatter of smoothed PC scores and trajectory lines
        num_points = sum(idx);
        original_indices = 1:num_points;
        new_indices = linspace(1, num_points, num_points * 10);

        % spline interpolation
        x_smooth = interp1(original_indices, PC1(idx), new_indices, 'spline');
        y_smooth = interp1(original_indices, PC2(idx), new_indices, 'spline');
        z_smooth = interp1(original_indices, PC3(idx), new_indices, 'spline');

        % 3D scatter of PC scores
        hLeg(mIdx) = scatter3(ax3, x_smooth, y_smooth, z_smooth, ...
            10, c, 'filled', 'MarkerFaceAlpha',0.3);
        % U.MarkerSize, c, 'filled', ...
        % 'MarkerFaceAlpha', U.MarkerAlpha, 'MarkerEdgeColor','none', ...
        % 'DisplayName', mv);
        hold on
        % trajectory lines
        plot3(ax3, x_smooth, y_smooth, z_smooth, 'Color',c, 'LineWidth', 3)
        % -----------------------------------------------------------------

        % Convex hull projections on XY, XZ, YZ planes
        plotHullProjections(ax3, PC1(idx), PC2(idx), PC3(idx), c, ...
            wallInfo, U.HullFaceAlpha, U.HullEdgeAlpha, U.HullLineWidth);
    end

    % Axes formatting
    xlabel(ax3,'PC1'); ylabel(ax3,'PC2'); zlabel(ax3,'PC3');
    if isKey(depthNames, dz)
        ttl = sprintf('3D PCA + Hull projections | %s', depthNames(dz));
    else
        ttl = sprintf('3D PCA + Hull projections | depth %s', dz);
    end
    title(ax3, ttl, 'Interpreter','none');

    view(ax3, U.View); % azimuth, elevation
    axis(ax3,'vis3d');  % preserve aspect ratio under rotation

    % center axes around data with small margin (0.1)
    xlim(ax3, limInfo(1,:)); ylim(ax3, limInfo(2,:)); zlim(ax3, limInfo(3,:));

    % Small padding so labels never touch figure edges
    ax3.LooseInset = max(ax3.TightInset, 0.05);

    % Legend
    ok = isgraphics(hLeg);
    if any(ok)
        legend(ax3, hLeg(ok), mvOrder(ok), 'Location','northeastoutside', 'FontSize', 11);
    end



    %% Cluster Eval:
    %  --- Silhouette (per unit) in PC space ---
    silTbl = table;

    if U.DoSilhouettePlot
        mvPerUnit = string(mvPerUnit(:));
        [clust,~,cl_idx] = unique(mvPerUnit, 'stable');

        % Need >=2 clusters and at least 2 points per cluster to be meaningful
        counts = accumarray(cl_idx, 1);
        okSil = (numel(clust) >= 2) && (sum(counts >= 2) >= 2) && (size([PC1 PC2 PC3],1) >= 4);

        if ~okSil
            silVals = nan(size([PC1 PC2 PC3],1),1);
            return;
        end
        try
            silVals = silhouette([PC1 PC2 PC3], mvPerUnit, U.SilhouetteMethod);
            okSil = true;
        catch
            silVals = nan(size([PC1 PC2 PC3],1),1);
            okSil = false;
        end

        if okSil
            silTbl = table(mvPerUnit, silVals, 'VariableNames', {'MoveType','Silhouette'});

            % silhouette plot (MATLAB's silhouette draws its own plot)
            axes(axS);
            cla(axS);
            silhouette([PC1 PC2 PC3], mvPerUnit, U.SilhouetteMethod);
            title(axS, sprintf('Silhouette | depth %s', dz), 'Interpreter','none');
            xlabel(axS, 'Silhouette value'); ylabel(axS, 'Cluster (MoveType)');

            % Summaries per MoveType
            for mIdx = 1:numel(mvOrder)
                mv = mvOrder{mIdx};
                mm = mvPerUnit == string(mv);
                if ~any(mm), continue; end
                silSummaryRows(end+1,:) = {dz, mv, nnz(mm), ...
                    mean(silVals(mm),'omitnan'), ...
                    median(silVals(mm),'omitnan'), ...
                    std(silVals(mm),'omitnan')};
            end
        else
            title(axS, sprintf('Silhouette not computed (insufficient groups) | depth %s', dz), 'Interpreter','none');
        end


        %% Stats:
        % Assess skewness
        x = silVals(isfinite(silVals));
        sk = skewness(x);
        m  = mean(x);
        md = median(x);
        q  = quantile(x,[0.25 0.5 0.75]);
        fprintf('Skewness = %.3f | mean=%.3f | median=%.3f | Q1=%.3f Q3=%.3f\n', ...
            sk, m, md, q(1), q(3));

        % Plot histograms per depth
        mvPerUnit = string(mvPerUnit(:)); % string array of MoveType per unit
        mvOrder = ["HAND OC","HAND PS","ARM EF"];   % enforce consistent ordering

        factors = categorical(mvPerUnit, mvOrder, 'Ordinal', false);
        cats = categories(factors);

        % --- Histogram per MoveType (color-coded) ---
        figure;
        tiledlayout(1, numel(cats), 'TileSpacing','compact','Padding','compact');

        for i = 1:numel(cats)
            nexttile;
            thisCat = cats{i};
            idx = (factors == thisCat);

            h = histogram(silVals(idx), 15);  % <-- capture handle

            % Color by MoveType
            mvKey = char(string(thisCat));
            if isKey(moveColorMap, mvKey)
                h.FaceColor = moveColorMap(mvKey);
            else
                h.FaceColor = [0.6 0.6 0.6]; % fallback
            end

            % Optional aesthetics
            h.EdgeColor = 'none';
            h.FaceAlpha = 0.85;

            xlim([-1 1]);
            title(mvKey, 'Interpreter','none');
            xlabel('Silhouette'); ylabel('Count');

            xline(median(silVals(idx), 'omitnan'), '-k', 'median');
        end


        %% Kruskal–Wallis + post-hoc tests per depth

        [p_kw, tbl_kw, stats_kw] = kruskalwallis(silVals, factors, 'off');
        fprintf('Kruskal–Wallis p = %.4g\n', p_kw);

        % convert cell array to table with variable names
        tblKW = cell2table(tbl_kw(2:end,:), ...
            'VariableNames', matlab.lang.makeValidName(string(tbl_kw(1,:))));

        % Add metadata
        tblKW.Depth = repmat(string(dz), height(tblKW), 1);
        tblKW = movevars(tblKW, 'Depth', 'Before', 1);
        tblKW.pValue_overall = repmat(p_kw, height(tblKW), 1);

        % Post-hoc pairwise comparisons (Bonferroni)
        mc_kw = multcompare(stats_kw, 'CType','bonferroni', 'Display','off');
        mcTbl_kw = array2table(mc_kw, 'VariableNames', {'Group1','Group2','LowerCI','Diff','UpperCI','pValue'});

        grpNames = categories(factors);
        mcTbl_kw.Group1 = grpNames(mcTbl_kw.Group1);
        mcTbl_kw.Group2 = grpNames(mcTbl_kw.Group2);
        mcTbl_kw.Depth  = repmat(string(dz), height(mcTbl_kw), 1);
        mcTbl_kw = movevars(mcTbl_kw, 'Depth', 'Before', 1);

        disp(mcTbl_kw);

        %% One-way ANOVA + post hoc tests per depth

        [p_aov, tbl_aov, stats_aov] = anova1(silVals, factors, 'off');
        fprintf('One-way ANOVA p = %.4g\n', p_aov);

        % convert cell array to table with variable names
        tblAOV = cell2table(tbl_aov(2:end,:), ...
            'VariableNames', matlab.lang.makeValidName(string(tbl_aov(1,:))));

        % Add metadata
        tblAOV.Depth = repmat(string(dz), height(tblAOV), 1);
        tblAOV = movevars(tblAOV, 'Depth', 'Before', 1);
        tblAOV.pValue_overall = repmat(p_aov, height(tblAOV), 1);

        % Post-hoc pairwise comparisons (Bonferroni)
        mc_aov = multcompare(stats_aov, 'CType','bonferroni', 'Display','off');
        mcTbl_aov = array2table(mc_aov, 'VariableNames', {'Group1','Group2','LowerCI','Diff','UpperCI','pValue'});

        grpNames = categories(factors);
        mcTbl_aov.Group1 = grpNames(mcTbl_aov.Group1);
        mcTbl_aov.Group2 = grpNames(mcTbl_aov.Group2);
        mcTbl_aov.Depth  = repmat(string(dz), height(mcTbl_aov), 1);
        mcTbl_aov = movevars(mcTbl_aov, 'Depth', 'Before', 1);

        disp(mcTbl_aov);
    end


    %% Save if requested
    if U.DoSave
        base = sprintf('%s_%s', char(string(U.FigPrefix)), dz);
        print(hFig, fullfile(U.SaveDir, [base '.png']), '-dpng', '-r300');
        savefig(hFig, fullfile(U.SaveDir, [base '.fig']));
    end

    % Store PCA and Silhouette outputs
    out.PerDepth(end+1).Depth      = dz;
    out.PerDepth(end).Coeff        = coeff;
    out.PerDepth(end).ExplainedPct = explained;
    out.PerDepth(end).PCA_Scores   = Tscore;
    out.PerDepth(end).Silhouette   = silTbl;

    % Append Stat outputs per-depth
    out.PerDepth(end).KW.tbl = tblKW;
    out.PerDepth(end).KW.mc  = mcTbl_kw;

    out.PerDepth(end).ANOVA.tbl = tblAOV;
    out.PerDepth(end).ANOVA.mc  = mcTbl_aov;


    % Aggregated stats
    out.Stats.KW.tbl  = [out.Stats.KW.tbl;  tblKW];
    out.Stats.KW.mc   = [out.Stats.KW.mc;   mcTbl_kw];
    out.Stats.ANOVA.tbl = [out.Stats.ANOVA.tbl; tblAOV];
    out.Stats.ANOVA.mc  = [out.Stats.ANOVA.mc;  mcTbl_aov];
end

% Summary table
if ~isempty(silSummaryRows)
    out.SilSummaryTbl = cell2table(silSummaryRows, ...
        'VariableNames', {'Depth','MoveType','N','SilMean','SilMedian','SilStd'});
else
    out.SilSummaryTbl = table;
end

end % main function


%% ============================
% Helpers
% ============================

function [X, tRef, mvPerUnit] = build_timeAlignedMatrix(Md, timeSource, interpMethod)
% Build X: time x units from ZETA_vecD interpolated to a shared timebase.
% mvPerUnit: MoveType label per unit/row.

% Determine reference time axis
switch lower(string(timeSource))
    case "psth_timecenters_s"
        tRef0 = Md.PSTH_TimeCenters_s{find(~cellfun(@isempty, Md.PSTH_TimeCenters_s), 1, 'first')};
    otherwise
        % use first ZETA_vecT as reference
        tRef0 = Md.ZETA_vecT{find(~cellfun(@isempty, Md.ZETA_vecT), 1, 'first')};
end

tRef0 = double(tRef0(:));
tRef0 = tRef0(isfinite(tRef0));
if numel(tRef0) < 6
    error('Reference time vector too short / invalid.');
end

% Enforce strictly increasing unique reference timebase
[tRef, ia] = unique(tRef0, 'stable');
tRef = tRef(:);

nU = height(Md);
X  = nan(numel(tRef), nU);
mvPerUnit = string(Md.MoveType);

for u = 1:nU
    dVec = Md.ZETA_vecD{u};
    tVec = Md.ZETA_vecT{u};

    if isempty(dVec) || isempty(tVec), continue; end
    dVec = double(dVec(:));
    tVec = double(tVec(:));

    good = isfinite(dVec) & isfinite(tVec);
    dVec = dVec(good);
    tVec = tVec(good);

    if numel(tVec) < 6 || numel(dVec) ~= numel(tVec)
        continue;
    end

    % Make sample points unique+increasing for interp1 safety
    [tU, iu] = unique(tVec, 'stable');
    dU = dVec(iu);

    if numel(tU) < 6, continue; end

    try
        X(:,u) = interp1(tU, dU, tRef, interpMethod, 'extrap');
    catch
        % If interpolation fails, leave NaNs
    end
end

% Drop units with all-NaN
keepU = any(isfinite(X), 1);
X = X(:, keepU);
mvPerUnit = mvPerUnit(keepU);

% Fill remaining NaNs per unit with that unit's mean (so PCA doesn't choke)
for j = 1:size(X,2)
    nanMask = ~isfinite(X(:,j));
    if any(nanMask)
        mu = mean(X(~nanMask,j), 'omitnan');
        X(nanMask,j) = mu;
    end
end

% Drop timepoints that are NaN across all units (rare after fill, but safe)
goodT = any(isfinite(X), 2);
X = X(goodT,:);
tRef = tRef(goodT);

end


function [wallInfo, limInfo] = computeWalls(x,y,z)
% Compute wall locations and axis limits with margins.
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
zmin = min(z); zmax = max(z);

mx = 0.10 * max(1e-6, xmax-xmin);
my = 0.10 * max(1e-6, ymax-ymin);
mz = 0.10 * max(1e-6, zmax-zmin);

% --- Walls ---
% x_wall = xmin - mx;   % YZ, front-right wall
x_wall = xmax + (1+mx); % flip --> YZ projection on BACK-right plane
% y_wall = ymax + my;   % XZ, back-left wall
y_wall = ymax + (1+my); % XZ projection on back-left plane
% z_wall = zmin - mz;   % XY, base/floor wall
z_wall = zmin - (1+mz); % XY projection on bottom plane

wallInfo = [x_wall, y_wall, z_wall];

% --- Axis limits ---
x_lim_min = xmin - (1+mx);
y_lim_min = ymin - (1+my);
z_lim_max = zmax + (1+mz);

limInfo = [...
    x_lim_min, x_wall;   % xlim
    y_lim_min, y_wall;   % ylim
    z_wall, z_lim_max];  % zlim
end


function plotHullProjections(ax, x, y, z, c, wallInfo, faceAlpha, edgeAlpha, lw)
% Plot convex hull projections onto XY (z=z_wall), XZ (y=y_wall), YZ (x=x_wall)
% Uses 2D convhull. Requires >=3 points; gracefully skips on failure.

x_wall = wallInfo(1);
y_wall = wallInfo(2);
z_wall = wallInfo(3);

% Ensure column vectors
x = x(:); y = y(:); z = z(:);

% Need at least 3 non-collinear points in each 2D plane; convhull will error if degenerate.
try
    if numel(x) >= 3
        k_xy = convhull(x, y);
        z_proj = repmat(z_wall, size(k_xy));
        p1 = patch(ax, x(k_xy), y(k_xy), z_proj, c, ...
            'FaceAlpha', faceAlpha, 'EdgeColor', c, 'LineWidth', lw);
        p1.EdgeAlpha = edgeAlpha;
        p1.FaceColor = 'none'; % outline-only; comment out for filled translucent
    end
catch
end

try
    if numel(x) >= 3
        k_xz = convhull(x, z);
        y_proj = repmat(y_wall, size(k_xz));
        p2 = patch(ax, x(k_xz), y_proj, z(k_xz), c, ...
            'FaceAlpha', faceAlpha, 'EdgeColor', c, 'LineWidth', lw);
        p2.EdgeAlpha = edgeAlpha;
        p2.FaceColor = 'none';
    end
catch
end

try
    if numel(y) >= 3
        k_yz = convhull(y, z);
        x_proj = repmat(x_wall, size(k_yz));
        p3 = patch(ax, x_proj, y(k_yz), z(k_yz), c, ...
            'FaceAlpha', faceAlpha, 'EdgeColor', c, 'LineWidth', lw);
        p3.EdgeAlpha = edgeAlpha;
        p3.FaceColor = 'none';
    end
catch
end

end


