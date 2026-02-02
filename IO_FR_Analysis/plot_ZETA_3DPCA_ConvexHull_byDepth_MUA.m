function out = plot_ZETA_3DPCA_ConvexHull_byDepth_MUA(MasterZETA_MUA, varargin)
% plot_ZETA_3DPCA_ConvexHull_byDepth_MUA
%
% MUA version of your SU pipeline:
%   1) Per depth (t/c/b): time-align MUA ZETA deviation vectors (vecD_MUA) across units
%   2) PCA (units as observations; timepoints as variables)
%   3) Plot 3D PC1–PC3 scatter by MoveType + convex hull projections on "walls"
%   4) Compute silhouette in PC space by MoveType + summary stats per MoveType
%   5) Kruskal–Wallis + ANOVA on silhouette values across MoveTypes, w/ post-hoc tests
%
% REQUIRED columns in MasterZETA_MUA:
%   - MoveType
%   - Depth (values: 't','c','b')
%   - vecD_MUA (cell array numeric vector)
%   - new_vecTime_MUA (cell array numeric vector)
%
% OUTPUT:
%   out.PerDepth(d).Depth
%   out.PerDepth(d).Coeff
%   out.PerDepth(d).ExplainedPct
%   out.PerDepth(d).PCA_Scores
%   out.PerDepth(d).Silhouette
%   out.PerDepth(d).KW.tbl / .KW.mc
%   out.PerDepth(d).ANOVA.tbl / .ANOVA.mc
%   out.Stats.KW.tbl / .KW.mc           (aggregated across depths)
%   out.Stats.ANOVA.tbl / .ANOVA.mc     (aggregated across depths)
%   out.SilSummaryTbl                   (Depth x MoveType summary)
%
% Example:
%   out = plot_ZETA_3DPCA_ConvexHull_byDepth_MUA(MasterZETA_MUA, ...
%       'DoSave', true, 'SaveDir', fullfile(pwd,'MUA_PCA'));

% -----------------------
% Parse inputs
% -----------------------
p = inputParser;

p.addParameter('DepthOrder', {'t','c','b'}, @(x) iscell(x) && ~isempty(x));
p.addParameter('ActiveMoveOrder', {'HAND OC','HAND PS','ARM EF'}, @(x) iscell(x) && ~isempty(x));

p.addParameter('SaveDir', '', @(x)ischar(x)||isstring(x));
p.addParameter('FigPrefix', 'ZETA_PCA3D_MUA', @(x)ischar(x)||isstring(x));

% Time alignment
p.addParameter('InterpMethod', 'linear', @(s)ischar(s)||isstring(s));

% PCA sign convention (optional)
p.addParameter('FlipPC1PostOnset', true, @islogical);
p.addParameter('PostOnsetWindow_s', [0 0.8], @(x)isnumeric(x)&&numel(x)==2);

% Hull projection appearance
p.addParameter('HullFaceAlpha', 0.10, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('HullEdgeAlpha', 0.80, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('HullLineWidth', 1.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% Silhouette + plotting
p.addParameter('DoSilhouettePlot', true, @islogical);
p.addParameter('SilhouetteMethod', 'euclidean', @(s)ischar(s)||isstring(s));

% If true, also make a separate silhouette histogram figure (overlay + stacked)
p.addParameter('DoSilhouetteHistFigs', false, @islogical);

% Plot geometry
p.addParameter('View', [-40 20], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('MarkerSize', 30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MarkerAlpha', 0.75, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);

p.parse(varargin{:});
U = p.Results;

if isempty(U.SaveDir)
    U.SaveDir = fullfile(pwd, '3D PCA Stats');
end
if ~exist(U.SaveDir,'dir')
    mkdir(U.SaveDir);
end

% -----------------------
% Sanity checks / required columns
% -----------------------
reqVars = {'MoveType','Depth','vecD_MUA','new_vecTime_MUA'};
missing = setdiff(reqVars, MasterZETA_MUA.Properties.VariableNames);
if ~isempty(missing)
    error('MasterZETA_MUA missing required variables: %s', strjoin(missing, ', '));
end

% Normalize types
MT = MasterZETA_MUA;
MT.MoveType = string(MT.MoveType);
MT.Depth    = string(MT.Depth);

% Restrict to active MoveTypes only (ignore REST etc.)
activeMoves = string(U.ActiveMoveOrder);
MT = MT(ismember(MT.MoveType, activeMoves), :);

% Require non-empty vectors
hasZ = ~cellfun(@isempty, MT.vecD_MUA) & ~cellfun(@isempty, MT.new_vecTime_MUA);
MT = MT(hasZ, :);

if isempty(MT)
    error('No usable rows after filtering to active MoveTypes with non-empty MUA vectors.');
end

% -----------------------
% Color scheme (match your SU)
% -----------------------
% greenShades = ([ ...
%     217-20,240-15,211-20;     % HAND OC
%     127-20,191-15,123-20;     % HAND PS
%     27-20, 120-15, 55-20] ... % ARM EF
%     ./ 255);
JNE_move = ([ ...
    128,128,128;    % REST  (grey)
    38,116,183;     % HAND OC  (blue)
    53,183,121;     % HAND PS  (green/teal)
    243,120,98] ...   % ARM EF   (coral)
    ./ 255);  % /255 = standard

moveColorMap = containers.Map( ...
    {'HAND OC','HAND PS','ARM EF'}, ...
    {JNE_move(2,:), JNE_move(3,:), JNE_move(4,:)} );

depthNames = containers.Map({'t','c','b'}, {'dorsal STN','central STN','ventral STN'});

% -----------------------
% Outputs
% -----------------------
depthOrder = string(U.DepthOrder);

out = struct;
out.PerDepth = struct([]);
silSummaryRows = {};

out.Stats = struct();
out.Stats.KW    = struct('tbl', table(), 'mc', table());
out.Stats.ANOVA = struct('tbl', table(), 'mc', table());

% ======================
% Main per-depth loop
% ======================
for d = 1:numel(depthOrder)
    dz = char(depthOrder(d));
    Md = MT(MT.Depth == string(dz), :);

    if height(Md) < 3
        continue;
    end

    % ---- Build time-aligned matrix X (time x units) ----
    [X, tRef, mvPerUnit] = build_timeAlignedMatrix_MUA(Md, U.InterpMethod);

    if size(X,2) < 3
        continue;
    end

    % ---- PCA (observations=units) ----
    [coeff, score, ~, ~, explained] = pca(X', 'Centered', true);

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

    Tscore = table(PC1,PC2,PC3, mvPerUnit, 'VariableNames',{'PC1','PC2','PC3','MoveType'});

    % ---- Figure: 3D PCA + silhouette tile ----
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

    % ---- 3D PCA scatter + hulls ----
    mvOrder = string(U.ActiveMoveOrder);

    [wallInfo, limInfo] = computeWalls(PC1,PC2,PC3); % includes "mirrored" YZ wall placement like your SU

    for mIdx = 1:numel(mvOrder)
        mv = mvOrder(mIdx);
        idx = (mvPerUnit == mv);
        if nnz(idx) < 3, continue; end

        c = moveColorMap(char(mv));

        % scatter (light) + trajectory (thicker)
        num_points = sum(idx);
        oi = 1:num_points;
        ni = linspace(1, num_points, num_points * 10);

        x_smooth = interp1(oi, PC1(idx), ni, 'spline');
        y_smooth = interp1(oi, PC2(idx), ni, 'spline');
        z_smooth = interp1(oi, PC3(idx), ni, 'spline');

        scatter3(ax3, x_smooth, y_smooth, z_smooth, 10, c, 'filled', 'MarkerFaceAlpha',0.3);
        plot3(ax3, x_smooth, y_smooth, z_smooth, 'Color', c, 'LineWidth', 3);

        % plotHullProjections(ax3, PC1(idx), PC2(idx), PC3(idx), c, ...
        %     wallInfo, U.HullFaceAlpha, U.HullEdgeAlpha, U.HullLineWidth);
    end

    xlabel(ax3,'PC1'); ylabel(ax3,'PC2'); zlabel(ax3,'PC3');
    title(ax3, sprintf('3D PCA | %s (MUA)', depthNames(dz)), 'Interpreter','none'); % 3D PCA + Hull projections 

    view(ax3, U.View);
    axis(ax3,'vis3d');
    xlim(ax3, limInfo(1,:)); ylim(ax3, limInfo(2,:)); zlim(ax3, limInfo(3,:));
    ax3.LooseInset = max(ax3.TightInset, 0.05);

    % Legend (dummy scatters so colors are opaque)
    hLeg2 = gobjects(numel(mvOrder),1);
    for mIdx = 1:numel(mvOrder)
        mv = char(mvOrder(mIdx));
        c  = moveColorMap(mv);
        hLeg2(mIdx) = scatter3(ax3, NaN, NaN, NaN, 60, c, 'filled', ...
            'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor','none', 'DisplayName', mv);
    end
    lgd = legend(ax3, hLeg2, cellstr(mvOrder), 'Location','northeastoutside', 'FontSize', 11);
    lgd.AutoUpdate = 'off';

    
    %% Silhouette + stats
    silTbl = table;
    tblKW = table; mcTbl_kw = table;
    tblAOV = table; mcTbl_aov = table;

    if U.DoSilhouettePlot
        mvPerUnit = string(mvPerUnit(:));
        [clust,~,cl_idx] = unique(mvPerUnit, 'stable');

        counts = accumarray(cl_idx, 1);
        okSil = (numel(clust) >= 2) && (sum(counts >= 2) >= 2) && (size([PC1 PC2 PC3],1) >= 4);

        if ~okSil
            silVals = nan(size([PC1 PC2 PC3],1),1);
        else
            try
                silVals = silhouette([PC1 PC2 PC3], mvPerUnit, U.SilhouetteMethod);
            catch
                silVals = nan(size([PC1 PC2 PC3],1),1);
                okSil = false;
            end
        end

        if okSil
            silTbl = table(mvPerUnit, silVals, 'VariableNames', {'MoveType','Silhouette'});

            axes(axS); cla(axS);
            silhouette([PC1 PC2 PC3], mvPerUnit, U.SilhouetteMethod);
            if isKey(depthNames, dz)
                title(axS, sprintf('Silhouette | %s (MUA)', depthNames(dz)), 'Interpreter','none');
            else
                title(axS, sprintf('Silhouette | depth %s (MUA)', dz), 'Interpreter','none');
            end
            xlabel(axS, 'Silhouette value'); ylabel(axS, 'Cluster (MoveType)');

            % Summary rows
            for mIdx = 1:numel(mvOrder)
                mv = mvOrder(mIdx);
                mm = mvPerUnit == mv;
                if ~any(mm), continue; end
                silSummaryRows(end+1,:) = { ...
                    string(depthNames(dz)), string(mv), nnz(mm), ...
                    mean(silVals(mm),'omitnan'), ...
                    median(silVals(mm),'omitnan'), ...
                    std(silVals(mm),'omitnan')};
            end

            % ---- Stats tables (KW + ANOVA) ----
            factors = categorical(mvPerUnit, mvOrder, 'Ordinal', false);

            [p_kw, tbl_kw_raw, stats_kw] = kruskalwallis(silVals, factors, 'off');
            tblKW = cell2table(tbl_kw_raw(2:end,:), ...
                'VariableNames', matlab.lang.makeValidName(string(tbl_kw_raw(1,:))));
            tblKW.Depth = repmat(string(depthNames(dz)), height(tblKW), 1);
            tblKW = movevars(tblKW, 'Depth', 'Before', 1);
            tblKW.pValue_overall = repmat(p_kw, height(tblKW), 1);

            mc_kw = multcompare(stats_kw, 'CType','bonferroni', 'Display','off');
            mcTbl_kw = array2table(mc_kw, 'VariableNames', {'Group1','Group2','LowerCI','Diff','UpperCI','pValue'});
            grpNames = categories(factors);
            mcTbl_kw.Group1 = grpNames(mcTbl_kw.Group1);
            mcTbl_kw.Group2 = grpNames(mcTbl_kw.Group2);
            mcTbl_kw.Depth  = repmat(string(depthNames(dz)), height(mcTbl_kw), 1);
            mcTbl_kw = movevars(mcTbl_kw, 'Depth', 'Before', 1);

            [p_aov, tbl_aov_raw, stats_aov] = anova1(silVals, factors, 'off');
            tblAOV = cell2table(tbl_aov_raw(2:end,:), ...
                'VariableNames', matlab.lang.makeValidName(string(tbl_aov_raw(1,:))));
            tblAOV.Depth = repmat(string(depthNames(dz)), height(tblAOV), 1);
            tblAOV = movevars(tblAOV, 'Depth', 'Before', 1);
            tblAOV.pValue_overall = repmat(p_aov, height(tblAOV), 1);

            mc_aov = multcompare(stats_aov, 'CType','bonferroni', 'Display','off');
            mcTbl_aov = array2table(mc_aov, 'VariableNames', {'Group1','Group2','LowerCI','Diff','UpperCI','pValue'});
            mcTbl_aov.Group1 = grpNames(mcTbl_aov.Group1);
            mcTbl_aov.Group2 = grpNames(mcTbl_aov.Group2);
            mcTbl_aov.Depth  = repmat(string(depthNames(dz)), height(mcTbl_aov), 1);
            mcTbl_aov = movevars(mcTbl_aov, 'Depth', 'Before', 1);

            % Optional: extra histogram figures (overlay + stacked) like your SU
            if U.DoSilhouetteHistFigs
                makeSilhouetteHistFigures(silVals, mvPerUnit, mvOrder, moveColorMap, depthNames(dz));
            end
        else
            if ~isempty(axS)
                title(axS, sprintf('Silhouette not computed (insufficient groups) | %s (MUA)', string(depthNames(dz))), 'Interpreter','none');
            end
        end
    end

    %% ---- Save if requested ----
   
        base = sprintf('%s_%s', char(string(U.FigPrefix)), dz);
        print(hFig, fullfile(U.SaveDir, [base '.png']), '-dpng', '-r300');
        savefig(hFig, fullfile(U.SaveDir, [base '.fig']));

    % ---- Store per-depth outputs ----
    out.PerDepth(end+1).Depth      = dz;
    out.PerDepth(end).Coeff        = coeff;
    out.PerDepth(end).ExplainedPct = explained;
    out.PerDepth(end).PCA_Scores   = Tscore;
    out.PerDepth(end).Silhouette   = silTbl;

    out.PerDepth(end).KW.tbl    = tblKW;
    out.PerDepth(end).KW.mc     = mcTbl_kw;
    out.PerDepth(end).ANOVA.tbl = tblAOV;
    out.PerDepth(end).ANOVA.mc  = mcTbl_aov;

    % ---- Aggregate stats ----
    out.Stats.KW.tbl     = [out.Stats.KW.tbl;     tblKW];
    out.Stats.KW.mc      = [out.Stats.KW.mc;      mcTbl_kw];
    out.Stats.ANOVA.tbl  = [out.Stats.ANOVA.tbl;  tblAOV];
    out.Stats.ANOVA.mc   = [out.Stats.ANOVA.mc;   mcTbl_aov];
end

% ---- Summary table ----
if ~isempty(silSummaryRows)
    out.SilSummaryTbl = cell2table(silSummaryRows, ...
        'VariableNames', {'Depth','MoveType','N','SilMean','SilMedian','SilStd'});
else
    out.SilSummaryTbl = table;
end

%% ============================
% Save Output Structs and Tables
% =============================

% Ensure SaveDir exists
if ~exist(U.SaveDir,'dir')
    mkdir(U.SaveDir);
end

% Use a stable base name for all exports
baseOut = sprintf('%s_%s', char(string(U.FigPrefix)));

% ---- Save the full output struct ----
outMatPath = fullfile(U.SaveDir, [baseOut '_OUT.mat']);
save(outMatPath, 'out', '-v7.3');  % -v7.3 safe for large structs
fprintf('[SAVE] out struct: %s\n', outMatPath);

% ---- Save aggregated stats tables (CSV) ----

if istable(out.Stats.KW.tbl) && ~isempty(out.Stats.KW.tbl)
    writetable(out.Stats.KW.tbl, fullfile(U.SaveDir, [baseOut '_KW_tbl.csv']));
end
if istable(out.Stats.KW.mc) && ~isempty(out.Stats.KW.mc)
    writetable(out.Stats.KW.mc,  fullfile(U.SaveDir, [baseOut '_KW_mc.csv']));
end
if istable(out.Stats.ANOVA.tbl) && ~isempty(out.Stats.ANOVA.tbl)
    writetable(out.Stats.ANOVA.tbl, fullfile(U.SaveDir, [baseOut '_ANOVA_tbl.csv']));
end
if istable(out.Stats.ANOVA.mc) && ~isempty(out.Stats.ANOVA.mc)
    writetable(out.Stats.ANOVA.mc,  fullfile(U.SaveDir, [baseOut '_ANOVA_mc.csv']));
end
if istable(out.SilSummaryTbl) && ~isempty(out.SilSummaryTbl)
    writetable(out.SilSummaryTbl, fullfile(U.SaveDir, [baseOut '_SilSummaryTbl.csv']));
end
fprintf('[SAVE] CSV tables written to: %s\n', U.SaveDir);

end % main function


%% ============================
% Helpers
% ============================

function [X, tRef, mvPerUnit] = build_timeAlignedMatrix_MUA(Md, interpMethod)
% Build X: time x units from vecD_MUA interpolated to shared timebase new_vecTime_MUA.

tRef0 = Md.new_vecTime_MUA{find(~cellfun(@isempty, Md.new_vecTime_MUA), 1, 'first')};
tRef0 = double(tRef0(:));
tRef0 = tRef0(isfinite(tRef0));
if numel(tRef0) < 6
    error('Reference time vector too short / invalid (MUA).');
end

[tRef, ~] = unique(tRef0, 'stable');
tRef = tRef(:);

nU = height(Md);
X  = nan(numel(tRef), nU);
mvPerUnit = string(Md.MoveType);

for u = 1:nU
    dVec = Md.vecD_MUA{u};
    tVec = Md.new_vecTime_MUA{u};

    if isempty(dVec) || isempty(tVec), continue; end
    dVec = double(dVec(:));
    tVec = double(tVec(:));

    good = isfinite(dVec) & isfinite(tVec);
    dVec = dVec(good);
    tVec = tVec(good);

    if numel(tVec) < 6 || numel(dVec) ~= numel(tVec)
        continue;
    end

    [tU, iu] = unique(tVec, 'stable');
    dU = dVec(iu);

    if numel(tU) < 6, continue; end

    try
        X(:,u) = interp1(tU, dU, tRef, interpMethod, 'extrap');
    catch
        % leave NaNs
    end
end

% Drop units with all NaN
keepU = any(isfinite(X), 1);
X = X(:, keepU);
mvPerUnit = mvPerUnit(keepU);

% Fill remaining NaNs per unit with that unit's mean
for j = 1:size(X,2)
    nanMask = ~isfinite(X(:,j));
    if any(nanMask)
        mu = mean(X(~nanMask,j), 'omitnan');
        X(nanMask,j) = mu;
    end
end

% Drop timepoints that are NaN across all units
goodT = any(isfinite(X), 2);
X = X(goodT,:);
tRef = tRef(goodT);

end


function [wallInfo, limInfo] = computeWalls(x,y,z)
% Mirrors your SU "good visibility" wall placement.

xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
zmin = min(z); zmax = max(z);

mx = 0.10 * max(1e-6, xmax-xmin);
my = 0.10 * max(1e-6, ymax-ymin);
mz = 0.10 * max(1e-6, zmax-zmin);

% --- Walls ---
x_wall = xmax + (0.2+mx); % YZ projection on BACK-right plane
y_wall = ymax + (0.2+my); % XZ projection on back-left plane
z_wall = zmin - (0.2+mz); % XY projection on bottom plane

wallInfo = [x_wall, y_wall, z_wall];

% --- Axis limits ---
x_lim_min = xmin - (0.2+mx);
y_lim_min = ymin - (0.2+my);
z_lim_max = zmax + (0.2+mz);

limInfo = [ ...
    x_lim_min, x_wall;   % xlim
    y_lim_min, y_wall;   % ylim
    z_wall,    z_lim_max]; % zlim
end


function plotHullProjections(ax, x, y, z, c, wallInfo, faceAlpha, edgeAlpha, lw)
% Plot convex hull projections onto XY (z=z_wall), XZ (y=y_wall), YZ (x=x_wall)

x_wall = wallInfo(1);
y_wall = wallInfo(2);
z_wall = wallInfo(3);

x = x(:); y = y(:); z = z(:);

try
    if numel(x) >= 3
        k_xy = convhull(x, y);
        z_proj = repmat(z_wall, size(k_xy));
        p1 = patch(ax, x(k_xy), y(k_xy), z_proj, c, ...
            'FaceAlpha', faceAlpha, 'EdgeColor', c, 'LineWidth', lw);
        p1.EdgeAlpha = edgeAlpha;
        p1.FaceColor = 'none';
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


function makeSilhouetteHistFigures(silVals, mvPerUnit, mvOrder, moveColorMap, depthLabel)
% Optional helper: make the same overlay + stacked silhouette histograms you already use.

mvPerUnit = string(mvPerUnit(:));
silVals   = silVals(:);

% Overlay
figure('Color','w'); axH1 = axes; hold(axH1,'on'); box(axH1,'on'); grid(axH1,'on');
edges = linspace(-1, 1, 16);
h = gobjects(numel(mvOrder),1);
for mIdx = 1:numel(mvOrder)
    mv = mvOrder(mIdx);
    idx = (mvPerUnit == mv) & isfinite(silVals);
    if ~any(idx), continue; end
    c = moveColorMap(char(mv));
    h(mIdx) = histogram(axH1, silVals(idx), 'BinEdges', edges, ...
        'DisplayStyle','bar', 'FaceColor', c, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
end
xlim(axH1, [-1 1]); xlabel(axH1,'Silhouette value'); ylabel(axH1,'Count');
title(axH1, sprintf('Silhouette distribution by MoveType | %s (MUA)', string(depthLabel)), 'Interpreter','none');
legend(axH1, cellstr(mvOrder), 'Location','northeast');
xline(axH1, median(silVals(isfinite(silVals))), '-k', 'median', 'LineWidth', 1.25);

% Stacked
figure('Color','w'); axH2 = axes; hold(axH2,'on'); box(axH2,'on'); grid(axH2,'on');
centers = edges(1:end-1) + diff(edges)/2;
counts = zeros(numel(centers), numel(mvOrder));
for mIdx = 1:numel(mvOrder)
    mv = mvOrder(mIdx);
    idx = (mvPerUnit == mv) & isfinite(silVals);
    counts(:,mIdx) = histcounts(silVals(idx), edges).';
end
b = bar(axH2, centers, counts, 'stacked');
for mIdx = 1:numel(mvOrder)
    c = moveColorMap(char(mvOrder(mIdx)));
    b(mIdx).FaceColor = c;
    b(mIdx).EdgeColor = 'none';
end
xlim(axH2, [-1 1]); xlabel(axH2,'Silhouette value'); ylabel(axH2,'Count');
title(axH2, sprintf('Silhouette (stacked) by MoveType | %s (MUA)', string(depthLabel)), 'Interpreter','none');
legend(axH2, cellstr(mvOrder), 'Location','northeast');

end
