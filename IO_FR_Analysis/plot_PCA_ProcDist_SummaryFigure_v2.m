function out = plot_PCA_ProcDist_SummaryFigure_v2(DistTbl, varargin)
% plot_PCA_ProcDist_SummaryFigure_v2
%
% Visualize Procrustes distances (bootstrapped) for each movement-pair,
% faceted by depth (rows). Adds per-depth ANOVA + Bonferroni posthoc
% significance bars.
%
% Expected DistTbl columns:
%   Depth   (string/categorical/char)  e.g., "dorsal STN"
%   MoveA   (string/categorical/char)  e.g., "HAND OC"
%   MoveB   (string/categorical/char)  e.g., "ARM EF"
%   BootIter (numeric)                e.g., 1..Nboot
%   ProcDist (numeric)                Procrustes distance
%
% Example:
%   plot_PCA_ProcDist_SummaryFigure_v2(out_SUtraj.PerDepth(1).ProcDist, ...)
%   plot_PCA_ProcDist_SummaryFigure_v2(DistTbl_allDepths, ...)

% -----------------------
% Inputs
% -----------------------
p = inputParser;
p.addParameter('SaveDir', '', @(x)ischar(x)||isstring(x));
p.addParameter('FigName', 'ProcDist_Summary', @(x)ischar(x)||isstring(x));
p.addParameter('DoSave', false, @islogical);

p.addParameter('DepthOrder', ["dorsal STN","central STN","ventral STN"], @(x)isstring(x)||iscellstr(x));
p.addParameter('MoveOrder',  ["HAND OC","HAND PS","ARM EF"], @(x)isstring(x)||iscellstr(x));

% aesthetics
p.addParameter('PointAlpha', 0.18, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('BoxAlpha',   0.35, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('JitterWidth', 0.18, @(x)isnumeric(x)&&isscalar(x)&&x>0);

p.parse(varargin{:});
U = p.Results;

if U.DoSave && isempty(U.SaveDir)
    U.SaveDir = fullfile(pwd, 'TrajectoryDistance_Summary');
end
if U.DoSave && ~exist(U.SaveDir,'dir')
    mkdir(U.SaveDir);
end

% -----------------------
% Validate + normalize
% -----------------------
req = {'Depth','MoveA','MoveB','ProcDist'};
missing = setdiff(req, DistTbl.Properties.VariableNames);
if ~isempty(missing)
    error('DistTbl missing required variables: %s', strjoin(missing, ', '));
end

T = DistTbl;

T.Depth   = string(T.Depth);
T.MoveA   = string(T.MoveA);
T.MoveB   = string(T.MoveB);
T.ProcDist = double(T.ProcDist);

% Build Pair label
% Use a consistent ordering for label display (MoveOrder-driven)
moveOrder = string(U.MoveOrder);

% helper: sort the pair names in a stable "MoveOrder" sense
ordIdxA = arrayfun(@(s) find(moveOrder==s,1,'first'), T.MoveA);
ordIdxB = arrayfun(@(s) find(moveOrder==s,1,'first'), T.MoveB);

swap = ordIdxA > ordIdxB;
A = T.MoveA; B = T.MoveB;
A(swap) = T.MoveB(swap);
B(swap) = T.MoveA(swap);

T.Pair = A + " vs " + B;

% Define canonical pair order: (OC vs PS), (OC vs EF), (PS vs EF)
pairOrder = [
    moveOrder(1) + " vs " + moveOrder(2)
    moveOrder(1) + " vs " + moveOrder(3)
    moveOrder(2) + " vs " + moveOrder(3)
];
T.Pair = categorical(T.Pair, pairOrder, 'Ordinal', true);

% drop NaNs
T = T(isfinite(T.ProcDist), :);

% -----------------------
% Color scheme (matches your SU palette)
% (Pairs get distinct, readable colors; tweak if you want)
% -----------------------
pairColors = [
    0.35 0.75 0.78  % OC vs PS (teal)
    0.70 0.55 0.75  % OC vs EF (purple-ish)
    0.75 0.72 0.45  % PS vs EF (olive)
];

% -----------------------
% Figure
% -----------------------
depthOrder = string(U.DepthOrder);

hFig = figure('Color','w','Units','pixels','Position',[80 80 1600 850]);
tlo = tiledlayout(hFig, numel(depthOrder), 1, 'TileSpacing','compact', 'Padding','compact');

out = struct;
out.PerDepth = struct([]);

for di = 1:numel(depthOrder)
    dz = depthOrder(di);
    ax = nexttile(tlo, di);
    hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    Td = T(T.Depth == dz, :);

    % If DistTbl passed was already single-depth, try graceful fallback:
    if isempty(Td) && di==1 && numel(unique(T.Depth))==1
        dz = unique(T.Depth);
        Td = T;
    end

    if isempty(Td)
        title(ax, sprintf('Procrustes distance | %s (no data)', dz), 'Interpreter','none');
        axis(ax,'off');
        continue;
    end

    % Keep only valid pair groups present
    Td = Td(~isundefined(Td.Pair), :);

    presentPairs = categories(removecats(Td.Pair));
    nGroups = numel(presentPairs);

    if nGroups < 2
        title(ax, sprintf('Procrustes distance | %s (insufficient groups)', dz), 'Interpreter','none');
        continue;
    end

    % --- Boxchart for each pair ---
    % Convert to numeric x positions 1..3 for stable layout
    pairCats = categories(T.Pair);
    [~,xpos] = ismember(string(Td.Pair), string(pairCats));
    xpos = xpos(:);

    % Draw boxcharts group-wise so we can color them
    for gi = 1:numel(pairOrder)
        thisPair = pairOrder(gi);
        mask = (string(Td.Pair) == thisPair);
        if ~any(mask), continue; end

        bc = boxchart(ax, repmat(gi, sum(mask), 1), Td.ProcDist(mask), ...
            'BoxFaceColor', pairColors(gi,:), ...
            'BoxFaceAlpha', U.BoxAlpha, ...
            'MarkerStyle','none', ...
            'WhiskerLineColor',[0.2 0.2 0.2], ...
            'LineWidth',1.1);

        % jittered points (bootstraps)
        xjit = gi + (rand(sum(mask),1)-0.5)*U.JitterWidth;
        sc = scatter(ax, xjit, Td.ProcDist(mask), 10, pairColors(gi,:), 'filled');
        sc.MarkerFaceAlpha = U.PointAlpha;
        sc.MarkerEdgeAlpha = 0;
    end

    ax.XLim = [0.5 3.5];
    ax.XTick = 1:3;
    ax.XTickLabel = cellstr(pairOrder);
    ax.XTickLabelRotation = 18;
    ylabel(ax, 'Procrustes distance (unitless)');
    title(ax, sprintf('Trajectory Dissimilarity Pairwise Comparisons | %s', dz), 'Interpreter','none');

    % --- Stats: one-way ANOVA across pairs (within this depth) ---
    grp = categorical(string(Td.Pair), pairOrder, 'Ordinal', true);
    y = Td.ProcDist;

    % Must have >=2 groups with >=2 samples each for multcompare
    counts = countcats(grp);
    okAOV = (sum(counts > 1) >= 2);

    statsBlock = struct('p_anova', nan, 'tbl', [], 'mc', []);
    if okAOV
        [p_aov, tbl_aov, stats_aov] = anova1(y, grp, 'off');

        statsBlock.p_anova = p_aov;

        % multcompare only if >=2 groups
        try
            mc = multcompare(stats_aov, 'CType','bonferroni', 'Display','off');
            % mc columns: [g1 g2 low diff up p]
            mcTbl = array2table(mc, 'VariableNames', {'G1','G2','LowerCI','Diff','UpperCI','pValue'});
            gnames = stats_aov.gnames;
            mcTbl.Group1 = string(gnames(mcTbl.G1));
            mcTbl.Group2 = string(gnames(mcTbl.G2));
            mcTbl = removevars(mcTbl, {'G1','G2'});

            statsBlock.tbl = tbl_aov;
            statsBlock.mc  = mcTbl;

            % Annotate significance bars on plot
            yTop = max(y) * 1.08;
            step = max(y) * 0.06;
            if ~isfinite(yTop) || yTop==0
                yTop = 1;
                step = 0.1;
            end

            % map group names -> x positions 1..3
            xMap = containers.Map(cellstr(pairOrder), num2cell(1:3));
            drawCount = 0;

            % Only draw significant comparisons
            sigRows = mcTbl(mcTbl.pValue < 0.05, :);
            % sort by p (smallest first) for clean stacking
            sigRows = sortrows(sigRows, 'pValue', 'ascend');

            for r = 1:height(sigRows)
                g1 = sigRows.Group1(r);
                g2 = sigRows.Group2(r);
                if ~isKey(xMap, char(g1)) || ~isKey(xMap, char(g2)), continue; end

                x1 = xMap(char(g1));
                x2 = xMap(char(g2));
                if x1==x2, continue; end
                drawCount = drawCount + 1;

                yy = yTop + (drawCount-1)*step;
                plot(ax, [x1 x1 x2 x2], [yy yy+step*0.15 yy+step*0.15 yy], ...
                    'k-', 'LineWidth', 1.0);

                % stars
                pval = sigRows.pValue(r);
                stars = p_to_stars(pval);
                text(ax, mean([x1 x2]), yy+step*0.18, stars, ...
                    'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                    'FontSize', 11, 'FontWeight','bold');
            end

        catch ME
            warning('Posthoc multcompare failed for %s: %s', dz, ME.message);
        end

        % report ANOVA p on panel
        txt = sprintf('ANOVA p = %.3g', p_aov);
        text(ax, 0.99, 0.93, txt, 'Units','normalized', ...
            'HorizontalAlignment','right', 'FontSize', 10);
    else
        text(ax, 0.99, 0.93, 'ANOVA not run (insufficient groups)', ...
            'Units','normalized', 'HorizontalAlignment','right', 'FontSize', 10);
    end

    out.PerDepth(end+1).Depth = dz;
    out.PerDepth(end).Stats  = statsBlock;
end

sgtitle(hFig, U.FigName, 'FontWeight','bold', 'Interpreter','none');

if U.DoSave
    pngPath = fullfile(U.SaveDir, U.FigName + ".png");
    figPath = fullfile(U.SaveDir, U.FigName + ".fig");
    exportgraphics(hFig, pngPath, 'Resolution', 300);
    savefig(hFig, figPath);
    fprintf('[SAVE] %s\n', pngPath);
    fprintf('[SAVE] %s\n', figPath);
end

end

% -----------------------
% helper: p-value to stars
% -----------------------
function s = p_to_stars(p)
if p < 0.001
    s = '***';
elseif p < 0.01
    s = '**';
elseif p < 0.05
    s = '*';
else
    s = 'n.s.';
end
end
