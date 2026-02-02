function out = plot_PCA_ProcDist_SummaryFigure(DistTbl, varargin)
% plot_PCA_ProcDist_SummaryFigure
% Figure 2: Procrustes distance summary
% - Boxplots of ProcDist per Pair, faceted by Depth (rows)
% - one-way ANOVA + Bonferroni post-hoc within each depth (ProcDist ~ Pair)
% - Significance brackets: * p<0.05, ** p<0.01, *** p<0.001
%
% Expected DistTbl columns:
%   Depth (string/categorical)
%   MoveA (string/categorical)
%   MoveB (string/categorical)
%   BootIter (numeric) [optional]
%   ProcDist (numeric)
%
% Example:
%   plot_PCA_ProcDist_SummaryFigure(DistTbl, 'SaveDir', SaveDir, ...
%       'FigName','SU_ProcrustesDistance_Summary','DoSave',true);

% -----------------------
% Parse inputs
% -----------------------
p = inputParser;
p.addParameter('SaveDir', '', @(x)ischar(x)||isstring(x));
p.addParameter('FigName', 'ProcDist_Summary', @(x)ischar(x)||isstring(x));
p.addParameter('DoSave', false, @islogical);

% Plot options
p.addParameter('ShowPoints', true, @islogical);
p.addParameter('PointAlpha', 0.12, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('BoxWidth', 0.55, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('YLim', [], @(x)isnumeric(x)&&isempty(x) || (numel(x)==2));

% Depth order (match your labels)
p.addParameter('DepthOrder', ["dorsal STN","central STN","ventral STN"], @(x)isstring(x)||iscell(x));
p.parse(varargin{:});
U = p.Results;

% -----------------------
% Sanity checks
% -----------------------
req = {'Depth','MoveA','MoveB','ProcDist'};
missing = setdiff(req, DistTbl.Properties.VariableNames);
if ~isempty(missing)
    error('DistTbl missing required columns: %s', strjoin(missing,', '));
end

T = DistTbl;
T.Depth   = string(T.Depth);
T.MoveA   = string(T.MoveA);
T.MoveB   = string(T.MoveB);
T.ProcDist = double(T.ProcDist);

% Drop non-finite
T = T(isfinite(T.ProcDist), :);

% Build Pair label
T.Pair = T.MoveA + " vs " + T.MoveB;

% Enforce depth order
depthOrder = string(U.DepthOrder(:));
T.Depth = categorical(T.Depth, depthOrder, 'Ordinal', true);

% Pair order (stable, but keep consistent across depths)
pairOrder = unique(T.Pair, 'stable');
T.Pair = categorical(T.Pair, pairOrder, 'Ordinal', false);

% -----------------------
% Colors (JNE_move scheme)
% -----------------------
JNE_move = ([ ...
    128,128,128;    % REST  (grey)
    38,116,183;     % HAND OC  (blue)
    53,183,121;     % HAND PS  (green/teal)
    243,120,98] ... % ARM EF   (coral)
    ./ 255);

moveColorMap = containers.Map( ...
    {'HAND OC','HAND PS','ARM EF'}, ...
    {JNE_move(2,:), JNE_move(3,:), JNE_move(4,:)} );

% Pair color = mean of the two move colors (works nicely for "distance between two moves")
pairCats = categories(T.Pair);
pairColors = zeros(numel(pairCats),3);
for i = 1:numel(pairCats)
    s = string(pairCats{i});
    toks = split(s, " vs ");
    if numel(toks) == 2 && isKey(moveColorMap, char(toks(1))) && isKey(moveColorMap, char(toks(2)))
        pairColors(i,:) = 0.5*(moveColorMap(char(toks(1))) + moveColorMap(char(toks(2))));
    else
        pairColors(i,:) = [0.6 0.6 0.6];
    end
end

% -----------------------
% Figure
% -----------------------
hFig = figure('Color','w','Units','pixels','Position',[120 80 900 950]);
tlo = tiledlayout(hFig, numel(depthOrder), 1, 'TileSpacing','compact','Padding','compact');

out = struct();
out.ANOVA = table();
out.MC    = table();

for d = 1:numel(depthOrder)
    dz = depthOrder(d);
    ax = nexttile(tlo, d);
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    Td = T(T.Depth == dz, :);

    title(ax, sprintf('Procrustes distance by Pair | %s', string(dz)), 'Interpreter','none');
    ylabel(ax, 'Procrustes distance');
    if d == numel(depthOrder)
        xlabel(ax, 'Movement Pair');
    else
        set(ax, 'XTickLabel', []);
    end

    if isempty(Td)
        text(ax, 0.5, 0.5, 'No data', 'Units','normalized', 'HorizontalAlignment','center');
        continue;
    end

    % Boxcharts by Pair
    cats = categories(Td.Pair);
    nC = numel(cats);
    xPos = 1:nC;

    for i = 1:nC
        thisCat = cats{i};
        yi = Td.ProcDist(Td.Pair == thisCat);

        if isempty(yi), continue; end

        % Boxchart expects a grouping; easiest is numeric x positions
        boxchart(ax, repmat(i, size(yi)), yi, ...
            'BoxWidth', U.BoxWidth, ...
            'MarkerStyle', 'none', ...
            'BoxFaceColor', pairColors(find(strcmp(pairCats, thisCat),1,'first'),:), ...
            'WhiskerLineColor', [0.2 0.2 0.2], ...
            'BoxFaceAlpha', 0.55);

        % Jittered points (bootstrap iterations)
        if U.ShowPoints
            xjit = i + 0.12*(rand(size(yi))-0.5);
            sc = scatter(ax, xjit, yi, 10, 'filled');
            sc.MarkerFaceColor = pairColors(find(strcmp(pairCats, thisCat),1,'first'),:);
            sc.MarkerFaceAlpha = U.PointAlpha;
            sc.MarkerEdgeColor = 'none';
        end
    end

    set(ax, 'XTick', xPos, 'XTickLabel', cats, 'XTickLabelRotation', 25);

    % Y limits (optional)
    if ~isempty(U.YLim)
        ylim(ax, U.YLim);
    else
        yl = ylim(ax);
        ylim(ax, [max(0, yl(1)), yl(2)]);
    end

    % -----------------------
    % Stats within depth: ANOVA + Bonferroni post-hoc across Pair
    % -----------------------
    pairGrp = Td.Pair;
    if numel(categories(pairGrp)) < 2
        % not enough groups to test
        continue;
    end

    try
        [p_aov, tbl_aov, stats_aov] = anova1(Td.ProcDist, pairGrp, 'off');

        tblAOV = cell2table(tbl_aov(2:end,:), ...
            'VariableNames', matlab.lang.makeValidName(string(tbl_aov(1,:))));
        tblAOV.Depth = repmat(string(dz), height(tblAOV), 1);
        tblAOV = movevars(tblAOV, 'Depth', 'Before', 1);
        tblAOV.pValue_overall = repmat(p_aov, height(tblAOV), 1);

        out.ANOVA = [out.ANOVA; tblAOV];

        mc = multcompare(stats_aov, 'CType','bonferroni', 'Display','off');
        mcTbl = array2table(mc, 'VariableNames', {'Group1','Group2','LowerCI','Diff','UpperCI','pValue'});
        grpNames = categories(pairGrp);
        mcTbl.Group1 = grpNames(mcTbl.Group1);
        mcTbl.Group2 = grpNames(mcTbl.Group2);
        mcTbl.Depth  = repmat(string(dz), height(mcTbl), 1);
        mcTbl = movevars(mcTbl, 'Depth', 'Before', 1);

        out.MC = [out.MC; mcTbl];

        % Add significance brackets to this axis
        addSigBrackets(ax, mcTbl, cats);

        % Add overall p-value to corner
        txt = sprintf('ANOVA p = %.2g', p_aov);
        text(ax, 0.98, 0.95, txt, 'Units','normalized', ...
            'HorizontalAlignment','right', 'VerticalAlignment','top', ...
            'FontSize', 10, 'Color', [0.15 0.15 0.15]);

    catch ME
        warning('Stats failed for %s: %s', string(dz), ME.message);
    end
end

sgtitle(tlo, char(string(U.FigName)), 'Interpreter','none', 'FontWeight','bold');

% -----------------------
% Save
% -----------------------
if U.DoSave
    if ~isempty(U.SaveDir) && ~exist(U.SaveDir,'dir')
        mkdir(U.SaveDir);
    end
    if isempty(U.SaveDir), U.SaveDir = pwd; end

    outPng = fullfile(U.SaveDir, char(string(U.FigName) + ".png"));
    outFig = fullfile(U.SaveDir, char(string(U.FigName) + ".fig"));
    exportgraphics(hFig, outPng, 'Resolution', 300);
    savefig(hFig, outFig);
    fprintf('[SAVE] %s\n[SAVE] %s\n', outPng, outFig);
end

end % main


% =====================================================================
% Helper: add significance brackets for pairwise comparisons
% =====================================================================
function addSigBrackets(ax, mcTbl, cats)
% mcTbl has Group1/Group2 as strings and pValue
% cats is cell array of category labels in x order

% only significant comparisons
mcTbl = mcTbl(isfinite(mcTbl.pValue), :);
mcTbl = sortrows(mcTbl, 'pValue', 'ascend');

if isempty(mcTbl), return; end

% Map category -> x position
xMap = containers.Map();
for i = 1:numel(cats)
    xMap(char(string(cats{i}))) = i;
end

% Y baseline
yl = ylim(ax);
y0 = yl(2);
dy = 0.06 * range(yl);
y = y0 - 0.02*range(yl);  % start a hair below top
nDrawn = 0;

for r = 1:height(mcTbl)
    p = mcTbl.pValue(r);
    if p >= 0.05
        continue;
    end

    g1 = char(string(mcTbl.Group1(r)));
    g2 = char(string(mcTbl.Group2(r)));
    if ~isKey(xMap,g1) || ~isKey(xMap,g2), continue; end

    x1 = xMap(g1);
    x2 = xMap(g2);
    if x1 == x2, continue; end
    if x2 < x1
        tmp = x1; x1 = x2; x2 = tmp;
    end

    % Limit how many brackets to avoid clutter (top 3 smallest p-values)
    nDrawn = nDrawn + 1;
    if nDrawn > 3, break; end

    y = y + dy;

    % bracket
    plot(ax, [x1 x1 x2 x2], [y-dy*0.2 y y y-dy*0.2], 'k-', 'LineWidth', 1.0);

    % stars
    if p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    else
        stars = '*';
    end
    text(ax, mean([x1 x2]), y + dy*0.05, stars, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontSize', 11, 'Color', 'k', 'FontWeight','bold');
end

% expand ylim if needed to show brackets
yl2 = ylim(ax);
if y > yl2(2)
    ylim(ax, [yl2(1) y + dy*0.6]);
end
end
