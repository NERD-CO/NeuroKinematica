function out = plot_ZETA_TrajectoryAnalysis_SU(MasterZETA, varargin)
% plot_ZETA_TrajectoryAnalysis_SU
%
% PURPOSE (Trajectory pipeline):
%   For each depth (t/c/b) and MoveType (HAND OC/PS/ARM EF),
%   compute a mean ZETA trajectory across units, convert it into a 3D
%   feature-trajectory over time, visualize it, and compute pairwise
%   trajectory distances using Procrustes (bootstrapped).
%
% KEY IDEA:
%   Each MoveType trajectory is a 3D curve over time:
%     f1(t)=zscore(mu(t))
%     f2(t)=zscore(dmu/dt)
%     f3(t)=zscore(cumtrapz(mu(t)))
%
% INPUT (required columns in MasterZETA):
%   - MoveType, Depth, ZETA_vecD, ZETA_vecT, PSTH_TimeCenters_s, PrettyLabel
%
% OUTPUT:
%   out.PerDepth(d):
%       .Depth
%       .TrajMean   table with time + 3D features per MoveType
%       .ProcDist   table of bootstrap Procrustes distances (pairwise MoveType)
%   out.ANOVA (optional aggregate across depths):
%       .tbl, .mc
%
% NOTE:
%   This is separate from your unit-PCA + silhouette pipeline.

% -----------------------
% Parse inputs
% -----------------------
p = inputParser;
p.addParameter('DepthOrder', {'t','c','b'}, @(x) iscell(x) && ~isempty(x));
p.addParameter('ActiveMoveOrder', {'HAND OC','HAND PS','ARM EF'}, @(x) iscell(x) && ~isempty(x));

p.addParameter('SaveDir', '', @(x)ischar(x)||isstring(x));
p.addParameter('FigPrefix', 'SU_ZETA_Traj', @(x)ischar(x)||isstring(x));

p.addParameter('Nboot', 500, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('BootMinUnits', 5, @(x)isnumeric(x)&&isscalar(x)&&x>=2);

% Trajectory feature options
p.addParameter('DoZscoreFeatures', true, @islogical);
p.addParameter('UseIntegralFeature', true, @islogical); % if false, uses 2nd derivative
p.addParameter('SmoothWin', 5, @(x)isnumeric(x)&&isscalar(x)&&x>=1); % movmean window on mu(t)

% Procrustes options
p.addParameter('DoProcrustes', true, @islogical);
p.addParameter('ProcrustesScaling', true, @islogical);
p.addParameter('ProcrustesReflection', false, @islogical);

% Plot
p.addParameter('View', [-40 20], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('LineWidth', 3, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('StartMarkerSize', 60, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('EndMarkerSize', 90, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('GradientN', 120, @(x)isnumeric(x)&&isscalar(x)&&x>=10); % downsample segments

p.parse(varargin{:});
U = p.Results;

if isempty(U.SaveDir)
    U.SaveDir = fullfile(pwd, 'Trajectory Stats');
end
if ~exist(U.SaveDir,'dir')
    mkdir(U.SaveDir);
end

% -----------------------
% Sanity checks
% -----------------------
reqVars = {'MoveType','Depth','ZETA_vecD','ZETA_vecT','PSTH_TimeCenters_s','PrettyLabel'};
missing = setdiff(reqVars, MasterZETA.Properties.VariableNames);
if ~isempty(missing)
    error('MasterZETA missing required variables: %s', strjoin(missing, ', '));
end

% Normalize
MasterZeta_T = MasterZETA;
MasterZeta_T.MoveType = string(MasterZeta_T.MoveType);
MasterZeta_T.Depth    = string(MasterZeta_T.Depth);

% Restrict MoveTypes
mvOrder = string(U.ActiveMoveOrder);
MasterZeta_T = MasterZeta_T(ismember(MasterZeta_T.MoveType, mvOrder), :);

% -----------------------
% Color scheme (your new scheme)
% -----------------------
JNE_move = ([ ...
    128,128,128;    % REST (unused here)
    38,116,183;     % HAND OC (blue)
    53,183,121;     % HAND PS (green)
    243,120,98] ... % ARM EF  (orange)
    ./ 255);

moveColorMap = containers.Map( ...
    {'HAND OC','HAND PS','ARM EF'}, ...
    {JNE_move(2,:), JNE_move(3,:), JNE_move(4,:)} );

depthNames = containers.Map({'t','c','b'}, {'dorsal STN','central STN','ventral STN'});

% -----------------------
% Main loop
% -----------------------
depthOrder = string(U.DepthOrder);
out = struct;
out.PerDepth = struct([]);

% Figure 1: trajectories per depth (tiled)
hFig1 = figure('Color','w','Units','pixels','Position',[120 80 900 1100]);
tlo1  = tiledlayout(hFig1, numel(depthOrder), 1, 'TileSpacing','compact','Padding','compact');

allProcRows = {}; % [Depth, Pair, dProc] for optional ANOVA

for d = 1:numel(depthOrder)
    dz = char(depthOrder(d));
    Md = MasterZeta_T(MasterZeta_T.Depth == string(dz), :);
    if height(Md) < 3
        continue;
    end

    % Build time x units matrix for this depth
    [X, t_use, mvPerUnit, ~] = build_SU_ZETA_timeAlignedMatrix(Md, 'DepthCode', dz);
    if isempty(X) || size(X,2) < 3
        continue;
    end

    ax = nexttile(tlo1, d);
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    view(ax, U.View);
    axis(ax,'vis3d');
    xlabel(ax,'F1'); ylabel(ax,'F2'); zlabel(ax,'F3');

    ttl = sprintf('SU Trajectory Features | %s', depthNames(dz));
    title(ax, ttl, 'Interpreter','none');

    % Compute non-bootstrapped mean trajectory per MoveType for plotting
    trajMean = struct();

    for mIdx = 1:numel(mvOrder)
        mv = char(mvOrder(mIdx));
        cols = find(mvPerUnit == string(mv));
        if numel(cols) < U.BootMinUnits
            continue;
        end

        mu = mean(X(:, cols), 2, 'omitnan');               % mean [time × units] ZETA deviation values (ZETA_vecD interpolated)
        mu = movmean(mu, U.SmoothWin, 'omitnan');          % gentle smoothing

        P = buildTrajectoryFeatures(t_use, mu, U.UseIntegralFeature, U.DoZscoreFeatures); % [T x 3]
        trajMean.(matlab.lang.makeValidName(mv)) = P;

        % Plot: gradient along time + start/end markers
        c = moveColorMap(mv);
        plotTrajectoryGradient3D(ax, P, c, U.LineWidth, U.GradientN, ...
            'SmoothMode','arclength', 'SmoothMethod','pchip', 'UpFactor',10, ...
            'UseCSAPS',true,'CSAPS_p',0.30); % Decrease CSAPS_p to smooth more

        plot3(ax, P(1,1), P(1,2), P(1,3), 'o', 'MarkerSize', sqrt(U.StartMarkerSize), ...
            'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'LineWidth', 0.5);
        plot3(ax, P(end,1), P(end,2), P(end,3), '^', 'MarkerSize', sqrt(U.EndMarkerSize), ...
            'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'LineWidth', 0.5);
    end

    % Legend (opaque dummy)
    hLeg = gobjects(numel(mvOrder),1);
    for mIdx = 1:numel(mvOrder)
        mv = char(mvOrder(mIdx));
        if ~isKey(moveColorMap, mv), continue; end
        c = moveColorMap(mv);
        hLeg(mIdx) = plot3(ax, NaN,NaN,NaN, '-', 'Color', c, 'LineWidth', 4);
    end
    legend(ax, hLeg, cellstr(mvOrder), 'Location','northeastoutside');

    
    % -----------------------
    % Procrustes distances (bootstrapped)
    % -----------------------
    procTbl = table;

    if U.DoProcrustes && U.Nboot > 0
        procRows = {};
        pairList = nchoosek(1:numel(mvOrder), 2);

        for b = 1:U.Nboot
            % bootstrap mean trajectory for each MoveType
            Pboot = cell(numel(mvOrder),1);
            ok = true;

            for mIdx = 1:numel(mvOrder)
                mv = char(mvOrder(mIdx));
                cols = find(mvPerUnit == string(mv));
                if numel(cols) < U.BootMinUnits
                    ok = false; break;
                end

                % sample units WITH replacement
                sampCols = cols(randi(numel(cols), numel(cols), 1));
                mu_b = mean(X(:, sampCols), 2, 'omitnan');
                mu_b = movmean(mu_b, U.SmoothWin, 'omitnan');

                Pboot{mIdx} = buildTrajectoryFeatures(t_use, mu_b, U.UseIntegralFeature, U.DoZscoreFeatures);
            end
            if ~ok, continue; end

            % compute pairwise Procrustes
            for pIdx = 1:size(pairList,1)
                iA = pairList(pIdx,1);
                iB = pairList(pIdx,2);

                Sa = Pboot{iA};
                Sb = Pboot{iB};

                dProc = procrustesDistance(Sa, Sb, U.ProcrustesScaling, U.ProcrustesReflection);

                procRows(end+1,:) = {string(depthNames(dz)), string(mvOrder(iA)), string(mvOrder(iB)), b, dProc}; %#ok<AGROW>
                allProcRows(end+1,:) = {string(depthNames(dz)), sprintf('%s vs %s', mvOrder(iA), mvOrder(iB)), dProc}; %#ok<AGROW>
            end
        end

        if ~isempty(procRows)
            procTbl = cell2table(procRows, 'VariableNames', {'Depth','MoveA','MoveB','BootIter','ProcDist'});
        end
    end

    % Store outputs per depth
    out.PerDepth(end+1).Depth   = dz;
    out.PerDepth(end).DepthName = string(depthNames(dz));
    out.PerDepth(end).t_use     = t_use;
    out.PerDepth(end).TrajMean  = trajMean;
    out.PerDepth(end).ProcDist  = procTbl;
end

% -----------------------
% Optional: ANOVA across depths for each pair (separately)
% -----------------------
out.ANOVA = struct('tbl', table(), 'mc', table());

if ~isempty(allProcRows)
    Dall = cell2table(allProcRows, 'VariableNames', {'Depth','Pair','ProcDist'});

    pairs = unique(Dall.Pair);
    Dall.Pair = categorical(string(Dall.Pair));
    anovaRows = table();
    mcAll = table();

    for i = 1:numel(pairs)
        thisPair = pairs(i);
        thisPair = categorical(string(thisPair));

        % T = Dall(Dall.Pair == thisPair, :);
        T = Dall(strcmp(string(Dall.Pair), string(thisPair)), :);


        if numel(unique(T.Depth)) < 2
            continue;
        end

        [p_aov, tbl_aov, stats_aov] = anova1(T.ProcDist, categorical(T.Depth), 'off');

        tblAOV = cell2table(tbl_aov(2:end,:), ...
            'VariableNames', matlab.lang.makeValidName(string(tbl_aov(1,:))));
        tblAOV.Pair = repmat(thisPair, height(tblAOV), 1);
        tblAOV.pValue_overall = repmat(p_aov, height(tblAOV), 1);
        tblAOV = movevars(tblAOV, {'Pair','pValue_overall'}, 'Before', 1);

        anovaRows = [anovaRows; tblAOV];

        mc = multcompare(stats_aov, 'CType','bonferroni', 'Display','off');
        mcTbl = array2table(mc, 'VariableNames', {'Group1','Group2','LowerCI','Diff','UpperCI','pValue'});
        grpNames = categories(categorical(T.Depth));
        mcTbl.Group1 = grpNames(mcTbl.Group1);
        mcTbl.Group2 = grpNames(mcTbl.Group2);
        mcTbl.Pair = repmat(thisPair, height(mcTbl), 1);
        mcTbl = movevars(mcTbl, 'Pair', 'Before', 1);

        mcAll = [mcAll; mcTbl];
    end

    out.ANOVA.tbl = anovaRows;
    out.ANOVA.mc  = mcAll;
end

% -----------------------
% Save
% -----------------------
base = char(string(U.FigPrefix));
figPathPng = fullfile(U.SaveDir, [base '_Fig1_Trajectories.png']);
figPathFig = fullfile(U.SaveDir, [base '_Fig1_Trajectories.fig']);
print(hFig1, figPathPng, '-dpng', '-r300');
savefig(hFig1, figPathFig);

outMatPath = fullfile(U.SaveDir, [base '_OUT.mat']);
save(outMatPath, 'out', '-v7.3');

% Save tables
if ~isempty(out.PerDepth)
    for k = 1:numel(out.PerDepth)
        if istable(out.PerDepth(k).ProcDist) && ~isempty(out.PerDepth(k).ProcDist)
            writetable(out.PerDepth(k).ProcDist, fullfile(U.SaveDir, sprintf('%s_ProcDist_%s.csv', base, out.PerDepth(k).Depth)));
        end
    end
end
if istable(out.ANOVA.tbl) && ~isempty(out.ANOVA.tbl)
    writetable(out.ANOVA.tbl, fullfile(U.SaveDir, [base '_ANOVA_ProcDist_tbl.csv']));
end
if istable(out.ANOVA.mc) && ~isempty(out.ANOVA.mc)
    writetable(out.ANOVA.mc, fullfile(U.SaveDir, [base '_ANOVA_ProcDist_mc.csv']));
end

fprintf('[DONE] Saved outputs to: %s\n', U.SaveDir);

end


% =====================================================================
% Helpers
% =====================================================================

function P = buildTrajectoryFeatures(t, mu, useIntegral, doZ)
% Return [T x 3] feature trajectory from mean ZETA mu(t)
t = double(t(:));
mu = double(mu(:));

dt = median(diff(t));
if ~isfinite(dt) || dt <= 0
    dt = 1;
end

f1 = mu;
f2 = gradient(mu, dt);

if useIntegral
    f3 = cumtrapz(t, mu);
else
    f3 = gradient(f2, dt); % second derivative
end

P = [f1(:) f2(:) f3(:)];

if doZ
    P = zscore(P, 0, 1);
end

% safety: remove any remaining NaNs by column mean
for j = 1:size(P,2)
    bad = ~isfinite(P(:,j));
    if any(bad)
        P(bad,j) = mean(P(~bad,j), 'omitnan');
    end
end
end


function dProc = procrustesDistance(A, B, doScaling, doReflection)
% Proper Procrustes distance wrapper.
% NOTE: MATLAB returns d as first output (not in transform struct).
A = double(A); B = double(B);
% equalize length if needed by interpolation onto common index
nA = size(A,1); nB = size(B,1);
n  = max(nA, nB);
A2 = interp1(linspace(0,1,nA), A, linspace(0,1,n), 'linear');
B2 = interp1(linspace(0,1,nB), B, linspace(0,1,n), 'linear');

[d, ~, ~] = procrustes(A2, B2, 'Scaling', doScaling, 'Reflection', doReflection);
dProc = d;
end



function plotTrajectoryGradient3D(ax, P, baseColor, lw, nSeg, varargin)
% plotTrajectoryGradient3D
% Draw a time-gradient line by plotting short segments with increasing brightness.
% Now includes smoothing via interpolation (pchip/spline) and optional smoothing spline (csaps).
%
% Inputs
%   ax        : axes handle
%   P         : [T x 3] trajectory points
%   baseColor : [1 x 3] RGB
%   lw        : line width
%   nSeg      : number of segments used for gradient drawing (also controls downsampling)
%
% Name-Value options
%   'SmoothMethod'  : 'pchip' (default), 'spline', 'makima'
%   'SmoothMode'    : 'arclength' (default) or 'time'
%   'UpFactor'      : upsample factor relative to nSeg (default 8)
%   'UseCSAPS'      : true/false, apply smoothing spline after interpolation (default false)
%   'CSAPS_p'       : smoothing parameter [0..1], higher = less smoothing (default 0.85)
%
% Example:
%   plotTrajectoryGradient3D(ax, P, c, 3, 120, 'SmoothMethod','pchip','UpFactor',10);

% -----------------------
% Parse options
% -----------------------
ip = inputParser;
ip.addParameter('SmoothMethod', 'pchip', @(s) ischar(s) || isstring(s));
ip.addParameter('SmoothMode', 'arclength', @(s) any(strcmpi(string(s),["arclength","time"])));
ip.addParameter('UpFactor', 8, @(x) isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('UseCSAPS', false, @islogical);
ip.addParameter('CSAPS_p', 0.85, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
ip.parse(varargin{:});
U = ip.Results;

P = double(P);

% -----------------------
% Basic cleaning
% -----------------------
if isempty(P) || size(P,2) ~= 3
    return;
end

good = all(isfinite(P),2);
P = P(good,:);
T = size(P,1);
if T < 2
    return;
end

% If repeated points exist, arc-length can be zero; we'll handle that.
dP = diff(P,1,1);
step = sqrt(sum(dP.^2,2));

% -----------------------
% Parameterize trajectory
% -----------------------
switch lower(string(U.SmoothMode))
    case "arclength"
        s = [0; cumsum(step)];
        if s(end) < eps
            % all points identical -> nothing to draw
            return;
        end
        tParam = s;
    otherwise
        tParam = (1:T)'; % time index
end

% Ensure strictly increasing parameter for interp
[tParamU, ia] = unique(tParam, 'stable');
PU = P(ia,:);

if size(PU,1) < 2
    return;
end

% -----------------------
% Smooth/upsample
% -----------------------
nUp = max( max(nSeg, 20) * U.UpFactor, 200 );   % enough points for a smooth curve
tFine = linspace(tParamU(1), tParamU(end), nUp)';

method = char(lower(string(U.SmoothMethod)));
% pchip is often the best "looks right" default (no overshoot)
xFine = interp1(tParamU, PU(:,1), tFine, method);
yFine = interp1(tParamU, PU(:,2), tFine, method);
zFine = interp1(tParamU, PU(:,3), tFine, method);

% Optional smoothing spline (great if jitter remains)
if U.UseCSAPS && exist('csaps','file') == 2
    xFine = fnval(csaps(tFine, xFine, U.CSAPS_p), tFine);
    yFine = fnval(csaps(tFine, yFine, U.CSAPS_p), tFine);
    zFine = fnval(csaps(tFine, zFine, U.CSAPS_p), tFine);
end

Pfine = [xFine(:), yFine(:), zFine(:)];

% -----------------------
% Gradient drawing
% -----------------------
% Choose nSeg+1 vertices along the smoothed curve
idx = round(linspace(1, size(Pfine,1), nSeg+1));
idx = unique(max(1, min(size(Pfine,1), idx)));

for i = 1:(numel(idx)-1)
    a = i/(numel(idx)-1);           % 0->1 progression
    col = baseColor*(0.7 + 0.3*a);  % brighten with time (keeps your original style)
    plot3(ax, Pfine(idx(i:i+1),1), Pfine(idx(i:i+1),2), Pfine(idx(i:i+1),3), ...
        '-', 'Color', col, 'LineWidth', lw, 'HandleVisibility','off');
end

end
