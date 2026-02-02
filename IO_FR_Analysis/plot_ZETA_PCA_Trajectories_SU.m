function out = plot_ZETA_PCA_Trajectories_SU(MasterZETA, varargin)
% plot_ZETA_PCA_Trajectories_SU
%
% Outputs:
%   Figure 1: 3 panels (t/c/b). For each depth: 3D PC trajectory per MoveType
%             with start marker, time-gradient path, end marker.
%   Figure 2: centroid + covariance ellipsoids in shared PC space per depth.
%   Procrustes distance tables (pairwise per depth).
%   Bootstrap distributions of Procrustes distances for ANOVA-style testing.
%
% Notes:
%   - PCA is computed per depth on X_all (time x units), with observations=time.
%   - Each MoveType trajectory is obtained by "mean-holding" non-move units.

% -----------------------
% Parse inputs
% -----------------------
p = inputParser;

p.addParameter('DepthOrder', {'t','c','b'}, @(x) iscell(x) && ~isempty(x));
p.addParameter('ActiveMoveOrder', {'HAND OC','HAND PS','ARM EF'}, @(x) iscell(x) && ~isempty(x));

p.addParameter('SaveDir', '', @(x)ischar(x)||isstring(x));
p.addParameter('FigPrefix', 'SU_ZETA_PCatraj', @(x)ischar(x)||isstring(x));

% Time alignment
p.addParameter('TimeSource', 'PSTH_TimeCenters_s', @(s) any(strcmpi(string(s), ["PSTH_TimeCenters_s","ZETA_vecT"])));
p.addParameter('InterpMethod', 'linear', @(s)ischar(s)||isstring(s));

% PCA / plotting
p.addParameter('NumPC', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=3);
p.addParameter('FlipPC1PostOnset', true, @islogical);
p.addParameter('PostOnsetWindow_s', [0 0.8], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('View', [-40 20], @(x)isnumeric(x)&&numel(x)==2);

% Trajectory styling
p.addParameter('StartMarkerSize', 60, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('EndMarkerSize', 80, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('LineWidth', 3, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('GradientN', 250, @(x)isnumeric(x)&&isscalar(x)&&x>=20);

% Figure 2 ellipsoids
p.addParameter('EllipsoidMode', 'chi2_95', @(s)ischar(s)||isstring(s)); % '1sigma' or 'chi2_95'
p.addParameter('EllipsoidAlpha', 0.10, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);

% Distances / bootstrap
p.addParameter('DoProcrustes', true, @islogical);
p.addParameter('ProcrustesScaling', true, @islogical);
p.addParameter('ProcrustesReflection', false, @islogical); % usually false for neural trajectories
p.addParameter('Nboot', 250, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('BootSeed', 1, @(x)isnumeric(x)&&isscalar(x));

p.parse(varargin{:});
U = p.Results;

if isempty(U.SaveDir)
    U.SaveDir = fullfile(pwd, '3D PCA Trajectories');
end
if ~exist(U.SaveDir,'dir')
    mkdir(U.SaveDir);
end

% -----------------------
% Required columns
% -----------------------
reqVars = {'MoveType','Depth','ZETA_vecD','ZETA_vecT'};
missing = setdiff(reqVars, MasterZETA.Properties.VariableNames);
if ~isempty(missing)
    error('MasterZETA missing required variables: %s', strjoin(missing, ', '));
end
if strcmpi(U.TimeSource,'PSTH_TimeCenters_s') && ~ismember('PSTH_TimeCenters_s', MasterZETA.Properties.VariableNames)
    warning('TimeSource=PSTH_TimeCenters_s requested but column missing. Falling back to ZETA_vecT.');
    U.TimeSource = 'ZETA_vecT';
end

% Normalize types + filter
MT = MasterZETA;
MT.MoveType = string(MT.MoveType);
MT.Depth    = string(MT.Depth);

activeMoves = string(U.ActiveMoveOrder);
MT = MT(ismember(MT.MoveType, activeMoves), :);

hasZ = ~cellfun(@isempty, MT.ZETA_vecD) & ~cellfun(@isempty, MT.ZETA_vecT);
MT = MT(hasZ, :);
if isempty(MT), error('No usable rows after filtering.'); end

% -----------------------
% Color scheme (your new one)
% -----------------------
JNE_move = ([ ...
    38,116,183;     % HAND OC  (blue)
    53,183,121;     % HAND PS  (green/teal)
    243,120,98] ... % ARM EF   (coral)
    ./ 255);

moveColorMap = containers.Map( ...
    {'HAND OC','HAND PS','ARM EF'}, ...
    {JNE_move(1,:), JNE_move(2,:), JNE_move(3,:)} );

depthNames = containers.Map({'t','c','b'}, {'dorsal STN','central STN','ventral STN'});

depthOrder = string(U.DepthOrder);

% Outputs
out = struct();
out.PerDepth = struct([]);
out.Procrustes = table();
out.ProcrustesBoot = table();
out.ANOVA = struct();

% -----------------------
% Figure 1: 3 panels (t/c/b)
% -----------------------
hF1 = figure('Color','w','Units','pixels','Position',[80 60 1100 950]);
tlo1 = tiledlayout(hF1, 3, 1, 'TileSpacing','compact','Padding','compact');

% -----------------------
% Figure 2: 3 panels (t/c/b)
% -----------------------
hF2 = figure('Color','w','Units','pixels','Position',[1240 60 1100 950]);
tlo2 = tiledlayout(hF2, 3, 1, 'TileSpacing','compact','Padding','compact');

rng(U.BootSeed);

for d = 1:numel(depthOrder)
    dz = char(depthOrder(d));
    Md = MT(MT.Depth == string(dz), :);
    if height(Md) < 3
        continue;
    end

    % Build time-aligned matrix: X_all (time x units)
    [X_all, tRef, mvPerUnit] = build_timeAlignedMatrix(Md, U.TimeSource, U.InterpMethod);
    if size(X_all,2) < 3, continue; end

    % --- PCA in unit-space with observations = time ---
    % score_time: (time x PCs), coeff_units: (units x PCs)
    [coeffU, scoreT, latent, ~, explained, muU] = pca(X_all, 'Centered', true);
    if size(scoreT,2) < U.NumPC
        % pad if needed
        scoreT(:,end+1:U.NumPC) = 0;
        coeffU(:,end+1:U.NumPC) = 0;
        explained(end+1:U.NumPC) = 0;
    end

    % flip PC1 based on post-onset window (using the TIME scores)
    if U.FlipPC1PostOnset
        postMask = (tRef >= U.PostOnsetWindow_s(1)) & (tRef <= U.PostOnsetWindow_s(2));
        if any(postMask)
            if mean(scoreT(postMask,1), 'omitnan') < 0
                scoreT(:,1) = -scoreT(:,1);
                coeffU(:,1) = -coeffU(:,1);
            end
        end
    end

    % Shared PCA basis per depth
    PCbasis = coeffU(:,1:3);      % units x 3
    mu = muU(:)';                 % 1 x units

    % Compute trajectories for each MoveType in shared PC space
    mvOrder = string(U.ActiveMoveOrder);
    traj = struct();

    for mIdx = 1:numel(mvOrder)
        mv = char(mvOrder(mIdx));
        cols = (mvPerUnit == string(mv));
        if nnz(cols) < 2
            traj.(matlab.lang.makeValidName(mv)).S = nan(numel(tRef),3);
            continue;
        end

        % Hold non-move units at mean so they contribute 0 after centering
        X_mv = repmat(mu, size(X_all,1), 1);
        X_mv(:, cols) = X_all(:, cols);

        S = (X_mv - mu) * PCbasis;   % time x 3
        traj.(matlab.lang.makeValidName(mv)).S = S;
    end

    % -----------------------
    % Figure 1 panel (trajectory + start/end + gradient)
    % -----------------------
    ax1 = nexttile(tlo1, d); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');

    % Plot each move trajectory
    for mIdx = 1:numel(mvOrder)
        mv = char(mvOrder(mIdx));
        S  = traj.(matlab.lang.makeValidName(mv)).S;
        if any(~isfinite(S(:))), continue; end

        c0 = moveColorMap(mv);

        % optional resample to make gradient smooth and consistent length
        Sg = resampleTrajectory(S, U.GradientN);

        % Gradient line (mixes toward white over time)
        plotGradient3(ax1, Sg, c0, U.LineWidth);

        % Start marker (dot)
        scatter3(ax1, Sg(1,1), Sg(1,2), Sg(1,3), ...
            U.StartMarkerSize, c0, 'filled', 'MarkerEdgeColor','none');

        % % End marker (triangle)
        % scatter3(ax1, Sg(end,1), Sg(end,2), Sg(end,3), ...
        %     U.EndMarkerSize, c0, '^', 'filled', 'MarkerEdgeColor','none');
    end

    xlabel(ax1,'PC1'); ylabel(ax1,'PC2'); zlabel(ax1,'PC3');
    ttl1 = sprintf('SU PCA trajectory | %s', depthNames(dz));
    title(ax1, ttl1, 'Interpreter','none');
    view(ax1, U.View);
    axis(ax1,'vis3d');

    % Legend: dummy handles (solid colors)
    hLeg = gobjects(numel(mvOrder),1);
    for mIdx = 1:numel(mvOrder)
        mv = char(mvOrder(mIdx));
        hLeg(mIdx) = scatter3(ax1, NaN,NaN,NaN, 80, moveColorMap(mv), 'filled', ...
            'MarkerEdgeColor','none');
    end
    lgd = legend(ax1, hLeg, cellstr(mvOrder), 'Location','northeastoutside');
    lgd.AutoUpdate = 'off';

    % -----------------------
    % Figure 2 panel (centroid + covariance ellipsoid over timepoints)
    % -----------------------
    ax2 = nexttile(tlo2, d); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');

    for mIdx = 1:numel(mvOrder)
        mv = char(mvOrder(mIdx));
        S  = traj.(matlab.lang.makeValidName(mv)).S;
        if any(~isfinite(S(:))), continue; end

        c0 = moveColorMap(mv);

        muS = mean(S, 1, 'omitnan');      % 1x3
        C   = cov(S, 'omitrows');         % 3x3

        % centroid
        scatter3(ax2, muS(1), muS(2), muS(3), 120, c0, 'filled', ...
            'MarkerEdgeColor','k', 'LineWidth', 0.5);

        % ellipsoid
        plotCovEllipsoid(ax2, muS, C, c0, U.EllipsoidMode, U.EllipsoidAlpha);
    end

    xlabel(ax2,'PC1'); ylabel(ax2,'PC2'); zlabel(ax2,'PC3');
    ttl2 = sprintf('SU centroid + dispersion | %s', depthNames(dz));
    title(ax2, ttl2, 'Interpreter','none');
    view(ax2, U.View);
    axis(ax2,'vis3d');

    % Save per-depth outputs
    out.PerDepth(end+1).Depth = dz;
    out.PerDepth(end).DepthName = depthNames(dz);
    out.PerDepth(end).tRef = tRef;
    out.PerDepth(end).PCbasis_units = PCbasis;
    out.PerDepth(end).mu_units = mu;
    out.PerDepth(end).ExplainedPct = explained(1:3);
    out.PerDepth(end).Traj = traj;

    % -----------------------
    % Procrustes distances (pairwise) per depth
    % -----------------------
    if U.DoProcrustes
        pairs = nchoosek(1:numel(mvOrder), 2);
        for pIdx = 1:size(pairs,1)
            a = char(mvOrder(pairs(pIdx,1)));
            b = char(mvOrder(pairs(pIdx,2)));

            Sa = traj.(matlab.lang.makeValidName(a)).S;
            Sb = traj.(matlab.lang.makeValidName(b)).S;
            if any(~isfinite(Sa(:))) || any(~isfinite(Sb(:))), continue; end

            Sa = resampleTrajectory(Sa, U.GradientN);
            Sb = resampleTrajectory(Sb, U.GradientN);

            dProc = procrustesDistance(Sa, Sb, U.ProcrustesScaling, U.ProcrustesReflection);

            out.Procrustes = [out.Procrustes; table(string(depthNames(dz)), string(a), string(b), dProc, ...
                'VariableNames', {'Depth','MoveA','MoveB','ProcrustesD'})];
        end
    end

    % -----------------------
    % Bootstrap Procrustes distances (for ANOVA framework)
    % -----------------------
    if U.Nboot > 0
        pairs = nchoosek(1:numel(mvOrder), 2);

        % Cache unit columns per MoveType (in X_all space)
        colIdx = struct();
        for mIdx = 1:numel(mvOrder)
            mv = char(mvOrder(mIdx));
            colIdx.(matlab.lang.makeValidName(mv)) = find(mvPerUnit == string(mv));
        end

        for bIter = 1:U.Nboot
            trajB = struct();

            for mIdx = 1:numel(mvOrder)
                mv = char(mvOrder(mIdx));
                cols = colIdx.(matlab.lang.makeValidName(mv));
                if numel(cols) < 2
                    trajB.(matlab.lang.makeValidName(mv)) = nan(numel(tRef),3);
                    continue;
                end

                % bootstrap sample columns WITH replacement
                colsB = cols(randi(numel(cols), [numel(cols) 1]));

                X_mv = repmat(mu, size(X_all,1), 1);
                X_mv(:, cols) = repmat(mu(1, cols), size(X_mv,1), 1);   % keep original set at mean
                X_mv(:, colsB) = X_all(:, colsB);

                S = (X_mv - mu) * PCbasis;
                trajB.(matlab.lang.makeValidName(mv)) = resampleTrajectory(S, U.GradientN);
            end

            % distances
            for pIdx = 1:size(pairs,1)
                a = char(mvOrder(pairs(pIdx,1)));
                b = char(mvOrder(pairs(pIdx,2)));
                Sa = trajB.(matlab.lang.makeValidName(a));
                Sb = trajB.(matlab.lang.makeValidName(b));
                if any(~isfinite(Sa(:))) || any(~isfinite(Sb(:))), continue; end

                dProc = procrustesDistance(Sa, Sb, U.ProcrustesScaling, U.ProcrustesReflection);

                out.ProcrustesBoot = [out.ProcrustesBoot; table( ...
                    bIter, string(depthNames(dz)), string(a), string(b), dProc, ...
                    'VariableNames', {'BootIter','Depth','MoveA','MoveB','ProcrustesD'})];
            end
        end
    end
end

% -----------------------
% Optional ANOVA framework on bootstrap distances
% -----------------------
if ~isempty(out.ProcrustesBoot)
    T = out.ProcrustesBoot;

    % Create a Pair label
    pair = strings(height(T),1);
    for i = 1:height(T)
        pair(i) = T.MoveA(i) + " vs " + T.MoveB(i);
    end
    T.Pair = categorical(pair);

    % Example: per-depth one-way ANOVA across movement pairs
    depths = categories(categorical(T.Depth));
    anovaRows = table();

    for i = 1:numel(depths)
        dd = depths{i};
        Td = T(T.Depth == dd, :);
        if height(Td) < 10, continue; end

        [pA, tblA, statsA] = anova1(Td.ProcrustesD, Td.Pair, 'off');

        tblA = cell2table(tblA(2:end,:), 'VariableNames', matlab.lang.makeValidName(string(tblA(1,:))));
        tblA.Depth = repmat(string(dd), height(tblA), 1);
        tblA = movevars(tblA, 'Depth', 'Before', 1);
        tblA.pValue_overall = repmat(pA, height(tblA), 1);

        mc = multcompare(statsA, 'CType','bonferroni', 'Display','off');
        mcTbl = array2table(mc, 'VariableNames', {'Group1','Group2','LowerCI','Diff','UpperCI','pValue'});
        grp = categories(Td.Pair);
        mcTbl.Group1 = grp(mcTbl.Group1);
        mcTbl.Group2 = grp(mcTbl.Group2);
        mcTbl.Depth  = repmat(string(dd), height(mcTbl), 1);
        mcTbl = movevars(mcTbl, 'Depth', 'Before', 1);

        anovaRows = [anovaRows; tblA]; %#ok<AGROW>

        out.ANOVA.(matlab.lang.makeValidName(char(dd))).tbl = tblA;
        out.ANOVA.(matlab.lang.makeValidName(char(dd))).mc  = mcTbl;
    end

    out.ANOVA_AllDepths_tbl = anovaRows;
end

% -----------------------
% Save figures + outputs
% -----------------------
baseOut = char(string(U.FigPrefix));

exportgraphics(hF1, fullfile(U.SaveDir, [baseOut '_Figure1_Trajectories.png']), 'Resolution',300);
savefig(hF1, fullfile(U.SaveDir, [baseOut '_Figure1_Trajectories.fig']));

exportgraphics(hF2, fullfile(U.SaveDir, [baseOut '_Figure2_CentroidEllipsoids.png']), 'Resolution',300);
savefig(hF2, fullfile(U.SaveDir, [baseOut '_Figure2_CentroidEllipsoids.fig']));

save(fullfile(U.SaveDir, [baseOut '_OUT.mat']), 'out', '-v7.3');

if istable(out.Procrustes) && ~isempty(out.Procrustes)
    writetable(out.Procrustes, fullfile(U.SaveDir, [baseOut '_Procrustes.csv']));
end
if istable(out.ProcrustesBoot) && ~isempty(out.ProcrustesBoot)
    writetable(out.ProcrustesBoot, fullfile(U.SaveDir, [baseOut '_ProcrustesBoot.csv']));
end
if isfield(out,'ANOVA_AllDepths_tbl') && istable(out.ANOVA_AllDepths_tbl) && ~isempty(out.ANOVA_AllDepths_tbl)
    writetable(out.ANOVA_AllDepths_tbl, fullfile(U.SaveDir, [baseOut '_BootDist_ANOVA_tbl.csv']));
end

fprintf('[DONE] Saved figures + outputs to: %s\n', U.SaveDir);

end % main


%% ==========================
% Helpers (local)
% ==========================

function [X, tRef, mvPerUnit] = build_timeAlignedMatrix(Md, timeSource, interpMethod)
switch lower(string(timeSource))
    case "psth_timecenters_s"
        tRef0 = Md.PSTH_TimeCenters_s{find(~cellfun(@isempty, Md.PSTH_TimeCenters_s), 1, 'first')};
    otherwise
        tRef0 = Md.ZETA_vecT{find(~cellfun(@isempty, Md.ZETA_vecT), 1, 'first')};
end

tRef0 = double(tRef0(:));
tRef0 = tRef0(isfinite(tRef0));
[tRef, ~] = unique(tRef0, 'stable');
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

    if numel(tVec) < 6 || numel(dVec) ~= numel(tVec), continue; end

    [tU, iu] = unique(tVec, 'stable');
    dU = dVec(iu);

    try
        X(:,u) = interp1(tU, dU, tRef, interpMethod, 'extrap');
    catch
    end
end

% drop all-NaN units
keepU = any(isfinite(X),1);
X = X(:,keepU);
mvPerUnit = mvPerUnit(keepU);

% fill remaining NaNs with unit mean
for j = 1:size(X,2)
    nanMask = ~isfinite(X(:,j));
    if any(nanMask)
        mu = mean(X(~nanMask,j), 'omitnan');
        X(nanMask,j) = mu;
    end
end

% drop any fully-bad timepoints
goodT = any(isfinite(X),2);
X = X(goodT,:);
tRef = tRef(goodT);
end


function Sg = resampleTrajectory(S, N)
% Resample trajectory to N points (by arc-length parameterization)
S = double(S);
d = sqrt(sum(diff(S,1,1).^2,2));
s = [0; cumsum(d)];
if s(end) < eps
    Sg = repmat(S(1,:), N, 1);
    return;
end
sq = linspace(0, s(end), N);
Sg = interp1(s, S, sq, 'linear');
end


function plotGradient3(ax, S, baseColor, lw)
% Plot 3D line with a time-gradient (baseColor -> lighter)
N = size(S,1);
w = linspace(0, 1, N);
for i = 1:(N-1)
    % blend toward white
    c = (1-w(i))*baseColor + w(i)*[1 1 1]*0.65;
    plot3(ax, S(i:i+1,1), S(i:i+1,2), S(i:i+1,3), 'Color', c, 'LineWidth', lw);
end
end


function plotCovEllipsoid(ax, muS, C, color, mode, alphaVal)
% Plot covariance ellipsoid for 3D data with mean muS and covariance C
if any(~isfinite(C(:))) || rank(C) < 3
    return;
end

% scale
switch lower(string(mode))
    case "1sigma"
        k = 1;
    otherwise
        % 95% chi-square for df=3
        k = sqrt(chi2inv(0.95, 3));
end

[V,D] = eig(C);
r = k * sqrt(max(diag(D), 0));  % radii

% unit sphere
[u,v,w] = sphere(28);
X = r(1)*u; Y = r(2)*v; Z = r(3)*w;

% rotate and translate
P = V * [X(:)'; Y(:)'; Z(:)'];
Xr = reshape(P(1,:), size(X)) + muS(1);
Yr = reshape(P(2,:), size(Y)) + muS(2);
Zr = reshape(P(3,:), size(Z)) + muS(3);

surf(ax, Xr, Yr, Zr, ...
    'FaceColor', color, 'FaceAlpha', alphaVal, ...
    'EdgeColor', 'none');
end


function d = procrustesDistance(A, B, doScaling, doReflection)
% Compute Procrustes distance (scalar dissimilarity) between two trajectories
%
% A, B: [T x 3] matrices (time-aligned PC trajectories)
% doScaling: true/false
% doReflection: 'best' | true | false

    % --- Procrustes ---
    % [d, Z, transform] = procrustes(X, Y, ...)
    % A and B must be Nx3 (or NxD) with same N
    [d, ~, ~] = procrustes(A, B, ...
        'Scaling',    logical(doScaling), ...
        'Reflection', reflectionFlag(doReflection));

end

function rf = reflectionFlag(doReflection)
    % procrustes expects: true/false OR 'best'
    if ischar(doReflection) || isstring(doReflection)
        rf = char(doReflection); % allow 'best'
    else
        rf = logical(doReflection); % true/false
    end
end