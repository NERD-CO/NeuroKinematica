function out = plot_ZETA_TrajectoryAnalysis_MUA(MasterZETA_MUA, varargin)
% plot_ZETA_TrajectoryAnalysis_MUA
%
% MUA analog of plot_ZETA_TrajectoryAnalysis_SU:
%   - Build mean MUA ZETA trajectories per MoveType per Depth
%   - Project into a common 3D space per depth (trajectory PCA space)
%   - Bootstrap trajectories by resampling units within each MoveType
%   - Compute pairwise Procrustes distances per boot
%   - Save figures + return ProcDist tables + ANOVA tables
%
% REQUIRED columns in MasterZETA_MUA:
%   MoveType, Depth, vecD_MUA, new_vecTime_MUA

% -----------------------
% Parse inputs
% -----------------------
p = inputParser;
p.addParameter('DepthOrder', {'t','c','b'}, @(x)iscell(x)&&~isempty(x));
p.addParameter('ActiveMoveOrder', {'HAND OC','HAND PS','ARM EF'}, @(x)iscell(x)&&~isempty(x));
p.addParameter('SaveDir', '', @(x)ischar(x)||isstring(x));
p.addParameter('FigPrefix', 'MUA_ZETA_TrajFeat', @(x)ischar(x)||isstring(x));
p.addParameter('InterpMethod', 'linear', @(x)ischar(x)||isstring(x));

p.addParameter('SmoothWin', 15, @(x)isnumeric(x)&&isscalar(x)&&x>=1); % in samples (match SU default)
p.addParameter('UseIntegralFeature', true, @islogical);
p.addParameter('DoZscoreFeatures', true, @islogical);
p.addParameter('BootMinUnits', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

p.addParameter('Nboot', 500, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('DoProcrustes', true, @islogical);
p.addParameter('ProcrustesScaling', true, @islogical);
p.addParameter('ProcrustesReflection', false, @islogical);

% Plot
p.addParameter('View', [-40 20], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('LineWidth', 3, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('GradientNSeg', 120, @(x)isnumeric(x)&&isscalar(x)&&x>=10);
p.addParameter('SmoothNSeg', 250, @(x)isnumeric(x)&&isscalar(x)&&x>=20);

p.addParameter('DoSave', true, @islogical);
p.parse(varargin{:});
U = p.Results;

if isempty(U.SaveDir)
    U.SaveDir = fullfile(pwd, 'Trajectory PCA Stats');
end
if ~exist(U.SaveDir,'dir'), mkdir(U.SaveDir); end

% -----------------------
% Sanity checks
% -----------------------
reqVars = {'MoveType','Depth','vecD_MUA','new_vecTime_MUA'};
missing = setdiff(reqVars, MasterZETA_MUA.Properties.VariableNames);
if ~isempty(missing)
    error('MasterZETA_MUA missing required variables: %s', strjoin(missing, ', '));
end

% Normalize types + filter active moves + non-empty
MZ_MUA = MasterZETA_MUA;
MZ_MUA.MoveType = string(MZ_MUA.MoveType);
MZ_MUA.Depth    = string(MZ_MUA.Depth);

activeMoves = string(U.ActiveMoveOrder);
MZ_MUA = MZ_MUA(ismember(MZ_MUA.MoveType, activeMoves), :);

hasZ = ~cellfun(@isempty, MZ_MUA.vecD_MUA) & ~cellfun(@isempty, MZ_MUA.new_vecTime_MUA);
MZ_MUA = MZ_MUA(hasZ, :);

if isempty(MZ_MUA), error('No usable MUA rows after filtering.'); end

% -----------------------
% Colors (match JNE move colors)
% -----------------------
JNE_move = ([ ...
    128,128,128;    % REST  (grey)
    38,116,183;     % HAND OC  (blue)
    53,183,121;     % HAND PS  (green)
    243,120,98] ... % ARM EF   (orange)
    ./ 255);

moveColorMap = containers.Map( ...
    {'HAND OC','HAND PS','ARM EF'}, ...
    {JNE_move(2,:), JNE_move(3,:), JNE_move(4,:)} );

depthNames = containers.Map({'t','c','b'}, {'dorsal STN','central STN','ventral STN'});

% -----------------------
% Output containers
% -----------------------
out = struct();
out.PerDepth = struct([]);
out.Stats = struct();
out.Stats.ANOVA_ProcDist = table();
out.Stats.ANOVA_ProcDist_mc = table();

% Collect all ProcDist values across depths for across-depth ANOVA (SU-style)
allProcRows = {};   % {DepthLabel, PairLabel, ProcDist}

% Per-depth loop
depthOrder = string(U.DepthOrder);
for d = 1:numel(depthOrder)
    dz = char(depthOrder(d));
    MZ_depth = MZ_MUA(MZ_MUA.Depth == string(dz), :);
    if isempty(MZ_depth), continue; end

    % ---------- Build time-aligned matrices per MoveType ----------
    % (need COMMON tRef per depth)
    mvOrder = string(U.ActiveMoveOrder);

    Md_all = MZ_MUA(MZ_MUA.Depth == string(dz) & ismember(MZ_MUA.MoveType, mvOrder), :);
    if isempty(Md_all)
        fprintf('[MUA Traj] Depth %s: no data.\n', dz);
        continue;
    end

    % pick reference time vector from the FULL depth set (not move-specific)
    idxRef = find(~cellfun(@isempty, Md_all.new_vecTime_MUA) & ~cellfun(@isempty, Md_all.vecD_MUA), 1, 'first');
    if isempty(idxRef)
        fprintf('[MUA Traj] Depth %s: no valid time reference.\n', dz);
        continue;
    end
    tRef_depth = Md_all.new_vecTime_MUA{idxRef};
    tRef_depth = double(tRef_depth(:));
    tRef_depth = tRef_depth(isfinite(tRef_depth));
    [tRef_depth,~] = unique(tRef_depth,'stable');

    % Ensure t_use is increasing
    [t_use, sortIdx] = sort(tRef_depth(:));
    t_use = t_use(isfinite(t_use));
    [t_use, ia] = unique(t_use, 'stable'); % define once per depth

    % Build temporal feature vector per MoveType (time x units)
    X_byMove = cell(numel(mvOrder),1);
    Y = nan(numel(t_use), numel(mvOrder));

    for m = 1:numel(mvOrder)
        mv = mvOrder(m);
        Md_mv = Md_all(Md_all.MoveType == mv, :);

        if isempty(Md_mv)
            fprintf('[MUA Traj] Depth %s: missing %s.\n', dz, mv);
            X_byMove{m} = [];
            continue;
        end

        % Align ALL units for this MoveType onto the SAME t_use
        Xmv = build_MUA_ZETA_timeAlignedMatrix_toRef(Md_mv, t_use, 'DepthCode', dz);

        % Xmv = time-aligned [time × units] matrix
        X_byMove{m} = Xmv;                  % store for bootstrapping
        Y(:,m) = mean(Xmv, 2, 'omitnan');   % mean trajectory
    end

    % Require all moves present
    if any(cellfun(@isempty, X_byMove))
        fprintf('[MUA Traj] Depth %s: missing one or more MoveTypes. Skipping.\n', dz);
        continue;
    end

    fprintf('[MUA traj] depth %s | tRef=%d | moves=%d\n', dz, numel(t_use), numel(mvOrder));


    % ---------- Compute non-bootstrapped mean trajectory per MoveType (SU-style) ----------
    trajMean = struct();
    traj = struct([]);

    hFig = figure('Color','w','Units','pixels','Position',[150 80 700 550]);
    ax = axes('Parent', hFig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    view(ax, U.View);
    xlabel(ax,'F1'); ylabel(ax,'F2'); zlabel(ax,'F3');
    title(ax, sprintf('MUA Trajectory Features | %s', depthNames(dz)), 'FontWeight','bold');

    for mIdx = 1:numel(mvOrder)
        mv = char(mvOrder(mIdx));
        Xmv = X_byMove{mIdx};               % [time x units]
        if isempty(Xmv) || size(Xmv,2) < U.BootMinUnits
            continue;
        end

        mu_timecourse = mean(Xmv, 2, 'omitnan');       % mean across units
        mu_timecourse = movmean(mu_timecourse, U.SmoothWin, 'omitnan');

        % Debug
        fprintf('%s | %s: mu range = %.3g (min=%.3g max=%.3g)\n', ...
            depthNames(dz), mv, range(mu_timecourse), min(mu_timecourse), max(mu_timecourse));

        % --- SU-matching trajectory feature mapping ---
        P = buildTrajectoryFeatures(t_use, mu_timecourse, U.UseIntegralFeature, U.DoZscoreFeatures); % [T x 3]

        trajMean.(matlab.lang.makeValidName(mv)) = P;

        % Plot
        c = moveColorMap(mv);
        plotTrajectoryGradient3D(ax, P, c, U.LineWidth, U.GradientNSeg);

        scatter3(ax, P(1,1), P(1,2), P(1,3), 80,  c, 'filled', 'o', 'HandleVisibility','off');
        scatter3(ax, P(end,1), P(end,2), P(end,3), 120, c, 'filled', '^', 'HandleVisibility','off');

        traj(mIdx).MoveType = string(mv);
        traj(mIdx).P = P;
    end

    % Legend
    hLeg = gobjects(numel(mvOrder),1);
    for mIdx = 1:numel(mvOrder)
        mv = char(mvOrder(mIdx));
        hLeg(mIdx) = scatter3(ax, NaN,NaN,NaN, 90, moveColorMap(mv), 'filled');
    end
    legend(ax, hLeg, cellstr(mvOrder), 'Location','northeast');


    % ---------- Bootstrap Procrustes distances ----------
    ProcDist = table();

    if U.DoProcrustes && U.Nboot > 0
        pairList = nchoosek(1:numel(mvOrder), 2);
        rows = cell(U.Nboot*size(pairList,1), 5);
        rr = 0;

        for b = 1:U.Nboot
            Pboot = cell(numel(mvOrder),1);
            ok = true;

            for mIdx = 1:numel(mvOrder)
                Xmv = X_byMove{mIdx};  % time x units
                if isempty(Xmv) || size(Xmv,2) < U.BootMinUnits
                    ok = false; break;
                end

                nU = size(Xmv,2);
                sampCols = randi(nU, [nU 1]); % resample units WITH replacement
                mu_b = mean(Xmv(:, sampCols), 2, 'omitnan');
                mu_b = movmean(mu_b, U.SmoothWin, 'omitnan');

                Pboot{mIdx} = buildTrajectoryFeatures(t_use, mu_b, U.UseIntegralFeature, U.DoZscoreFeatures);
            end
            if ~ok, continue; end

            for pIdx = 1:size(pairList,1)
                iA = pairList(pIdx,1);
                iB = pairList(pIdx,2);

                dProc = procrustesDistanceSafe(Pboot{iA}, Pboot{iB}, ...
                    U.ProcrustesScaling, U.ProcrustesReflection);

                rr = rr + 1;
                rows(rr,:) = {string(depthNames(dz)), string(mvOrder(iA)), string(mvOrder(iB)), b, dProc};
                allProcRows(end+1,:) = {string(depthNames(dz)), sprintf('%s vs %s', mvOrder(iA), mvOrder(iB)), dProc}; 
            end
        end

        rows = rows(1:rr,:); % trim
        ProcDist = cell2table(rows, 'VariableNames', {'Depth','MoveA','MoveB','BootIter','ProcDist'});
    end

    % ---------- Save ----------
    if U.DoSave
        base = sprintf('%s_%s', char(string(U.FigPrefix)), dz);
        exportgraphics(hFig, fullfile(U.SaveDir, [base '.png']), 'Resolution', 300);
        savefig(hFig, fullfile(U.SaveDir, [base '.fig']));
    end

    % ---------- Store outputs ----------
    out.PerDepth(end+1).Depth = dz;
    out.PerDepth(end).DepthLabel = string(depthNames(dz));
    out.PerDepth(end).t = t_use;
    out.PerDepth(end).Traj = traj;
    out.PerDepth(end).TrajMean  = trajMean;
    out.PerDepth(end).ProcDist = ProcDist;

end % end per-depth loop


% -----------------------
% Optional: ANOVA across depths for each movement-pair (SU-style)
% -----------------------
out.ANOVA = struct('tbl', table(), 'mc', table());

if ~isempty(allProcRows)
    Dall = cell2table(allProcRows, 'VariableNames', {'Depth','Pair','ProcDist'});

    pairs = unique(Dall.Pair);
    anovaRows = table();
    mcAll = table();

    for i = 1:numel(pairs)
        thisPair = pairs(i);

        Tpair = Dall(strcmp(string(Dall.Pair), string(thisPair)), :);

        if numel(unique(Tpair.Depth)) < 2
            continue;
        end

        [p_aov, tbl_aov, stats_aov] = anova1(Tpair.ProcDist, categorical(Tpair.Depth), 'off');

        tblAOV = cell2table(tbl_aov(2:end,:), ...
            'VariableNames', matlab.lang.makeValidName(string(tbl_aov(1,:))));
        tblAOV.Pair = repmat(string(thisPair), height(tblAOV), 1);
        tblAOV.pValue_overall = repmat(p_aov, height(tblAOV), 1);
        tblAOV = movevars(tblAOV, {'Pair','pValue_overall'}, 'Before', 1);

        anovaRows = [anovaRows; tblAOV]; %#ok<AGROW>

        % Bonferroni posthoc
        mc = multcompare(stats_aov, 'CType','bonferroni', 'Display','off');
        mcTbl = array2table(mc, 'VariableNames', {'Group1','Group2','LowerCI','Diff','UpperCI','pValue'});

        grpNames = categories(categorical(Tpair.Depth));
        mcTbl.Group1 = grpNames(mcTbl.Group1);
        mcTbl.Group2 = grpNames(mcTbl.Group2);

        mcTbl.Pair = repmat(string(thisPair), height(mcTbl), 1);
        mcTbl = movevars(mcTbl, 'Pair', 'Before', 1);

        mcAll = [mcAll; mcTbl]; %#ok<AGROW>
    end

    out.ANOVA.tbl = anovaRows;
    out.ANOVA.mc  = mcAll;

    % also populate your existing out.Stats fields for backwards compatibility
    out.Stats.ANOVA_ProcDist    = out.ANOVA.tbl;
    out.Stats.ANOVA_ProcDist_mc = out.ANOVA.mc;
end

% -----------------------
% Save OUT + tables (SU-style)
% -----------------------
if U.DoSave
    base = char(string(U.FigPrefix));

    % Save whole output struct
    outMatPath = fullfile(U.SaveDir, [base '_OUT.mat']);
    save(outMatPath, 'out', '-v7.3');

    % Save per-depth ProcDist CSVs
    if ~isempty(out.PerDepth)
        for k = 1:numel(out.PerDepth)
            if isfield(out.PerDepth(k),'ProcDist') && istable(out.PerDepth(k).ProcDist) && ~isempty(out.PerDepth(k).ProcDist)
                writetable(out.PerDepth(k).ProcDist, ...
                    fullfile(U.SaveDir, sprintf('%s_ProcDist_%s.csv', base, out.PerDepth(k).Depth)));
            end
        end
    end

    % Save ANOVA tables
    if isfield(out,'ANOVA') && istable(out.ANOVA.tbl) && ~isempty(out.ANOVA.tbl)
        writetable(out.ANOVA.tbl, fullfile(U.SaveDir, [base '_ANOVA_ProcDist_tbl.csv']));
    end
    if isfield(out,'ANOVA') && istable(out.ANOVA.mc) && ~isempty(out.ANOVA.mc)
        writetable(out.ANOVA.mc, fullfile(U.SaveDir, [base '_ANOVA_ProcDist_mc.csv']));
    end

    fprintf('[DONE] Saved outputs to: %s\n', U.SaveDir);
end
end

% =========================================================
% Helpers
% =========================================================

function [X, mtPerUnit, labelsUnit] = build_MUA_ZETA_timeAlignedMatrix_toRef(M, tRef, varargin)
% Align vecD_MUA to a PROVIDED time reference tRef (column vector).
% Returns X [time x units].

p = inputParser;
p.addParameter('DepthCode','', @(x)ischar(x) || isstring(x));
p.parse(varargin{:});
depthCode = char(p.Results.DepthCode);

reqVars = {'MoveType','vecD_MUA','new_vecTime_MUA'};
missing = setdiff(reqVars, M.Properties.VariableNames);
if ~isempty(missing)
    error('build_MUA_ZETA_timeAlignedMatrix_toRef: missing: %s', strjoin(missing,', '));
end

tRef = double(tRef(:));
tRef = tRef(isfinite(tRef));
[tRef,~] = unique(tRef,'stable');

nUnits    = height(M);
X         = nan(numel(tRef), nUnits);
keepUnit  = false(1,nUnits);
mtPerUnit = strings(nUnits,1);
labelsUnit= strings(nUnits,1);

% label column if present
varNames = M.Properties.VariableNames;
if ismember('PrettyLabel', varNames)
    labelVar = 'PrettyLabel';
elseif ismember('ChannelLabel', varNames)
    labelVar = 'ChannelLabel';
else
    labelVar = '';
end

for u = 1:nUnits
    dVec = M.vecD_MUA{u};
    tVec = M.new_vecTime_MUA{u};

    if isempty(dVec) || isempty(tVec), continue; end
    dVec = double(dVec(:));
    tVec = double(tVec(:));

    good = isfinite(dVec) & isfinite(tVec);
    dVec = dVec(good);
    tVec = tVec(good);

    if numel(dVec) < 6 || numel(tVec) ~= numel(dVec)
        continue;
    end

    % interp safety
    [tU, iu] = unique(tVec, 'stable');
    dU = dVec(iu);
    if numel(tU) < 6, continue; end

    try
        X(:,u) = interp1(tU, dU, tRef, 'linear', 'extrap');
        keepUnit(u)  = true;
        mtPerUnit(u) = string(M.MoveType(u));
        if ~isempty(labelVar), labelsUnit(u) = string(M.(labelVar)(u));
        else, labelsUnit(u) = "MUA_" + u;
        end
    catch ME
        warning('interp1 failed for MUA unit %d at depth %s: %s', u, depthCode, ME.message);
    end
end

X = X(:,keepUnit);
goodT = any(isfinite(X),2);
X = X(goodT,:);
mtPerUnit = mtPerUnit(keepUnit);
labelsUnit = labelsUnit(keepUnit);

if isempty(X)
    error('No valid MUA units after alignment at depth "%s".', depthCode);
end

% Fill remaining NaNs with column means
for j = 1:size(X,2)
    nanMask = ~isfinite(X(:,j));
    if any(nanMask)
        mu = mean(X(~nanMask,j),'omitnan');
        if ~isfinite(mu), mu = 0; end
        X(nanMask,j) = mu;
    end
end
end


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


function Psm = smoothTrajectory3D(P, nSeg)
% Smooth trajectory by resampling in index-space with PCHIP (no overshoot)
P = double(P);
T = size(P,1);
if T < 5
    Psm = P;
    return;
end
xi = 1:T;
xq = linspace(1, T, nSeg);
Psm = [ ...
    interp1(xi, P(:,1), xq, 'pchip')', ...
    interp1(xi, P(:,2), xq, 'pchip')', ...
    interp1(xi, P(:,3), xq, 'pchip')'];
end

function plotTrajectoryGradient3D(ax, P, baseColor, lw, nSeg)
% time-gradient polyline (smoothed P assumed)
T = size(P,1);
idx = round(linspace(1, T, nSeg));
idx = unique(max(1, min(T, idx)));

for i = 1:(numel(idx)-1)
    a = i/(numel(idx)-1);
    col = baseColor*(0.7 + 0.3*a);
    plot3(ax, P(idx(i:i+1),1), P(idx(i:i+1),2), P(idx(i:i+1),3), ...
        '-', 'Color', col, 'LineWidth', lw, 'HandleVisibility','off');
end
end

function d = procrustesDistanceSafe(A, B, doScale, doReflect)
% Safe wrapper around procrustes output fields
% MATLAB returns: [d, Z, transform] = procrustes(X, Y, ...)
A = double(A); B = double(B);

% match length by interpolation (index-space)
n = min(size(A,1), size(B,1));
if size(A,1) ~= n
    A = smoothTrajectory3D(A, n);
end
if size(B,1) ~= n
    B = smoothTrajectory3D(B, n);
end

[d,~,~] = procrustes(A, B, ...
    'Scaling', doScale, ...
    'Reflection', doReflect);
end
