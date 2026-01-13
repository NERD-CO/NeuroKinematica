function C = runLMM_PlannedContrasts(lme, varargin)

% runIFR_PlannedContrasts_SU
%
% Planned contrasts for SU LME: DV ~ MoveType*Depth + (1|Subject)
% Reference: MoveType=REST, Depth=b (ventral)
% Intercept: REST @ Depth_b (ventral)
%
% Returns a table with estimates, SE, t, df, p,
% 95% CI, plus optional multiple-comparisons adjustment.

%% ---- Options ----
p = inputParser;
p.addParameter('DVLabel', "", @(x) isstring(x) || ischar(x));
p.addParameter('Alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('DoActiveMeanVsRest', true, @islogical);   % mean(active) vs REST within depth
p.addParameter('DoRestVsEachActive', true, @islogical);   % REST vs each active within depth
p.addParameter('DoActivePairs', true, @islogical);        % active-active pairwise within depth
p.addParameter('DoDepthDiffOfDiff', true, @islogical);    % (move-rest) depth modulation vs refernce depth
p.addParameter('DoDepthPairs', true, @islogical);         % depth modulation of ActivePairs (A-B) across depths vs b
p.addParameter('DoSameMoveAcrossDepth', true, @islogical); % mv@t-vs-b and mv@c-vs-b (includes REST optionally)

p.addParameter('PAdjust', "bonferroni", @(x) isstring(x) || ischar(x)); % "none","holm","bonferroni","fdr"
p.addParameter('AdjustScope', "all", @(x) isstring(x) || ischar(x)); % "all" or "family"
p.parse(varargin{:});
U = p.Results;

DVLabel = string(U.DVLabel);

%% ---- Pull coefficient names + fixed effects + covariance ----
coefNames = string(lme.CoefficientNames(:));
beta      = fixedEffects(lme);              % column vector
CovB      = lme.CoefficientCovariance;      % covariance of fixed effects

% Basic checks
mustHave = ["(Intercept)", ... % REST @ Depth_b (ventral)
    "MoveType_HAND OC","MoveType_HAND PS","MoveType_ARM EF", ...
    "Depth_t","Depth_c"];

missing  = setdiff(mustHave, coefNames);
if ~isempty(missing)
    warning('Some expected coefficients missing: %s', strjoin(missing, ", "));
end

%% ---- Define canonical levels ----

depthList  = ["t","c","b"];                  % dorsal, central, ventral
activeMT   = ["HAND OC","HAND PS","ARM EF"];
allMT      = [activeMT, "REST"];

% For nice labels
depthPretty = containers.Map({'t','c','b'}, {'dorsal STN','central STN','ventral STN'});

%% ---- Storage ----

rows = [];

%% ============================================================
%  A) Within-depth contrasts
% ============================================================

% Helper to add a contrast row
    function addContrastRow(family, depthCode, contrastName, Lrow)
        Lrow = Lrow(:)'; % ensure row
        est  = Lrow * beta;
        se   = sqrt(max(0, Lrow * CovB * Lrow')); % numeric safety
        if se==0
            tStat = NaN;
        else
            tStat = est / se;
        end

        % Use coefTest for p-value + df
        % For 1-df tests, F = t^2 (sign from est)
        [pVal, F, df1, df2] = coefTest(lme, Lrow, 0); %#ok<ASGLU>
        if df1 == 1 && ~isnan(tStat)
            % keep signed tStat from est/se
        end

        % CI using df2 (residual df from coefTest)
        alpha = U.Alpha;
        if isfinite(df2) && se>0
            tCrit = tinv(1 - alpha/2, df2);
            ciLo  = est - tCrit*se;
            ciHi  = est + tCrit*se;
        else
            ciLo = NaN; ciHi = NaN;
        end

        if isKey(depthPretty, char(depthCode))
            depthLabel = string(depthPretty(char(depthCode)));
        else
            depthLabel = "depth " + depthCode;
        end

        newRow = table( ...
            DVLabel, string(family), string(depthCode), depthLabel, ...
            string(contrastName), est, se, tStat, df2, pVal, ciLo, ciHi, ...
            'VariableNames', {'DV','Family','Depth','DepthLabel','Contrast','Estimate','SE','t','DF','p','CI_Lower','CI_Upper'} );

        rows = [rows; newRow];
    end

% Build L vector from coefficient *names*
    function L = makeL_map(nameWeightMap)
        % nameWeightMap: containers.Map where keys are coefficient names (strings/chars)
        L = zeros(1, numel(coefNames));

        k = nameWeightMap.keys;
        for ii = 1:numel(k)
            nm = string(k{ii});
            w  = nameWeightMap(k{ii});

            idx = find(coefNames == nm, 1);
            if isempty(idx)
                error('Coefficient "%s" not found in lme.CoefficientNames.', nm);
            end
            L(idx) = w;
        end
    end

% Convenience: coefficient name builders
    function nm = mtMain(mt)  % MoveType_<mt>
        nm = "MoveType_" + mt;
    end
    function nm = depthMain(depthCode) % Depth_t / Depth_c (b is reference)
        nm = "Depth_" + depthCode;
    end
    function nm = mtDepth(mt, depthCode) % MoveType_<mt>:Depth_t / :Depth_c
        nm = "MoveType_" + mt + ":Depth_" + depthCode;
    end


% A1) REST vs each Active MoveType within each depth
if U.DoRestVsEachActive
    for d = 1:numel(depthList)
        dz = depthList(d);

        for m = 1:numel(activeMT)
            mv = activeMT(m);

            % Active - REST at depth
            w = containers.Map('KeyType','char','ValueType','double');
            w(char(mtMain(mv))) = 1;

            if dz == "t"
                w(char(mtDepth(mv,"t"))) = 1;
            elseif dz == "c"
                w(char(mtDepth(mv,"c"))) = 1;
            % elseif dz == "b" -> nothing (b is reference)
            end

            L = makeL_map(w);
            addContrastRow("RestVsEachActive", dz, mv + " vs REST", L);     % Active - REST
            addContrastRow("RestVsEachActive", dz, "REST vs " + mv, -L);    % REST - Active (often redundant, both directions)
        end
    end
end


% A2) ActiveMean vs REST within depth
if U.DoActiveMeanVsRest
    for d = 1:numel(depthList)
        dz = depthList(d);

        w = containers.Map('KeyType','char','ValueType','double');

        % mean of the three MoveType main effects
        for m = 1:numel(activeMT)
            mv = activeMT(m);
            w(char(mtMain(mv))) = 1/3;

            if dz == "t"
                w(char(mtDepth(mv,"t"))) = 1/3;
            elseif dz == "c"
                w(char(mtDepth(mv,"c"))) = 1/3;
            % elseif dz == "b" -> nothing (b is reference)
            end
        end

        L = makeL_map(w);
        addContrastRow("ActiveMeanVsRest", dz, "Mean(Active) vs REST", L);
        addContrastRow("ActiveMeanVsRest", dz, "REST vs Mean(Active)", -L);
    end
end


% A3) Active pairs within depth
if U.DoActivePairs
    pairs = [ ...
        "HAND OC","HAND PS"; ...
        "HAND OC","ARM EF"; ...
        "HAND PS","ARM EF"  ...
        ];

    for d = 1:numel(depthList)
        dz = depthList(d);

        for k = 1:size(pairs,1)
            A = pairs(k,1);
            B = pairs(k,2);

            w = containers.Map('KeyType','char','ValueType','double');
            w(char(mtMain(A))) =  1;
            w(char(mtMain(B))) = -1;

            if dz == "t"
                w(char(mtDepth(A,"t"))) =  1;
                w(char(mtDepth(B,"t"))) = -1;
            elseif dz == "c"
                w(char(mtDepth(A,"c"))) =  1;
                w(char(mtDepth(B,"c"))) = -1;
            end

            L = makeL_map(w);
            addContrastRow("ActivePairs", dz, A + " vs " + B, L);
            addContrastRow("ActivePairs", dz, B + " vs " + A, -L);
        end
    end
end


%% =====================================================================
%  B) Across-depth contrasts (depth-dependent modulation tests) 
%     "difference-of-differences" 
% ======================================================================

% B1) (Move - REST), Interaction terms vs. Reference depth (ventral)
if U.DoDepthDiffOfDiff
    for mv = activeMT
        w = containers.Map('KeyType','char','ValueType','double');
        w(char(mtDepth(mv,"t"))) = 1;
        addContrastRow("DepthModulation", "t", "(" + mv + "-REST) t vs b", makeL_map(w));

        w = containers.Map('KeyType','char','ValueType','double');
        w(char(mtDepth(mv,"c"))) = 1;
        addContrastRow("DepthModulation", "c", "(" + mv + "-REST) c vs b", makeL_map(w));
    end

    % REST baseline shift across depths (Reference: Depth_b)
    % (REST@t - REST@b = Depth_t; REST@c - REST@b = Depth_c)
    w = containers.Map('KeyType','char','ValueType','double');
    w(char(depthMain("t"))) = 1;
    addContrastRow("DepthModulation", "t", "REST t vs b", makeL_map(w));

    w = containers.Map('KeyType','char','ValueType','double');
    w(char(depthMain("c"))) = 1;
    addContrastRow("DepthModulation", "c", "REST c vs b", makeL_map(w));
end



% B2) Depth modulation of ActivePairs: (A-B)@d - (A-B)@b, depth ref=b
if U.DoDepthPairs
    pairs = [ ...
        "HAND OC","HAND PS"; ...
        "HAND OC","ARM EF"; ...
        "HAND PS","ARM EF"  ...
    ];

    for k = 1:size(pairs,1)
        A = pairs(k,1);
        B = pairs(k,2);

        % t vs b
        w = containers.Map('KeyType','char','ValueType','double');
        w(char(mtDepth(A,"t"))) =  1;
        w(char(mtDepth(B,"t"))) = -1;
        addContrastRow("DepthMod_ActivePair", "t", "(" + A + "-" + B + ") t vs b", makeL_map(w));

        % c vs b
        w = containers.Map('KeyType','char','ValueType','double');
        w(char(mtDepth(A,"c"))) =  1;
        w(char(mtDepth(B,"c"))) = -1;
        addContrastRow("DepthMod_ActivePair", "c", "(" + A + "-" + B + ") c vs b", makeL_map(w));
    end
end

% ----------------------------------------------------------
% B3) Same MoveType across depths: mv@d - mv@b
%     For active mv:
%       mv@t - mv@b = Depth_t + (mv:Depth_t)
%       mv@c - mv@b = Depth_c + (mv:Depth_c)
%     For REST:
%       REST@t - REST@b = Depth_t
%       REST@c - REST@b = Depth_c
% ----------------------------------------------------------
if U.DoSameMoveAcrossDepth
    allMovesForThis = [activeMT, "REST"];  % include REST; remove if you prefer active-only

    for mv = allMovesForThis
        % t vs b
        w = containers.Map('KeyType','char','ValueType','double');
        w(char(depthMain("t"))) = 1;
        if mv ~= "REST"
            w(char(mtDepth(mv,"t"))) = 1;
        end
        addContrastRow("SameMove_AcrossDepth", "t", mv + " t vs b", makeL_map(w));

        % c vs b
        w = containers.Map('KeyType','char','ValueType','double');
        w(char(depthMain("c"))) = 1;
        if mv ~= "REST"
            w(char(mtDepth(mv,"c"))) = 1;
        end
        addContrastRow("SameMove_AcrossDepth", "c", mv + " c vs b", makeL_map(w));
    end
end

%% ---- Assemble table ----

C = rows;

%% ---- Multiple-comparisons adjustment ----

C.p_adj = nan(height(C),1);

method = lower(string(U.PAdjust));
scope  = lower(string(U.AdjustScope));

if method == "none"
    C.p_adj = C.p;
else
    if scope == "family"
        fams = unique(C.Family);
        for f = 1:numel(fams)
            idx = (C.Family == fams(f));
            C.p_adj(idx) = adjustP(C.p(idx), method);
        end
    else
        C.p_adj = adjustP(C.p, method);
    end
end

%% ---- Add a compact "Estimate [CI]" string for manuscript tables ----

C.Est_CI = strings(height(C),1);
for i = 1:height(C)
    if isfinite(C.CI_Lower(i))
        C.Est_CI(i) = sprintf('%.3f [%.3f, %.3f]', C.Estimate(i), C.CI_Lower(i), C.CI_Upper(i));
    else
        C.Est_CI(i) = sprintf('%.3f', C.Estimate(i));
    end
end

%% ---- Sort for readability ----

C = sortrows(C, {'DV','Family','Depth','Contrast'});

end

%% ===================== Local helper: p-adjust =====================

function pAdj = adjustP(pVals, method)
pVals = double(pVals(:));
n = numel(pVals);
pAdj = nan(size(pVals));

switch method
    case "bonferroni"
        pAdj = min(1, pVals * n);

    case "holm"
        [ps, ord] = sort(pVals);
        adj = (n - (1:n)' + 1) .* ps;           % (n-k+1)*p_k
        adj = cummax(adj);                      % enforce monotonicity
        adj = min(1, adj);
        pAdj(ord) = adj;

    case "fdr"  % Benjamini-Hochberg
        [ps, ord] = sort(pVals);
        adj = ps .* n ./ (1:n)';
        adj = cummin(flipud(adj));
        adj = flipud(adj);
        adj = min(1, adj);
        pAdj(ord) = adj;

    otherwise
        error('Unknown PAdjust method: %s', method);
end
end
