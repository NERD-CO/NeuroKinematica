function C = run_PlannedContrasts_NoREST(lme, varargin)

% Planned contrasts for SU LME without REST:
%   DV ~ MoveType*Depth + (1|Subject)
% Reference: MoveType = HAND OC, Depth = b (ventral)
% Intercept: HAND OC @ b
%
% Contrasts:
%   A) Within each depth: Active pairs (HAND OC vs HAND PS, HAND OC vs ARM EF, HAND PS vs ARM EF)
%   B) Across depths (ref=b):
%       - Same MoveType across depths (mv@t vs mv@b; mv@c vs mv@b)
%       - Depth modulation of ActivePairs: (A-B)@t - (A-B)@b and (A-B)@c - (A-B)@b

%% ---- Options ----
p = inputParser;
p.addParameter('DVLabel', "", @(x) isstring(x) || ischar(x));
p.addParameter('Alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('DoActivePairs', true, @islogical);
p.addParameter('DoSameMoveAcrossDepth', true, @islogical);
p.addParameter('DoDepthPairs', true, @islogical);
p.addParameter('PAdjust', "bonferroni", @(x) isstring(x) || ischar(x));
p.addParameter('AdjustScope', "all", @(x) isstring(x) || ischar(x));
p.parse(varargin{:});
U = p.Results;

DVLabel = string(U.DVLabel);

coefNames = string(lme.CoefficientNames(:));
beta      = fixedEffects(lme);
CovB      = lme.CoefficientCovariance;

activeMT  = ["HAND OC","HAND PS","ARM EF"];
depthList = ["b","t","c"]; % b is ref, but we still label outputs by b/t/c

depthPretty = containers.Map({'t','c','b'}, {'dorsal STN','central STN','ventral STN'});

rows = [];

    function addContrastRow(family, depthCode, contrastName, Lrow)
        Lrow = Lrow(:)'; 
        est  = Lrow * beta;
        se   = sqrt(max(0, Lrow * CovB * Lrow'));
        tStat = est / se;

        [pVal,~,df1,df2] = coefTest(lme, Lrow, 0); %#ok<ASGLU>

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

    function L = makeL_map(nameWeightMap)
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

    function nm = mtMain(mt), nm = "MoveType_" + mt; end
    function nm = depthMain(depthCode), nm = "Depth_" + depthCode; end % Depth_t/Depth_c (b ref)
    function nm = mtDepth(mt, depthCode), nm = "MoveType_" + mt + ":Depth_" + depthCode; end

%% A) Within-depth active pairs
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

            % (A - B) at depth dz
            % At b: depends only on MoveType main effects (since Depth ref)
            % At t/c: add Depth interactions for non-reference MoveTypes.
            w = containers.Map('KeyType','char','ValueType','double');

            % Main effects: HAND OC is reference (no coefficient)
            if A ~= "HAND OC", w(char(mtMain(A))) =  1; end
            if B ~= "HAND OC", w(char(mtMain(B))) = -1; end

            if dz == "t"
                if A ~= "HAND OC", w(char(mtDepth(A,"t"))) =  1; end
                if B ~= "HAND OC", w(char(mtDepth(B,"t"))) = -1; end
            elseif dz == "c"
                if A ~= "HAND OC", w(char(mtDepth(A,"c"))) =  1; end
                if B ~= "HAND OC", w(char(mtDepth(B,"c"))) = -1; end
            end

            L = makeL_map(w);
            addContrastRow("ActivePairs_NoREST", dz, A + " vs " + B, L);
            addContrastRow("ActivePairs_NoREST", dz, B + " vs " + A, -L);
        end
    end
end

%% B1) Same MoveType across depths: mv@t - mv@b; mv@c - mv@b
if U.DoSameMoveAcrossDepth
    for mv = activeMT
        % mv@t - mv@b
        w = containers.Map('KeyType','char','ValueType','double');
        w(char(depthMain("t"))) = 1;           % Depth_t always enters
        if mv ~= "HAND OC"
            w(char(mtDepth(mv,"t"))) = 1;      % interaction exists for non-reference MoveType
        end
        addContrastRow("SameMove_AcrossDepth_NoREST", "t", mv + " t vs b", makeL_map(w));

        % mv@c - mv@b
        w = containers.Map('KeyType','char','ValueType','double');
        w(char(depthMain("c"))) = 1;
        if mv ~= "HAND OC"
            w(char(mtDepth(mv,"c"))) = 1;
        end
        addContrastRow("SameMove_AcrossDepth_NoREST", "c", mv + " c vs b", makeL_map(w));
    end
end

%% B2) Depth modulation of ActivePairs: (A-B)@d - (A-B)@b
if U.DoDepthPairs
    pairs = [ ...
        "HAND OC","HAND PS"; ...
        "HAND OC","ARM EF"; ...
        "HAND PS","ARM EF"  ...
    ];

    for k = 1:size(pairs,1)
        A = pairs(k,1);
        B = pairs(k,2);

        % (A-B) t vs b
        w = containers.Map('KeyType','char','ValueType','double');
        if A ~= "HAND OC", w(char(mtDepth(A,"t"))) =  1; end
        if B ~= "HAND OC", w(char(mtDepth(B,"t"))) = -1; end
        addContrastRow("DepthMod_ActivePair_NoREST", "t", "(" + A + "-" + B + ") t vs b", makeL_map(w));

        % (A-B) c vs b
        w = containers.Map('KeyType','char','ValueType','double');
        if A ~= "HAND OC", w(char(mtDepth(A,"c"))) =  1; end
        if B ~= "HAND OC", w(char(mtDepth(B,"c"))) = -1; end
        addContrastRow("DepthMod_ActivePair_NoREST", "c", "(" + A + "-" + B + ") c vs b", makeL_map(w));
    end
end

%% Assemble + p-adjust + formatting
C = rows;

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

C.Est_CI = strings(height(C),1);
for i = 1:height(C)
    if isfinite(C.CI_Lower(i))
        C.Est_CI(i) = sprintf('%.3f [%.3f, %.3f]', C.Estimate(i), C.CI_Lower(i), C.CI_Upper(i));
    else
        C.Est_CI(i) = sprintf('%.3f', C.Estimate(i));
    end
end

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
