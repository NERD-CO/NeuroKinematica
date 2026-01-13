function OUT = run_LMMs(T, dvList, varargin)

% runIFR_LMMs
% Fits LMM per DV: DV ~ 1 + MoveType*Depth + (1|Subject)
% Uses reference coding by default.
%
% Optional varargin name-value:
%   'FitMethod'        : 'ML' (default) or 'REML'
%   'DropUnusedCats'   : true (default) -> removecats after DV NaN-drop
%   'MinRows'          : 5 (default)
%   'MinSubjects'      : 2 (default)
%   'VerboseCoefs'     : false (default) -> print coefficient names

% varargin name-value:
p = inputParser;
p.addParameter('FitMethod','ML', @(x)ischar(x)||isstring(x));
p.parse(varargin{:});
FitMethod = char(p.Results.FitMethod);

% Ensure categoricals
% (Order should be set upstream in your pipeline to enforce references.)
if ~iscategorical(T.Subject),  T.Subject  = categorical(T.Subject);  end
if ~iscategorical(T.MoveType), T.MoveType = categorical(T.MoveType); end
if ~iscategorical(T.Depth),    T.Depth    = categorical(T.Depth);    end

% Model formulas
% f = 'DV ~ MoveType * Depth + (1|Subject)';

% update: accept a list of DV column names and return struct of LME outputs by DV.
OUT = struct();

for i = 1:numel(dvList)
    dv = string(dvList(i));

    if ~ismember(dv, string(T.Properties.VariableNames))
        warning('DV "%s" not found in table. Skipping.', dv);
        continue;
    end

    % Drop missing for *this DV only*
    Ti = T(~isnan(T.(dv)), :);

    if height(Ti) < 5
        warning('DV "%s": too few rows after NaN drop (N=%d). Skipping.', dv, height(Ti));
        continue;
    end

    % Model formula (explicit intercept (1) is fine)
    formula = sprintf('%s ~ 1 + MoveType*Depth + (1|Subject)', dv);

    % Fit linear mixed-effects model
    lme = fitlme(Ti, formula, 'FitMethod', FitMethod, 'DummyVarCoding', 'reference');

    % Store fitted model and ANOVA results in the output structure
    OUT.(dv).lme   = lme;
    OUT.(dv).anova = anova(lme);
    OUT.(dv).N     = height(Ti);

    fprintf('[LME] DV=%s | N=%d | AIC=%.2f\n', dv, height(Ti), lme.ModelCriterion.AIC);
end

end