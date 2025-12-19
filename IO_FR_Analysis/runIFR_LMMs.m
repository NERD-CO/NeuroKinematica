function OUT = runIFR_LMMs(T, dvList)

% dvList example: ["IFR_norm","IFR_baselineNorm","IFR_mean_Hz","logIFR_mean_Hz"]

if nargin < 2 || isempty(dvList)
    dvList = ["IFR_norm","IFR_baselineNorm"];
end

% Ensure categoricals
T.Subject  = categorical(T.Subject);
T.MoveType = categorical(T.MoveType);
T.Depth    = categorical(T.Depth);

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

    formula = sprintf('%s ~ 1 + MoveType*Depth + (1|Subject)', dv);
    lme = fitlme(T, formula, 'FitMethod','ML');

    OUT.(dv).lme   = lme;
    OUT.(dv).anova = anova(lme);

    fprintf('[LME] DV=%s | N=%d | AIC=%.2f\n', dv, height(T), lme.ModelCriterion.AIC);
end

end