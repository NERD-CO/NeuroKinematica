function T = build_ZETA_LMM_Table_MUA(MasterZETA_MUA)

% ---- Required vars (match SU pipeline) ----
need = {'SubjectNum','MoveType','Depth','ZetaZ_MUA','ZetaP_MUA'};
missing = setdiff(need, MasterZETA_MUA.Properties.VariableNames);
if ~isempty(missing)
    error('Missing required columns in MasterZETA_MUA: %s', strjoin(missing, ', '));
end

T = table();

% Random effect: Subject as categorical
T.Subject = categorical(string(MasterZETA_MUA.SubjectNum));

% Fixed effects as categorical with explicit reference levels
T.MoveType = categorical(string(MasterZETA_MUA.MoveType), {'REST','HAND OC','HAND PS','ARM EF'}); % REST: reference
T.Depth = categorical(string(MasterZETA_MUA.Depth), {'b','c','t'}); % b (ventral): reference

% Dependent variables (DVs, MUA)
T.ZetaZ_MUA    = double(MasterZETA_MUA.ZetaZ_MUA);
T.ZetaD_MUA    = double(MasterZETA_MUA.ZetaD_MUA);

% Useful for filtering or later reporting (not a DV)
T.ZetaP_MUA = double(MasterZETA_MUA.ZetaP_MUA);

% minimal cleanup:
T = T(~isundefined(T.MoveType) & ~isundefined(T.Depth), :);

end