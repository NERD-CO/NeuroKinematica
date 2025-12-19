function T = build_IFR_LMM_Table_MUA(MasterZETA_MUA)

% Required vars
need = {'SubjectNum','MoveType','Depth','IFR_norm','IFR_baselineNorm'};
missing = setdiff(need, MasterZETA_MUA.Properties.VariableNames);
if ~isempty(missing)
    error('Missing required columns: %s', strjoin(missing, ', '));
end

T = table();

% Subject as categorical
T.Subject = categorical(string(MasterZETA_MUA.SubjectNum));

% Fixed effects as categorical, with explicit reference levels
T.MoveType = categorical(string(MasterZETA_MUA.MoveType), ...
    {'REST','HAND OC','HAND PS','ARM EF'});  % REST first = reference
T.Depth = categorical(string(MasterZETA_MUA.Depth), {'t','c','b'}); % t as reference

% DVs
T.IFR_norm = double(MasterZETA_MUA.IFR_norm);
T.IFR_baselineNorm = double(MasterZETA_MUA.IFR_baselineNorm);

% Drop missing DV rows (do separately per model later if you prefer)
T = T(~isnan(T.IFR_norm) & ~isnan(T.IFR_baselineNorm), :);

end
