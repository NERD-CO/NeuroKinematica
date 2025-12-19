function T = build_IFR_LMM_Table(MasterZETA)

% Required vars
need = {'SubjectNum','MoveType','Depth', 'IFR_mean_Hz', 'IFR_norm','IFR_baselineNorm'};
missing = setdiff(need, MasterZETA.Properties.VariableNames);
if ~isempty(missing)
    error('Missing required columns: %s', strjoin(missing, ', '));
end

T = table();

% Subject as categorical
T.Subject = categorical(string(MasterZETA.SubjectNum));

% Fixed effects as categorical, with explicit reference levels
T.MoveType = categorical(string(MasterZETA.MoveType), ...
    {'REST','HAND OC','HAND PS','ARM EF'});  % update: HAND OC first = reference, no longer 'REST'
T.Depth = categorical(string(MasterZETA.Depth), {'b','c','t'}); % b as reference, no longer 't'

% Dep.Vars
T.IFR_mean_Hz = double(MasterZETA.IFR_mean_Hz);
T.logIFR_mean_Hz = log(T.IFR_mean_Hz + eps); % log transformed IFR_mean_Hz
T.IFR_norm = double(MasterZETA.IFR_norm);
T.IFR_baselineNorm = double(MasterZETA.IFR_baselineNorm);

% Drop missing DV rows (do separately per model later if you prefer)
T = T(~isnan(T.IFR_mean_Hz) & ~isnan(T.IFR_norm) & ~isnan(T.IFR_baselineNorm) & ~isnan(T.logIFR_mean_Hz), :);

end
