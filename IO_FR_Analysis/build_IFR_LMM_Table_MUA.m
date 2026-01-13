function T = build_IFR_LMM_Table_MUA(MasterZETA_MUA)
% build_IFR_LMM_Table_MUA
% Builds the LMM input table for MUA IFR outcomes from MasterZETA_MUA.

% ---- Required vars (match SU pipeline) ----
need = {'SubjectNum','MoveType','Depth', ...
    'IFR_mean_baselineNorm','IFR_max_baselineNorm'};                        % exclude: 'IFR_mean_Znorm', 'IFR_mean_Hz', 'IFR_max_Hz'

missing = setdiff(need, MasterZETA_MUA.Properties.VariableNames);
if ~isempty(missing)
    error('Missing required columns in MasterZETA_MUA: %s', strjoin(missing, ', '));
end

T = table();

% Random effect: Subject as categorical
T.Subject = categorical(string(MasterZETA_MUA.SubjectNum));

% Fixed effects as categorical with explicit reference levels
% (Reference is the FIRST level in the list.)
T.MoveType = categorical(string(MasterZETA_MUA.MoveType), ...
    {'REST','HAND OC','HAND PS','ARM EF'});  % REST reference (matches contrasts function)
T.Depth = categorical(string(MasterZETA_MUA.Depth), {'b','c','t'}); % ventral b reference

% ---- Dependent variables ----
T.IFR_mean_Hz       = double(MasterZETA_MUA.IFR_mean_Hz);
T.IFR_max_Hz        = double(MasterZETA_MUA.IFR_max_Hz);
T.IFR_mean_Znorm    = double(MasterZETA_MUA.IFR_mean_Znorm);
T.IFR_mean_baselineNorm = double(MasterZETA_MUA.IFR_mean_baselineNorm);
T.IFR_max_baselineNorm  = double(MasterZETA_MUA.IFR_max_baselineNorm);

% Optional (if you later decide to include baseline Hz as an additional DV)
if ismember('IFR_baseline_Hz', MasterZETA_MUA.Properties.VariableNames)
    T.IFR_baseline_Hz = double(MasterZETA_MUA.IFR_baseline_Hz);
else
    T.IFR_baseline_Hz = nan(height(T),1);
end

% ---- Row filtering ----
% For now, minimal cleanup:
T = T(~isundefined(T.MoveType) & ~isundefined(T.Depth), :);

end