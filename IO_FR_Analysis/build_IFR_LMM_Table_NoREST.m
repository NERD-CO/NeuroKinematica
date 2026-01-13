function T = build_IFR_LMM_Table_NoREST(MasterZETA)

% Required vars
need = {'SubjectNum','MoveType','Depth', ...
    'IFR_mean_baselineNorm','IFR_max_baselineNorm', ...                     % exclude 'IFR_mean_Znorm', 'IFR_mean_Hz', 'IFR_max_Hz'...
    'IFR_peakLatency','IFR_peakOnset_Latency'};

missing = setdiff(need, MasterZETA.Properties.VariableNames);
if ~isempty(missing)
    error('Missing required columns: %s', strjoin(missing, ', '));
end

% Exclude REST
maskNoRest = string(strtrim(MasterZETA.MoveType)) ~= "REST";
M = MasterZETA(maskNoRest, :);

T = table();

% Random Effect: Subject as categorical
T.Subject = categorical(string(M.SubjectNum));

% Fixed effects as categorical, with explicit reference levels
T.MoveType = categorical(string(M.MoveType), ...
    {'HAND OC','HAND PS','ARM EF'});  % update: HAND OC first = reference, no longer 'REST'
T.Depth = categorical(string(M.Depth), {'b','c','t'}); % b as reference, no longer 't'

% ---- Dependent variables ----
T.IFR_mean_Hz = double(M.IFR_mean_Hz);
T.IFR_max_Hz = double(M.IFR_max_Hz);

T.IFR_mean_Znorm = double(M.IFR_mean_Znorm);
T.IFR_mean_baselineNorm = double(M.IFR_mean_baselineNorm);
T.IFR_max_baselineNorm = double(M.IFR_max_baselineNorm);

T.IFR_peakLatency       = double(M.IFR_peakLatency);
T.IFR_peakOnset_Latency = double(M.IFR_peakOnset_Latency);

% For now, minimal cleanup:
T = T(~isundefined(T.MoveType) & ~isundefined(T.Depth), :);

end
