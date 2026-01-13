function T = build_IFR_LMM_Table(MasterZETA)

% Required vars
need = {'SubjectNum','MoveType','Depth', ...
    'IFR_mean_baselineNorm','IFR_max_baselineNorm', ...                     % exclude 'IFR_mean_Znorm', 'IFR_mean_Hz', 'IFR_max_Hz'...
    'IFR_peakLatency','IFR_peakOnset_Latency'};

missing = setdiff(need, MasterZETA.Properties.VariableNames);
if ~isempty(missing)
    error('Missing required columns: %s', strjoin(missing, ', '));
end

T = table();

% Random Effect: Subject as categorical
T.Subject = categorical(string(MasterZETA.SubjectNum));

% Fixed effects as categorical, with explicit reference levels
T.MoveType = categorical(string(MasterZETA.MoveType), ...
    {'REST','HAND OC','HAND PS','ARM EF'});                     % REST = reference,
T.Depth = categorical(string(MasterZETA.Depth), {'b','c','t'}); % b as reference, no longer 't'

% ---- Dependent variables ----
T.IFR_mean_Hz = double(MasterZETA.IFR_mean_Hz);
T.IFR_max_Hz = double(MasterZETA.IFR_max_Hz);

T.IFR_mean_Znorm = double(MasterZETA.IFR_mean_Znorm);
T.IFR_mean_baselineNorm = double(MasterZETA.IFR_mean_baselineNorm);
T.IFR_max_baselineNorm = double(MasterZETA.IFR_max_baselineNorm);

T.IFR_peakLatency       = double(MasterZETA.IFR_peakLatency);
T.IFR_peakOnset_Latency = double(MasterZETA.IFR_peakOnset_Latency);

% For now, minimal cleanup:
T = T(~isundefined(T.MoveType) & ~isundefined(T.Depth), :);

% Drop missing DV rows (do separately per model later if you prefer)
% T = T(~isnan(T.IFR_mean_Hz) & ~isnan(T.IFR_norm) & ~isnan(T.IFR_baselineNorm) & ~isnan(T.logIFR_mean_Hz), :);

end
