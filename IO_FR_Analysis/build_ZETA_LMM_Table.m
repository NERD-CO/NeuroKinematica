function T = build_ZETA_LMM_Table(MasterZETA)

% Required vars
need = {'SubjectNum','MoveType','Depth','ZetaZ', ...
        'ZETA_peakLatency','ZETA_invSign_Latency', 'dblZetaP'};

missing = setdiff(need, MasterZETA.Properties.VariableNames);
if ~isempty(missing)
    error('Missing required columns: %s', strjoin(missing, ', '));
end

T = table();

% Random effect: Subject as categorical
T.Subject = categorical(string(MasterZETA.SubjectNum));

% Fixed effects as categorical with explicit reference levels: first = ref
T.MoveType = categorical(string(MasterZETA.MoveType), {'REST','HAND OC','HAND PS','ARM EF'});  % ref: 'REST'
T.Depth = categorical(string(MasterZETA.Depth), {'b','c','t'}); % ref: b (ventral)

% Dependent variables (DVs, SU)
T.ZetaZ = double(MasterZETA.ZetaZ);
T.ZetaD = double(MasterZETA.ZetaD);

T.ZETA_peakLatency     = double(MasterZETA.ZETA_peakLatency);
T.ZETA_invSign_Latency = double(MasterZETA.ZETA_invSign_Latency);

% Useful for filtering or later reporting (not a DV)
T.ZetaP = double(MasterZETA.dblZetaP);

% minimal cleanup:
T = T(~isundefined(T.MoveType) & ~isundefined(T.Depth), :);

end
