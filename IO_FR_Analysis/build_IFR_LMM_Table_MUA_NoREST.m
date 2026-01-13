function T = build_IFR_LMM_Table_MUA_NoREST(MasterZETA_MUA)

% build_IFR_LMM_Table_MUA_NoREST
% Build LMM-ready table for MUA IFR outcomes, excluding REST.
%
% Reference levels:
%   MoveType ref = HAND OC
%   Depth ref    = b

need = {'SubjectNum','MoveType','Depth', ...
    'IFR_mean_baselineNorm','IFR_max_baselineNorm'};                        % exclude: 'IFR_mean_Znorm', 'IFR_mean_Hz', 'IFR_max_Hz'

missing = setdiff(need, MasterZETA_MUA.Properties.VariableNames);
if ~isempty(missing)
    error('Missing required columns in MasterZETA_MUA: %s', strjoin(missing, ', '));
end

% Exclude REST
maskNoRest = string(strtrim(MasterZETA_MUA.MoveType)) ~= "REST";
M = MasterZETA_MUA(maskNoRest, :);

T = table();

% Random effect: Subject as categorical
T.Subject = categorical(string(M.SubjectNum));

% Fixed effects (explicit levels so coefficient names match planned-contrast code)
T.MoveType = categorical(string(M.MoveType), {'HAND OC','HAND PS','ARM EF'}); % HAND OC ref
T.Depth    = categorical(string(M.Depth),    {'b','c','t'});                  % b ref

% DVs
T.IFR_mean_Hz      = double(M.IFR_mean_Hz);
T.IFR_max_Hz = double(M.IFR_max_Hz);
T.IFR_mean_Znorm = double(M.IFR_mean_Znorm);
T.IFR_mean_baselineNorm = double(M.IFR_mean_baselineNorm);
T.IFR_max_baselineNorm = double(M.IFR_max_baselineNorm);

% ---- Row filtering ----
% For now, minimal cleanup:
T = T(~isundefined(T.MoveType) & ~isundefined(T.Depth), :);

end
