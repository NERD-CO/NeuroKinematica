function [X, t_use, mtPerUnit, labelsUnit] = build_SU_ZETA_timeAlignedMatrix(M, varargin)

% build_SU_ZETA_timeAlignedMatrix
%
% For single-unit ZETA:
%   Takes a subset of MasterZETA rows (M) that all share a common
%   PSTH_TimeCenters_s grid (e.g., one depth or depth×MoveType subset),
%   and returns a time-aligned data matrix X:
%
%       X      : [time × units] ZETA deviation values (ZETA_vecD interpolated)
%       t_use  : [time × 1] time axis used after removing all-NaN rows
%       mtPerUnit  : [units × 1] MoveType (string) for each column in X
%       labelsUnit : [units × 1] PrettyLabel for each column in X
%
% Usage example:
%   [X, t_use, mtPerUnit, labelsUnit] = build_SU_ZETA_timeAlignedMatrix(catRows);
%
% Optional Name–Value:
%   'DepthCode'  – for nicer warning messages (e.g. 't','c','b')

p = inputParser;
p.addParameter('DepthCode', '', @(s) isstring(s) || ischar(s));
p.parse(varargin{:});
depthCode = string(p.Results.DepthCode);

% Sanity: required columns
reqVars = {'MoveType','ZETA_vecD','ZETA_vecT','PSTH_TimeCenters_s','PrettyLabel'};
missing = setdiff(reqVars, M.Properties.VariableNames);
if ~isempty(missing)
    error('build_SU_ZETA_timeAlignedMatrix: missing columns: %s', strjoin(missing,', '));
end

if isempty(M)
    X = [];
    t_use = [];
    mtPerUnit = strings(0,1);
    labelsUnit = strings(0,1);
    return;
end

%% --- Reference time axis from first valid row ---

time_ref = M.PSTH_TimeCenters_s{1};
if isempty(time_ref) || ~isnumeric(time_ref)
    error('build_SU_ZETA_timeAlignedMatrix: Reference PSTH_TimeCenters_s is empty or non-numeric.');
end
time_ref = time_ref(:);   % column vector

nUnits    = height(M);
X         = nan(numel(time_ref), nUnits);
keepUnit  = false(1, nUnits);
mtPerUnit = strings(nUnits,1);      % MoveType per unit
labelsUnit = strings(nUnits,1);     % PrettyLabel per unit

%% --- Loop units and interpolate onto time_ref ---

for u = 1:nUnits
    dVec = M.ZETA_vecD{u};
    tVec = M.ZETA_vecT{u};

    if isempty(dVec) || isempty(tVec)
        continue;
    end

    dVec = double(dVec(:));
    tVec = double(tVec(:));

    if numel(dVec) < 3 || numel(tVec) ~= numel(dVec)
        continue;
    end

    try
        X(:,u) = interp1(tVec, dVec, time_ref, 'linear', 'extrap');
        keepUnit(u)  = true;
        mtPerUnit(u) = string(M.MoveType(u));
        labelsUnit(u)= string(M.PrettyLabel(u));
    catch ME
        warning('build_SU_ZETA_timeAlignedMatrix: interp1 failed for unit %d at depth %s: %s', ...
            u, depthCode, ME.message);
    end
end

%% --- Keep only good units ---

X          = X(:, keepUnit);
mtPerUnit  = mtPerUnit(keepUnit);
labelsUnit = labelsUnit(keepUnit);

if isempty(X)
    t_use = time_ref(:);
    return;
end

%% --- Remove timepoints that are NaN across all units ---

goodTime = any(isfinite(X), 2);
X        = X(goodTime, :);
t_use    = time_ref(goodTime);

%% --- Fill remaining NaNs with unit-wise means ---

for j = 1:size(X,2)
    nanMask = ~isfinite(X(:,j));
    if any(nanMask)
        mu_j = mean(X(~nanMask,j), 'omitnan');
        X(nanMask,j) = mu_j;
    end
end

end
