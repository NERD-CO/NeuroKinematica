function [X, t_use, mtPerUnit, labelsUnit] = build_MUA_ZETA_timeAlignedMatrix(M, varargin)
% build_MUA_ZETA_timeAlignedMatrix
%
% Aligns MUA temporal deviation vectors (vecD_MUA) to a common time axis
% (new_vecTime_MUA) and returns a matrix X [time x units] suitable for PCA.
%
% INPUT
%   M : subtable of MasterZETA_MUA for a single depth, containing:
%       - MoveType        (categorical / string / char)
%       - vecD_MUA        (cell, each = numeric vector)
%       - new_vecTime_MUA (cell, each = numeric vector, same size as vecD_MUA)
%
% Optional Name–Value:
%   'DepthCode' : used only for warnings (e.g., 't','c','b')
%
% OUTPUT
%   X         : [time x nUnits] matrix of aligned temporal deviation values
%   t_use     : column vector of time points (seconds)
%   mtPerUnit : MoveType label (string) per column of X
%   labelsUnit: unit / channel label (best effort; may be generic)

p = inputParser;
p.addParameter('DepthCode','', @(x)ischar(x) || isstring(x));
p.parse(varargin{:});
depthCode = char(p.Results.DepthCode);

% ---- Find a reference time axis ----
hasTime = ~cellfun(@isempty, M.new_vecTime_MUA);
hasDev  = ~cellfun(@isempty, M.vecD_MUA);
idxRef  = find(hasTime & hasDev, 1, 'first');

if isempty(idxRef)
    error('No non-empty vecD_MUA + new_vecTime_MUA rows for depth "%s".', depthCode);
end

time_ref = M.new_vecTime_MUA{idxRef};
if isempty(time_ref) || ~isnumeric(time_ref)
    error('Reference new_vecTime_MUA is empty or non-numeric at depth "%s".', depthCode);
end
time_ref = time_ref(:);   % column

nUnits   = height(M);
X        = nan(numel(time_ref), nUnits);
keepUnit = false(1, nUnits);
mtPerUnit = strings(nUnits,1);
labelsUnit = strings(nUnits,1);

% Try to use a nice label column if present
varNames = M.Properties.VariableNames;
if ismember('PrettyLabel', varNames)
    labelVar = 'PrettyLabel';
elseif ismember('ChannelLabel', varNames)
    labelVar = 'ChannelLabel';
else
    labelVar = '';
end

for u = 1:nUnits
    dVec = M.vecD_MUA{u};
    tVec = M.new_vecTime_MUA{u};

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

        if ~isempty(labelVar)
            labelsUnit(u) = string(M.(labelVar)(u));
        else
            labelsUnit(u) = "MUA_" + u;
        end
    catch ME
        warning('interp1 failed for MUA unit %d at depth %s: %s', ...
            u, depthCode, ME.message);
    end
end

% Keep only good units
X         = X(:, keepUnit);
mtPerUnit = mtPerUnit(keepUnit);
labelsUnit = labelsUnit(keepUnit);

if size(X,2) < 1
    error('No valid MUA units after alignment at depth "%s".', depthCode);
end

% Remove timepoints that are NaN across all units
goodTime = any(isfinite(X),2);
t_use    = time_ref(goodTime);
X        = X(goodTime,:);

% Fill remaining NaNs with unit-wise means
for j = 1:size(X,2)
    nanMask = ~isfinite(X(:,j));
    if any(nanMask)
        mu_j = mean(X(~nanMask,j),'omitnan');
        X(nanMask,j) = mu_j;
    end
end
end
