function [vecSpikeTimes_sec, matEventTimes_sec, useMaxDur_s] = makeZetaInputsFromTbl(move_tbl, AO_spike_fs, varargin)
% Build ZETA inputs by stitching trials onto a continuous time axis.
%
% INPUTS
%   move_tbl         : table with columns:
%                      - C1: cell array of spike sample indices (AO units)
%                      - TTL_spk_idx_Start: scalar start index (AO units) per trial
%                      - (optional) TTL_spk_idx_End or segment end samples, if available
%   AO_spike_fs      : spike sampling rate (e.g., 44000)
% VARARGIN (Name-Value):
%   'UseMaxDur_s'    : override useMaxDur (seconds). If empty, computed as
%                      robust max segment length across trials.
%   'StimDur_s'      : if provided, we output [T×2] event on/off; otherwise [T×1].
%   'PadITI_s'       : pad between stitched trials (default 0.005 s)
%
% OUTPUTS
%   vecSpikeTimes_sec   : [S×1] spike times (seconds), stitched continuous time axis
%   matEventTimes_sec   : [T×1] or [T×2] event on(/off) times (seconds) in same axis
%   useMaxDur_s         : scalar window used for ZETA

% ---- parse args
p = inputParser;
p.addParameter('UseMaxDur_s', [], @(x) isempty(x) || isscalar(x));
p.addParameter('StimDur_s',   [], @(x) isempty(x) || isscalar(x));
p.addParameter('PadITI_s', 0.005, @(x) isscalar(x) && x>=0);
p.parse(varargin{:});
UseMaxDur_s = p.Results.UseMaxDur_s;
StimDur_s   = p.Results.StimDur_s;
PadITI_s    = p.Results.PadITI_s;

% ---- keep only valid trials with at least one spike in-window
nTrials = height(move_tbl);
hasSpikes = false(nTrials,1);
segLen_s   = nan(nTrials,1);

for i = 1:nTrials
    spk = move_tbl.C1{i};
    if ~isempty(spk) && numel(spk) >= 1 && ~isnan(move_tbl.TTL_spk_idx_Start(i))
        hasSpikes(i) = true;
        % estimate segment duration from max spike index relative to start
        maxRelSamp = max(spk - move_tbl.TTL_spk_idx_Start(i));
        if ~isempty(maxRelSamp) && isfinite(maxRelSamp) && maxRelSamp>0
            segLen_s(i) = maxRelSamp / AO_spike_fs;
        else
            segLen_s(i) = 0;
        end
    end
end
move_tbl = move_tbl(hasSpikes,:);
nTrials  = height(move_tbl);
if nTrials==0
    vecSpikeTimes_sec = [];
    matEventTimes_sec = [];
    useMaxDur_s = [];
    return
end

% ---- decide useMaxDur
if isempty(UseMaxDur_s)
    % robust max across trials (ignore empty/NaN)
    useMaxDur_s = max(segLen_s(isfinite(segLen_s)));
    if isempty(useMaxDur_s) || useMaxDur_s<=0
        % fallback to 0.4 s if nothing usable is found
        useMaxDur_s = 0.4;
    end
else
    useMaxDur_s = UseMaxDur_s;
end

% ---- stitch trials:
% each trial gets a base offset so events are spaced by (useMaxDur_s + PadITI_s)
trialOffsets = (0:nTrials-1)' .* (useMaxDur_s + PadITI_s);

% ---- event times (onsets, and optional offsets)
eventOn_sec  = trialOffsets;                         % onsets start each stitched trial
if isempty(StimDur_s)
    matEventTimes_sec = eventOn_sec;                 % [T×1]
else
    eventOff_sec = eventOn_sec + StimDur_s;          % [T×2] if you want mean-rate stats
    matEventTimes_sec = [eventOn_sec, eventOff_sec];
end

% ---- spikes: convert each trial's spike samples to (seconds since onset),
% then add that trial's offset to stitch
vecSpikeTimes_sec = [];
for i = 1:nTrials
    spkSamp = move_tbl.C1{i};
    relSamp = spkSamp - move_tbl.TTL_spk_idx_Start(i);
    relSec  = double(relSamp) ./ AO_spike_fs;
    % keep only spikes within [0, useMaxDur_s]
    keep    = relSec >= 0 & relSec <= useMaxDur_s;
    relSec  = relSec(keep);
    % stitch
    vecSpikeTimes_sec = [vecSpikeTimes_sec; relSec + trialOffsets(i)]; %#ok<AGROW>
end

% ensure sorted (ZETA sorts anyway, but good hygiene)
vecSpikeTimes_sec = sort(vecSpikeTimes_sec(:));
end
