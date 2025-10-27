function [vecSpikeTimes_sec, matEventTimes_sec, useMaxDur_s, perTrialStimDur_s, kept_tbl] = makeZetaInputs_fromAOStartStopTimes(move_tbl, AO_spike_fs, varargin)

% Helper function to compute per-trial durations from TTL_spk_idx_Start/TTL_spk_idx_Stop 
% Returns [T×2] onset/offset to build ZETA inputs by stitching trials, using TRUE per-trial stop samples.
%
% REQUIRED move_tbl columns:
%   - C1                  : cell array; AO spike sample indices per trial
%   - TTL_spk_idx_Start   : scalar AO sample (trial onset) per row
%   - TTL_spk_idx_Stop    : scalar AO sample (trial end)   per row
%
% Name-Value (optional):
%   'UseMaxDur_s' : analysis window after each onset (sec). If [], uses max(perTrialStimDur_s). Default [].
%   'PadITI_s'    : pad between stitched trials (sec). Default 0.005.
%   'MinStimDur_s': clamp minimum duration (sec). Default 0.
%   'MaxStimDur_s': clamp maximum duration (sec). Default inf. %%%% adjust
%
% OUTPUTS:
%   vecSpikeTimes_sec   : [S×1] stitched spike times (sec)
%   matEventTimes_sec   : [T×2] [on off] per trial (sec) on stitched axis
%   useMaxDur_s         : scalar ZETA window (sec)
%   perTrialStimDur_s   : [T×1] durations actually used (sec)
%   kept_tbl            : filtered move_tbl rows that contributed

% ---- parse args
p = inputParser;
p.addParameter('UseMaxDur_s', [], @(x) isempty(x) || isscalar(x));
p.addParameter('PadITI_s', 0.005, @(x) isscalar(x) && x>=0);
p.addParameter('MinStimDur_s', 0, @(x) isscalar(x) && x>=0);
p.addParameter('MaxStimDur_s', inf, @(x) isscalar(x) && x>0);
p.parse(varargin{:});
U = p.Results;

% ---- validate columns
needCols = {'C1','TTL_spk_idx_Start','TTL_spk_idx_Stop'};
assert(all(ismember(needCols, move_tbl.Properties.VariableNames)), ...
    'Table must contain columns: %s', strjoin(needCols, ', '));

% ---- keep usable trials
keep = ~cellfun(@isempty, move_tbl.C1) ...
     & ~isnan(move_tbl.TTL_spk_idx_Start) ...
     & ~isnan(move_tbl.TTL_spk_idx_Stop) ...
     & (move_tbl.TTL_spk_idx_Stop > move_tbl.TTL_spk_idx_Start);
kept_tbl = move_tbl(keep,:);
nT = height(kept_tbl);
if nT==0
    vecSpikeTimes_sec = []; matEventTimes_sec = []; useMaxDur_s = []; perTrialStimDur_s = []; return;
end

% ---- per-trial durations (sec) from AO start/stop
perTrialStimDur_s = (double(kept_tbl.TTL_spk_idx_Stop) - double(kept_tbl.TTL_spk_idx_Start)) ./ AO_spike_fs;

% clamp/sanitize
perTrialStimDur_s = max(perTrialStimDur_s, U.MinStimDur_s);
perTrialStimDur_s = min(perTrialStimDur_s, U.MaxStimDur_s);
perTrialStimDur_s(~isfinite(perTrialStimDur_s)) = U.MinStimDur_s;
perTrialStimDur_s(perTrialStimDur_s<=0) = max(1/AO_spike_fs, 1e-4);

% ---- choose useMaxDur_s (spacing/analysis window)
if isempty(U.UseMaxDur_s)
    useMaxDur_s = max(perTrialStimDur_s);              % safe default; ZETA ignores spikes > useMaxDur_s
    if ~isfinite(useMaxDur_s) || useMaxDur_s<=0
        useMaxDur_s = 0.4; % conservative fallback
    end
else
    useMaxDur_s = U.UseMaxDur_s;
end

% ---- stitch offsets
trialOffsets = (0:nT-1)' .* (useMaxDur_s + U.PadITI_s);

% [T×2] event matrix on stitched axis
eventOn_sec  = trialOffsets;
eventOff_sec = trialOffsets + perTrialStimDur_s;
matEventTimes_sec = [eventOn_sec, eventOff_sec];

% ---- convert spikes to stitched axis; clip to [0,useMaxDur_s] per trial for ZETA window
vecSpikeTimes_sec = [];
for i = 1:nT
    spkSamp = double(kept_tbl.C1{i});
    relSec  = (spkSamp - double(kept_tbl.TTL_spk_idx_Start(i))) ./ AO_spike_fs;
    keepSpk = relSec >= 0 & relSec <= useMaxDur_s;
    if any(keepSpk)
        vecSpikeTimes_sec = [vecSpikeTimes_sec; relSec(keepSpk) + trialOffsets(i)]; %#ok<AGROW>
    end
end
vecSpikeTimes_sec = sort(vecSpikeTimes_sec(:));
end
