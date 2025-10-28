function [vecSpikeTimes_sec, matEventTimes_sec, useMaxDur_s, perTrialStimDur_s, kept_tbl] = ...
    makeZetaInputs_fromAOStartStopTimes(move_tbl, AO_spike_fs, varargin)

% Build ZETA inputs by stitching trials with TRUE per-trial durations.
% Uses user-provided field names for spikes, start, and stop.
%
% OUTPUTS
%   vecSpikeTimes_sec   [S×1]  stitched spike times (s)
%   matEventTimes_sec   [T×2]  [on off] per trial (s) on stitched axis (true durations)
%   useMaxDur_s         scalar ZETA analysis window (s) after onset
%   perTrialStimDur_s   [T×1]  true per-trial durations (s)
%   kept_tbl            table  filtered trials actually used

% -------------------- params --------------------
p = inputParser;
p.addParameter('UseMaxDur_s', [],        @(x) isempty(x) || isscalar(x));
p.addParameter('PadITI_s',     0.005,    @(x) isscalar(x) && x>=0);
p.addParameter('MinStimDur_s', 0,        @(x) isscalar(x) && x>=0);
p.addParameter('MaxStimDur_s', inf,      @(x) isscalar(x) && x>0);
p.addParameter('SpikeField',   'C1',     @(s) ischar(s) || isstring(s));
p.addParameter('StartField',   'TTL_spk_idx_Start', @(s) ischar(s) || isstring(s));
p.addParameter('StopField',    'TTL_spk_idx_End',   @(s) ischar(s) || isstring(s));
% NEW: optional lead-in / lead-out to keep around each event (seconds)
p.addParameter('PreWindow_s',  0.0,      @(x) isscalar(x) && x>=0);
p.addParameter('PostWindow_s', 0.0,      @(x) isscalar(x) && x>=0);
p.parse(varargin{:});
U = p.Results;

Sf = char(U.SpikeField);
St = char(U.StartField);
Ed = char(U.StopField);

% -------------------- validate columns --------------------
needCols = {Sf, St, Ed};
assert(all(ismember(needCols, move_tbl.Properties.VariableNames)), ...
    'Table must contain columns: %s', strjoin(needCols, ', '));

% -------------------- keep usable trials --------------------
keep = ~cellfun(@isempty, move_tbl.(Sf)) ...
     & ~isnan(move_tbl.(St)) ...
     & ~isnan(move_tbl.(Ed)) ...
     & (move_tbl.(Ed) > move_tbl.(St));

kept_tbl = move_tbl(keep,:);
nT = height(kept_tbl);
if nT==0
    vecSpikeTimes_sec = [];
    matEventTimes_sec = [];
    useMaxDur_s       = [];
    perTrialStimDur_s = [];
    return;
end

% -------------------- per-trial true durations (s) --------------------
perTrialStimDur_s = (double(kept_tbl.(Ed)) - double(kept_tbl.(St))) ./ AO_spike_fs;
perTrialStimDur_s = max(perTrialStimDur_s, U.MinStimDur_s);
perTrialStimDur_s = min(perTrialStimDur_s, U.MaxStimDur_s);
perTrialStimDur_s(~isfinite(perTrialStimDur_s)) = U.MinStimDur_s;
perTrialStimDur_s(perTrialStimDur_s <= 0)       = max(1/AO_spike_fs, 1e-4);

% -------------------- choose analysis window --------------------
if isempty(U.UseMaxDur_s)
    useMaxDur_s = max(perTrialStimDur_s);
    if ~isfinite(useMaxDur_s) || useMaxDur_s <= 0
        useMaxDur_s = 0.4; % conservative fallback
    end
else
    useMaxDur_s = U.UseMaxDur_s;
end

% -------------------- stitch offsets --------------------
trialOffsets = (0:nT-1)' .* (useMaxDur_s + U.PadITI_s);

% [T×2] event matrix on stitched axis (true on/off; does NOT include pre/post windows)
eventOn_sec  = trialOffsets;
eventOff_sec = trialOffsets + perTrialStimDur_s;
matEventTimes_sec = [eventOn_sec, eventOff_sec];

% -------------------- spikes → stitched axis --------------------
% Keep spikes from (-PreWindow_s) before onset through (StimDur + PostWindow_s)
vecSpikeTimes_sec = [];
for i = 1:nT
    spkSamp = double(kept_tbl.(Sf){i});                         % spikes in AO samples
    relSec  = (spkSamp - double(kept_tbl.(St)(i))) ./ AO_spike_fs; % relative to onset (s)

    hi      = perTrialStimDur_s(i) + U.PostWindow_s;
    keepSpk = (relSec >= -U.PreWindow_s) & (relSec <= hi);

    if any(keepSpk)
        vecSpikeTimes_sec = [vecSpikeTimes_sec; relSec(keepSpk) + trialOffsets(i)]; %#ok<AGROW>
    end
end

vecSpikeTimes_sec = sort(vecSpikeTimes_sec(:));
end
