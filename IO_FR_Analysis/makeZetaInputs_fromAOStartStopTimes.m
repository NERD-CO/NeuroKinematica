function [vecSpikeTimes_sec, matEventTimes_sec, useMaxDur_s, perTrialDur_s, kept_tbl] = ...
    makeZetaInputs_fromAOStartStopTimes(move_tbl, AO_spike_fs, varargin)

% Build ZETA inputs by stitching spike data trials with per-trial durations.
% Uses user-provided field names for spike data, start, and stop.
%
% OUTPUTS (for zetatest.m and getIFR)
%   vecSpikeTimes_sec   [S×1],  vector of stitched spike times (s)
%   matEventTimes_sec   [T×2],  event matrix [start stop] per trial (s) on stitched axis (true durations)
%   useMaxDur_s         scalar, ZETA analysis window (s) after onset
%   perTrialDur_s       [T×1],  true per-trial durations (s)
%   kept_tbl            table,  filtered trials actually used


% -------------------- varargin params --------------------
% Read name–value options with defaults:
p = inputParser;
p.addParameter('UseMaxDur_s', [],        @(x) isempty(x) || isscalar(x));   
p.addParameter('PadITI_s',     0.005,    @(x) isscalar(x) && x>=0);         % gap inserted between trials on the stitched timeline (prevents overlap)
p.addParameter('MinDur_s',     0,        @(x) isscalar(x) && x>=0);
p.addParameter('MaxDur_s',     inf,      @(x) isscalar(x) && x>0);
p.addParameter('SpikeField',   'C1',     @(s) ischar(s) || isstring(s));
p.addParameter('StartField',   'TTL_spk_idx_Start', @(s) ischar(s) || isstring(s));
p.addParameter('StopField',    'TTL_spk_idx_End',   @(s) ischar(s) || isstring(s)); % goes from start of current rep to start of next rep (minus 5 frames), a full move rep
% lead-in / lead-out to keep around each event (seconds)
p.addParameter('PreWindow_s',  0.050,      @(x) isscalar(x) && x>=0);
p.addParameter('PostWindow_s', 0.000,      @(x) isscalar(x) && x>=0);
p.parse(varargin{:}); % parse inputs
ParsedInputs = p.Results; % holds parsed inputs

Spkf = char(ParsedInputs.SpikeField);
Start_idx = char(ParsedInputs.StartField);
End_idx = char(ParsedInputs.StopField);

% -------------------- validate columns --------------------
needCols = {Spkf, Start_idx, End_idx};
assert(all(ismember(needCols, move_tbl.Properties.VariableNames)), ...
    'Table must contain columns: %s', strjoin(needCols, ', '));

% -------------------- keep usable trials --------------------
keep = ~cellfun(@isempty, move_tbl.(Spkf)) ...
     & ~isnan(move_tbl.(Start_idx)) ...
     & ~isnan(move_tbl.(End_idx)) ...
     & (move_tbl.(End_idx) > move_tbl.(Start_idx));

kept_tbl = move_tbl(keep,:);
nT = height(kept_tbl);
if nT==0
    vecSpikeTimes_sec = [];
    matEventTimes_sec = [];
    useMaxDur_s       = [];
    perTrialDur_s = [];
    return;
end

% -------------------- per-trial durations (s) --------------------
perTrialDur_s = (double(kept_tbl.(End_idx)) - double(kept_tbl.(Start_idx))) ./ AO_spike_fs;
perTrialDur_s = max(perTrialDur_s, ParsedInputs.MinDur_s);
perTrialDur_s = min(perTrialDur_s, ParsedInputs.MaxDur_s);
perTrialDur_s(~isfinite(perTrialDur_s)) = ParsedInputs.MinDur_s;
perTrialDur_s(perTrialDur_s <= 0) = max(1/AO_spike_fs, 1e-4);

% -------------------- choose analysis window --------------------
if isempty(ParsedInputs.UseMaxDur_s)
    useMaxDur_s = max(perTrialDur_s);
    if ~isfinite(useMaxDur_s) || useMaxDur_s <= 0
        useMaxDur_s = 0.4; % conservative fallback
    end
else
    useMaxDur_s = ParsedInputs.UseMaxDur_s;
end

% -------------------- stitch offsets --------------------
trialOffsets = (0:nT-1)' .* (useMaxDur_s + ParsedInputs.PadITI_s); % stitch trials on one timeline

% [T×2] event matrix on stitched axis (true start/stop; does NOT include pre/post windows)
eventOn_sec  = trialOffsets;
eventOff_sec = trialOffsets + perTrialDur_s;
matEventTimes_sec = [eventOn_sec, eventOff_sec];

% -------------------- spikes → stitched axis --------------------
% Keep spikes from (-PreWindow_s) before onset through (trialDur + PostWindow_s)
vecSpikeTimes_sec = [];
% Convert spike samples to seconds relative to that trial's onset.
for i = 1:nT
    spkSamp = double(kept_tbl.(Spkf){i}); % spikes in AO samples
    relSec  = (spkSamp - double(kept_tbl.(Start_idx)(i))) ./ AO_spike_fs; % each spike time relative to onset (s)

    % Decide which spikes to keep [−PreWindow_s​, perTrialDurs​(i) + PostWindow_s​]
    hi = perTrialDur_s(i) + ParsedInputs.PostWindow_s;  % upper bound = true duration + any post-extension
    keepSpk = (relSec >= -ParsedInputs.PreWindow_s) & (relSec <= hi);

    % Shift kept spikes forward by per-trial trialOffsets(i) to place each trial in correct stitched block
    if any(keepSpk)
        vecSpikeTimes_sec = [vecSpikeTimes_sec; relSec(keepSpk) + trialOffsets(i)];
    end
end

vecSpikeTimes_sec = sort(vecSpikeTimes_sec(:)); % return a single, sorted vector of stitched spike times in seconds
end
