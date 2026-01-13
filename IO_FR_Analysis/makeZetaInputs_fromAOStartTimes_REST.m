function [vecSpikeTimes_sec, matEventTimes_sec, useMaxDur_s, perTrialDur_s, kept_tbl] = ...
    makeZetaInputs_fromAOStartTimes_REST(move_tbl, AO_spike_fs, varargin)

% REST-only stitcher for getIFR / zetatest inputs.
%
% Key behavior:
% - Effective trial start = StartField + StartOffset_s (default 0.250s)
% - Trial duration >= MinDur_s (default 0.750s)
% - PreWindow_s default 0 (no pre-onset)
% - DropEmptySpikeTrials default false (keep trials w/ zero spikes)
%
% Optional: if EndField exists and yields a longer duration than MinDur_s,
%          use it (after applying start offset).
% Modes:
%   DropEmptySpikeTrials = true  -> "unit-centric" (only trials where this SpikeField fired)
%   DropEmptySpikeTrials = false -> "rep-centric"  (all valid reps; empty-spike trials kept)
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
p.addParameter('PadITI_s',     0.005,    @(x) isscalar(x) && x>=0);         % 5 ms gap inserted between trials on the stitched timeline (prevents overlap)
p.addParameter('SpikeField',   'C1',     @(s) ischar(s) || isstring(s));

% Start field name (in AO samples)
p.addParameter('StartField',   'TTL_spk_idx_Start', @(s) ischar(s) || isstring(s)); 

% REST-specific timing
p.addParameter('StartOffset_s', 0.250,   @(x) isscalar(x) && x>=0);         % add 250 ms post-onset (to buffer out noise)
p.addParameter('MinDur_s',      0.750,   @(x) isscalar(x) && x>=0);         % enforces a min duration of 750 ms after shifted start
p.addParameter('MaxDur_s',      1.250,   @(x) isscalar(x) && x>0);          % enforces a max duration of 1.25 s after shifted start

% Optional end field name (in AO samples)
p.addParameter('EndField',    'TTL_spk_idx_End',   @(s) ischar(s) || isstring(s));

% lead-in / lead-out to keep around each event (seconds)
p.addParameter('PreWindow_s',  0.000,      @(x) isscalar(x) && x>=0);       % remove 50 ms pre-onset for REST
p.addParameter('PostWindow_s', 0.000,      @(x) isscalar(x) && x>=0);

% REST default: keep empty-spike trials
p.addParameter('DropEmptySpikeTrials', false, @(x) isscalar(x) && (islogical(x) || isnumeric(x))); % do NOT drop trials drop for REST

p.parse(varargin{:}); % parse inputs
U = p.Results; % holds parsed inputs


Spkf      = char(U.SpikeField);
Start_idx = char(U.StartField);
End_idx   = char(U.EndField);


% -------------------- validate columns --------------------
needCols = {Spkf, Start_idx};
assert(all(ismember(needCols, move_tbl.Properties.VariableNames)), ...
    'Table must contain columns: %s', strjoin(needCols, ', '));

hasEndField = ~isempty(End_idx) && ismember(End_idx, move_tbl.Properties.VariableNames);


% -------------------- keep usable trials --------------------
hasStart = ~isnan(move_tbl.(Start_idx));

% Spike presence definition (robust to non-cell / missing)
hasSpikes = true(height(move_tbl),1);
if ismember(Spkf, move_tbl.Properties.VariableNames)
    col = move_tbl.(Spkf);
    if iscell(col)
        hasSpikes = ~cellfun(@isempty, col);
    else
        % if non-cell, treat non-missing/non-empty as "has spikes"
        try
            hasSpikes = ~ismissing(col);
        catch
            hasSpikes = true(height(move_tbl),1);
        end
    end
end

if logical(U.DropEmptySpikeTrials)
    keep = hasStart & hasSpikes; % unit centric
else
    keep = hasStart; % rep centric, keep all REST trials with a valid start
end

kept_tbl = move_tbl(keep,:);
nT = height(kept_tbl);

if nT==0
    vecSpikeTimes_sec = [];
    matEventTimes_sec = [];
    useMaxDur_s       = [];
    perTrialDur_s     = [];
    return;
end

% -------------------- build per-trial effective start (AO samples) --------------------
startSamp_raw = double(kept_tbl.(Start_idx));
startSamp_eff = startSamp_raw + round(U.StartOffset_s * AO_spike_fs);

% -------------------- per-trial durations (s) --------------------
% perTrialDur_s = (double(kept_tbl.(End_idx)) - double(kept_tbl.(Start_idx))) ./ AO_spike_fs;
% perTrialDur_s = max(perTrialDur_s, U.MinDur_s);
% perTrialDur_s = min(perTrialDur_s, U.MaxDur_s);

% Default duration = MinDur_s
perTrialDur_s = repmat(double(U.MinDur_s), nT, 1);

% If EndField exists, use it when it yields 
 
% If EndField exists, use it when it yields a valid duration AFTER offset
% longer than MinDur_s, shorter than MaxDur_s:
if hasEndField
    endSamp_raw = double(kept_tbl.(End_idx));

    % duration after shifted start
    dur_s = (endSamp_raw - startSamp_eff) ./ AO_spike_fs;

    % valid if finite and end is after shifted start
    ok = isfinite(dur_s) & (dur_s > 0);

    % clamp to [MinDur_s, MaxDur_s]
    dur_s(ok) = max(dur_s(ok), double(U.MinDur_s));
    dur_s(ok) = min(dur_s(ok), double(U.MaxDur_s));

    % apply only where ok; otherwise keep default MinDur_s
    perTrialDur_s(ok) = dur_s(ok);
end


% Ensure finite + positive
perTrialDur_s(~isfinite(perTrialDur_s)) = double(U.MinDur_s);
perTrialDur_s(perTrialDur_s <= 0) = max(1/AO_spike_fs, 1e-4);


% -------------------- choose analysis window --------------------
if isempty(U.UseMaxDur_s)
    useMaxDur_s = max(perTrialDur_s);
    if ~isfinite(useMaxDur_s) || useMaxDur_s <= 0
        useMaxDur_s = double(U.MinDur_s);
    end
else
    useMaxDur_s = double(U.UseMaxDur_s);
end

% -------------------- stitch offsets --------------------
trialOffsets = (0:nT-1)' .* (useMaxDur_s + U.PadITI_s);

% Event matrix (true start/stop; does not include pre/post)
eventOn_sec  = trialOffsets;
eventOff_sec = trialOffsets + perTrialDur_s;
matEventTimes_sec = [eventOn_sec, eventOff_sec];

% -------------------- spikes → stitched axis --------------------
vecSpikeTimes_sec = [];

for i = 1:nT
    spkCell = kept_tbl.(Spkf){i};
    if isempty(spkCell)
        continue; % keep trial, just contributes no spikes
    end

    spkSamp = double(spkCell); % spikes in AO samples
    % relSec  = (spkSamp - double(kept_tbl.(Start_idx)(i))) ./ AO_spike_fs;
    relSec  = (spkSamp - startSamp_eff(i)) ./ AO_spike_fs;

    hi = perTrialDur_s(i) + U.PostWindow_s;
    keepSpk = (relSec >= -U.PreWindow_s) & (relSec <= hi);

    if any(keepSpk)
        vecSpikeTimes_sec = [vecSpikeTimes_sec; relSec(keepSpk) + trialOffsets(i)];
    end
end

vecSpikeTimes_sec = sort(vecSpikeTimes_sec(:)); % return a single, sorted vector of stitched spike times in seconds

end
