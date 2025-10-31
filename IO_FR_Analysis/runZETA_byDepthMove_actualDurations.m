function [ZETA_Summary, all_sZETA, all_sRate, all_sLatencies] = runZETA_byDepthMove_actualDurations( ...
    All_SpikesPerMove_Tbl, AO_spike_fs, varargin)

% Run ZETA for each MoveType × STN depth using true per-trial durations.
% Returns a summary table and per-subset sZETA/sRate/sLatencies collections.
%
% OUTPUTS (for zetatest.m and getIFR)
% ZETA_Summary: one row per (MoveType × Depth) with key metrics,
% all_sZETA:    the full sZETA structs (one per row),
% all_sRate:    the IFR/peak metadata structs (one per row),
% all_sLatencies: the latency structs (one per row).


% -------------------- varargin params --------------------
% Read name–value options with defaults:
p = inputParser;
p.addParameter('UseMaxDur_s', [], @(x) isempty(x) || isscalar(x));
p.addParameter('PadITI_s', 0.005, @(x) isscalar(x) && x>=0);
p.addParameter('Resamples', 2000, @(x) isscalar(x) && x>0);
p.addParameter('PlotFlag', 0, @(x) isscalar(x) && ismember(x,0:4));
p.addParameter('RestrictRange', [-inf inf], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('DirectQuantile', false, @islogical);
p.addParameter('JitterSize', 2, @(x) isscalar(x) && x>0);
p.addParameter('Stitch', true, @islogical);
p.addParameter('DepthIDs', {'t','c','b'});
p.addParameter('MoveTypeOrder', {'HAND OC','HAND PS','ARM EF','REST'});
% table defaults:
% schema: keep backward-compat 'SpikeField' but prefer 'SpikeFields'
p.addParameter('SpikeField', 'C1', @(s) ischar(s) || isstring(s));   % if single electrode {'C1'}
p.addParameter('SpikeFields', [], @(x) iscellstr(x) || isempty(x));  % if multiple electrodes {'C1','C2',...}
% start/stop bounds for spike field segments:
p.addParameter('StartField', 'TTL_spk_idx_Start', @(s) ischar(s) || isstring(s));
p.addParameter('StopField',  'TTL_spk_idx_End',   @(s) ischar(s) || isstring(s)); %%% 
% optional lead-in/out (sec) padding:
p.addParameter('PreWindow_s',  0.050, @(x) isscalar(x) && x>=0);  % 50 ms pre-trial offset
p.addParameter('PostWindow_s', 0.000, @(x) isscalar(x) && x>=0);
p.parse(varargin{:});
ParsedInputs = p.Results; % holds parsed inputs


% Normalize spike fields: prefer SpikeFields; else use SpikeField; else auto-detect C\d+ with any spikes
spikeFields = ParsedInputs.SpikeFields;
if isempty(spikeFields)
    if ~isempty(ParsedInputs.SpikeField)
        spikeFields = {char(ParsedInputs.SpikeField)};
    else
        cand = All_SpikesPerMove_Tbl.Properties.VariableNames;
        isCnum = ~cellfun('isempty', regexp(cand,'^C\d+$','once'));
        spikeFields = cand(isCnum);
        % keep only C# columns that are cell arrays with any non-empty cell
        keep = false(size(spikeFields));
        for k = 1:numel(spikeFields)
            col = All_SpikesPerMove_Tbl.(spikeFields{k});
            keep(k) = iscell(col) && any(~cellfun(@isempty,col));
        end
        spikeFields = spikeFields(keep);
    end
end


move_types = intersect(ParsedInputs.MoveTypeOrder, unique(All_SpikesPerMove_Tbl.MoveType),'stable');

rows = {};            % will hold one cell row per (MoveType × Depth) to later convert into a table
all_sZETA = {};       % cell array aligned to rows of ZETA_Summary
all_sRate = {};
all_sLatencies = {};

% Outer loop over spike fields (C#s) 
for SpkF = spikeFields
    curSF = SpkF{1};

    % Loop over MoveType × STN depth
    for moveT_i = 1:numel(move_types)
        for depth_i = 1:numel(ParsedInputs.DepthIDs)
            mv = move_types{moveT_i};
            dz = ParsedInputs.DepthIDs{depth_i};

            move_tbl = All_SpikesPerMove_Tbl( ...
                strcmp(All_SpikesPerMove_Tbl.MoveType, mv) & ...
                contains(All_SpikesPerMove_Tbl.move_trial_ID, dz), :);

            if isempty(move_tbl), continue; end

            % Run helper function: makeZetaInputs_fromAOStartStopTimes
            % Build stitched inputs with true [start stop] and optional pre/post spikes
            [spkT, evTimes, useMaxDur, trialDur, kept_tbl] = makeZetaInputs_fromAOStartStopTimes( ...
                move_tbl, AO_spike_fs, ...
                'UseMaxDur_s', ParsedInputs.UseMaxDur_s, ...
                'PadITI_s',    ParsedInputs.PadITI_s, ...
                'SpikeField',  curSF, ...
                'StartField',  ParsedInputs.StartField, ...
                'StopField',   ParsedInputs.StopField, ...
                'PreWindow_s', ParsedInputs.PreWindow_s, ...
                'PostWindow_s',ParsedInputs.PostWindow_s);

            if isempty(spkT) || isempty(evTimes), continue; end

            % Run ZETA
            % (https://github.com/JorritMontijn/zetatest/blob/main/zetatest.m)
            [pZETA, sZETA, sRate, sLat] = zetatest( ...
                spkT, evTimes, useMaxDur, ...
                ParsedInputs.Resamples, ParsedInputs.PlotFlag, ParsedInputs.RestrictRange, ...
                ParsedInputs.DirectQuantile, ParsedInputs.JitterSize, ParsedInputs.Stitch);

            % Pull optional mean-rate stats (field names differ by version)
            MeanD = getFieldOr(sZETA,'dblMeanD', NaN);   % doc says dblMeanD
            MeanZ = getFieldOr(sZETA,'dblMeanZ', NaN);   % code sometimes provides dblMeanZ
            MeanP = getFieldOr(sZETA,'dblMeanP', NaN);

            % Collect in table row (include big vectors as cell entries)
            rows(end+1,:) = { ...
                curSF, mv, dz, height(kept_tbl), useMaxDur, ...
                pZETA, sZETA.dblZETA, sZETA.dblD, sZETA.dblP, ...
                sZETA.dblZetaT, sZETA.dblD_InvSign, sZETA.dblZetaT_InvSign, ...
                getFieldOr(sRate,'dblPeakTime',nan), getFieldOr(sRate,'dblOnset',nan), ...
                mean(trialDur,'omitnan'), std(trialDur,'omitnan'), ...
                MeanD, MeanZ, MeanP, ...
                sZETA.vecSpikeT(:)', ...            % spike timestamps at ZETA sampling grid
                sZETA.dblUseMaxDur, ...             % ZETA's internal record of window length
                sZETA.vecLatencies(:)'};            % different latency estimates, number determined by intLatencyPeaks. If no peaks are detected, it returns NaNs:

            % Keep full structures for saving
            all_sZETA{end+1,1} = sZETA;
            all_sRate{end+1,1} = sRate;
            all_sLatencies{end+1,1} = sLat;
        end
    end
end


varNames = { ...
    'SpikeField','MoveType','Depth','nTrials','UseMaxDur_s', ...
    'dblZetaP','ZetaZ','ZetaD','ZetaP', ...
    'ZetaTime','ZetaD_Inv','ZetaTime_Inv', ...
    'IFR_PeakTime','IFR_OnsetTime', ...
    'MeanStimDur_s','StdStimDur_s', ...
    'MeanD','MeanZ','MeanP', ...
    'VecSpikeT','ZetaUseMaxDur','VecLatencies'};

if isempty(rows)
    ZETA_Summary = cell2table(cell(0,numel(varNames)), 'VariableNames', varNames);
else
    ZETA_Summary = cell2table(rows, 'VariableNames', varNames);
end
end

% ---------- utils ----------
function v = getFieldOr(S, fld, defaultVal)
if isempty(S) || ~isfield(S,fld) || isempty(S.(fld)), v = defaultVal; else, v = S.(fld); end
end
