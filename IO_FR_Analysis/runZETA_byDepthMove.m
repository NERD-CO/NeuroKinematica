function ZETA_Summary = runZETA_byDepthMove(All_SpikesPerMove_Tbl, AO_spike_fs, varargin)

% Quickly run ZETA for each MoveType × STN depth ('t','c','b')

p = inputParser;
p.addParameter('UseMaxDur_s', [], @(x) isempty(x) || isscalar(x)); % optional fixed window
p.addParameter('StimDur_s',   [], @(x) isempty(x) || isscalar(x)); % optional fixed stim duration
p.addParameter('Resamples', 1000, @(x) isscalar(x) && x>0);
p.addParameter('PlotFlag', 0, @(x) isscalar(x) && ismember(x,0:4));
p.addParameter('RestrictRange', [-inf inf], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('DirectQuantile', false, @islogical);
p.addParameter('JitterSize', 2, @(x) isscalar(x) && x>0);
p.addParameter('Stitch', true, @islogical);
p.parse(varargin{:});
u = p.Results;

depth_ids    = {'t','c','b'};
move_types   = intersect({'HAND OC','HAND PS','ARM EF','REST'}, unique(All_SpikesPerMove_Tbl.MoveType),'stable');

rows = [];
for m = 1:numel(move_types)
    for d = 1:numel(depth_ids)
        move_type = move_types{m};
        depth_id  = depth_ids{d};
        % subset one neuron's trials at this depth & move
        move_tbl = All_SpikesPerMove_Tbl(strcmp(All_SpikesPerMove_Tbl.MoveType, move_type) & ...
                                         contains(All_SpikesPerMove_Tbl.move_trial_ID, depth_id), :);
        if isempty(move_tbl); continue; end

        % build ZETA inputs
        [vecSpikeTimes_sec, matEventTimes_sec, useMaxDur_s] = ...
            makeZetaInputsFromTbl(move_tbl, AO_spike_fs, ...
                                  'UseMaxDur_s', u.UseMaxDur_s, ...
                                  'StimDur_s',   u.StimDur_s, ...
                                  'PadITI_s',    0.005);

        if isempty(vecSpikeTimes_sec) || numel(matEventTimes_sec)==0
            continue
        end

        % call ZETA
        [pZ, sZ, sRate, sLat] = zetatest( ...
            vecSpikeTimes_sec, ...
            matEventTimes_sec, ...
            useMaxDur_s, ...
            u.Resamples, ...
            u.PlotFlag, ...
            u.RestrictRange, ...
            u.DirectQuantile, ...
            u.JitterSize, ...
            u.Stitch);

        % collect
        rows = [rows; {move_type, depth_id, useMaxDur_s, pZ, sZ.dblZETA, sZ.dblD, sZ.dblP, ...
                       sZ.dblZetaT, sZ.dblD_InvSign, sZ.dblZetaT_InvSign, ...
                       getfield(sRate,'dblPeakTime',nan), getfield(sRate,'dblOnset',nan)}]; %#ok<GFLD,AGROW>
    end
end

ZETA_Summary = cell2table(rows, 'VariableNames', { ...
    'MoveType','Depth','UseMaxDur_s','ZetaP','ZetaZ','ZetaD','ZetaP_again', ...
    'ZetaTime','ZetaD_Inv','ZetaTime_Inv','IFR_PeakTime','IFR_OnsetTime'});
end
