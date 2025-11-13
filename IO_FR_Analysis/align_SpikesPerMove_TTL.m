function All_SpikesPerMove_Tbl = align_SpikesPerMove_TTL(Subject_AO, AO_spike_fs, TTL_fs, ProcDataDir, ClustSpkTimesDir, Move_CaseDir, pre_offset_ms, useOffset, SpikesPerMove_Dir)

% Align IO spike data segments with corresponding movement data segments per trial


%% Define case-specific directory for movement indices per trial

% Move_CaseDir data subfolders:
Move_CaseMats = [Move_CaseDir, filesep, 'mat folder'];      % contains processed DLC timeseries data (csv-to-mat)
Move_CaseVideos = [Move_CaseDir, filesep, 'video folder'];  % contains DLC-labeled videos and Movement Index CSVs
cd(Move_CaseVideos)

% % for kinematic analyses
% cd(Move_CaseMats)
% moveMat = dir('*.mat');
% moveMat_names = {moveMat.name};


%% define offset duration

offset_seconds = pre_offset_ms / 1000; % seconds

% Calculate the number of TTL samples
offset_TTLs = round(TTL_fs * offset_seconds); % ensure value is integer

% for future function input: useOffset; when = 1,  offset = offset; when = 0, offset = 0

%% Run useOffset_spikes function

[offset_spike_samples, meta_Offset_spk] = useOffset_spikes(TTL_fs, AO_spike_fs, pre_offset_ms, useOffset);


%% Main function

% go to ClustSpkTimesDir
cd(ClustSpkTimesDir)

% list of filename
SPKmatfiles = dir('*.mat');
SPKmatnames = {SPKmatfiles.name};

%All_spikesPermove = {length(SPKmatnames)};
All_spikesPermove = cell(length(SPKmatnames), 1);

for spk_mat_name = 1:length(SPKmatnames)

    cd(ClustSpkTimesDir)
    load(SPKmatnames{spk_mat_name},'spikeClInfo')

    fileparts = split(SPKmatnames{spk_mat_name},'_');
    ProcName = fileparts{3};

    cd(ProcDataDir)
    matfiles = dir('*.mat');
    matnames = {matfiles.name};
    ProcFile = matnames{contains(matnames, ProcName)};
    load(ProcFile, 'ProcEphys')

    % isolate fields of interest
    TTL_Down = ProcEphys.TTL.Down; % TTL signal down voltage deflection     % 1 cell represents frame#, values within the cell represent sample#
    TTL_clockStart = ProcEphys.TTL.startTime; % start-time of TTL clock     (in seconds wrt AO system start)

    % number of TTL frames (works for cell or numeric array)
    nFramesTTL = numel(TTL_Down);

    % Find row of ao_MAT_file that corresponds with trial
    SubjectAO_row = Subject_AO(contains(Subject_AO.ao_MAT_file,ProcName),:);

    switch SubjectAO_row.stn_loc{1}(1)
        case 'd'
            depthName = 't';
        case 'c'
            depthName = 'c';
        case 'v'
            depthName = 'b';
    end

    motor_trial_ID = [depthName, num2str(SubjectAO_row.trialNum),'_', SubjectAO_row.depth{1}];

    % Generate list of Motor Index CSVs
    cd(Move_CaseVideos)
    moveCSV = dir('*.csv');
    moveCSV_names = {moveCSV.name};

    % Generate list of Motor Index CSVs (filters for CSVs that contain 'Move' string)
    moveCSV = moveCSV_names(contains(moveCSV_names,'Move'));

    if ~any(contains(moveCSV, motor_trial_ID))  % checking if logical is false
        continue
    end

    moveTbl_name = moveCSV{contains(moveCSV, motor_trial_ID)};
    moveTbl = readtable(moveTbl_name);
    SpkMoveTbl = moveTbl(cellfun(@(X) ~isempty(X), moveTbl.MoveType, 'UniformOutput', true),:); % clean moveTbl (remove zeros)

    % initialize arrays
    spike_trial_ID = repmat({ProcFile}, height(SpkMoveTbl), 1);
    move_trial_ID = repmat({motor_trial_ID}, height(SpkMoveTbl), 1);
    TTL_spk_idx_Start = nan(height(SpkMoveTbl), 1);
    TTL_spk_idx_End = nan(height(SpkMoveTbl), 1);
    trial_seconds = nan(height(SpkMoveTbl), 1);


    % loop through trials and pull out BeginF and EndF indices in MoveIndex csv per corresponding trial
    for move_i = 1:height(SpkMoveTbl)

        % define frame indices based on task recording context (frame indices in TTL_Down)
        % ---------- get these from movement index csv ----------
        frame_startTime = SpkMoveTbl.BeginF(move_i); % frame index of move initiation (open, pronated, extended)
        frame_endTime1 = SpkMoveTbl.EndF(move_i); % frame index of move termination (closed, supinated, flexed)

        % -----------------------------------------------------------

        % redefine frame_endTime dynamically based on next moveRep BeginF
        % ---------- choose frame_endTime by MoveType ----------

        frameLead_buffer = 10;  % frames before next movement onset
        isLastRow = (move_i == height(SpkMoveTbl));

        % Clamp start into valid range, always enforce start+1
        frame_startTime = min(max(1, frame_startTime), max(1, nFramesTTL-1));

        mvAll = string(SpkMoveTbl.MoveType);
        mv = mvAll(move_i);
        activeSet = ["HAND OC","HAND PS","ARM EF"];

        if any(mv == activeSet)
            % Find the next row of the SAME MoveType
            rowIdx = (1:height(SpkMoveTbl))';
            nextIdx = find(mvAll == mv & rowIdx > move_i, 1, 'first');

            if ~isempty(nextIdx) && isfinite(SpkMoveTbl.BeginF(nextIdx))
                % Candidate 1: next BeginF - buffer (same MoveType)
                frame_endTime = round(SpkMoveTbl.BeginF(nextIdx)) - frameLead_buffer;
            else
                % Last repetition of this MoveType -> Candidate 2: this row's EndF
                frame_endTime = frame_endTime1;
            end
        else
            % REST (or any other non-active type) -> Candidate 2: this row's EndF
            frame_endTime = frame_endTime1;
        end

        % Fallback if EndF is missing
        if ~isfinite(frame_endTime) || isnan(frame_endTime) || isLastRow
            frame_endTime = nFramesTTL;
        end

        % Final clamps: at least one frame long, within TTL range, integer
        frame_endTime = max(frame_endTime, frame_startTime + 1);
        frame_endTime = min(frame_endTime, nFramesTTL);
        frame_endTime = round(frame_endTime);

        % -------------------------------------------------------------------------

        % mv = string(SpkMoveTbl.MoveType(move_i));
        %
        % if any(mv == ["HAND OC","HAND PS","ARM EF"])
        %     % Candidate 1: next BeginF - buffer (if there is a next row)
        %     if ~isnan(SpkMoveTbl.BeginF(move_i + 1))
        %         frame_endTime = SpkMoveTbl.BeginF(move_i + 1) - frameLead_buffer;
        %     else
        %         % If the next BeginF is missing, fall back to last TTL frame
        %         frame_endTime = nFramesTTL;
        %     end
        %
        % elseif mv == "REST"
        %     % Candidate 2: this row's EndF
        %     frame_endTime = frame_endTime1;
        %
        % elseif isLastRow
        %     % Candidate 3: last TTL frame for the last row
        %     frame_endTime = nFramesTTL;
        %
        % else
        %     % Fallback for any other MoveType: behave like REST (use EndF)
        %     frame_endTime = frame_endTime1;
        % end
        %
        % % Safety clamps: finite, within [start+1, nFramesTTL], and integer
        % if ~isfinite(frame_endTime) || isnan(frame_endTime)
        %     frame_endTime = nFramesTTL;
        % end
        % % Ensure at least one frame long
        % frame_endTime = max(frame_endTime, frame_startTime + 1);
        % % Don't exceed last TTL frame
        % frame_endTime = min(frame_endTime, nFramesTTL);
        % frame_endTime = round(frame_endTime);  % frames should be integer indices

        % -----------------------------------------------------------

        % % Ensure ≥1 frame long
        % if frame_endTime <= frame_startTime
        %     frame_endTime = min(frame_startTime + 1, nFramesTTL);
        % end

        % --------------------------------------------------------------

        % --- fetch TTL sample indices for those frames (cell- or numeric-safe) ---
        if iscell(TTL_Down)
            TTL_samp_taskStart = TTL_Down{frame_startTime};  % sample index (integer)
            TTL_samp_taskEnd   = TTL_Down{frame_endTime};    % sample index (integer)
        else
            TTL_samp_taskStart = TTL_Down(frame_startTime);
            TTL_samp_taskEnd   = TTL_Down(frame_endTime);
        end
        % --------------------------------------------------------------


        % define AO recording times
        AO_startTime = spikeClInfo.AOstartTS;

        % calculate difference in clock startTimes
        time_offset = TTL_clockStart - AO_startTime; % time (seconds) wrt TTL clock

        % convert time offset to sample offset
        sample_offset = round(time_offset*AO_spike_fs); % number of samples

        % convert TTL sample indices by the sample offset to transform into AO_spike clock / spike sample domain
        TTL_spk_idx_Start(move_i) = TTL_samp_taskStart + sample_offset; % number of samples wrt AO clock
        TTL_spk_idx_End(move_i) = TTL_samp_taskEnd + sample_offset;

        % incorporate offset
        % TTL_spk_idx_Start(move_i) =  TTL_spk_idx_Start(move_i) - offset_TTLs; % start [50 ms] before
        TTL_spk_idx_Start(move_i) =  TTL_spk_idx_Start(move_i) - offset_spike_samples;

        % check # of clusters in spike file
        if numel(unique(spikeClInfo.clusterIDS)) == 1
            % find spike cluster indices
            clust_time = spikeClInfo.SpikeTSindex; % samples wrt AO clock

            % determine spikes in cluster within movement block
            spikes_in_move1 = clust_time > TTL_spk_idx_Start(move_i) & clust_time < TTL_spk_idx_End(move_i); % logical indicating spikes w/in moveblock in AO_time

            % get spike times in cluster within movement block
            clustered_spikeTimes = clust_time(spikes_in_move1);

            SpkMoveTbl.('C1'){move_i} = clustered_spikeTimes;
            %((temp_spks - temp_row.TTL_spk_idx_Start)/AO_spike_fs) - 0.05)

            trial_seconds = ((clustered_spikeTimes - TTL_spk_idx_Start(move_i))/AO_spike_fs) - offset_seconds; % incorporate this in All_SpikesPerMove_Tbl
            SpkMoveTbl.('C1_ts'){move_i} = trial_seconds;

        else
            % determine cluster IDs
            clust_IDs = unique(spikeClInfo.clusterIDS);

            for cii = 1:length(clust_IDs)
                % find unique cluster ID indices
                clust_index = ismember(spikeClInfo.clusterIDS,clust_IDs(cii));
                clust_time = spikeClInfo.SpikeTSindex(clust_index); % samples wrt AO clock

                % determine spikes in distinct clusters within movement block
                spikes_in_move1 = clust_time > TTL_spk_idx_Start(move_i) & clust_time < TTL_spk_idx_End(move_i); % logical indicating spikes w/in moveblock in AO_time

                % get spike times in distinct clusters within movement block
                clustered_spikeTimes = clust_time(spikes_in_move1);

                SpkMoveTbl.(['C', num2str(cii)]){move_i} = clustered_spikeTimes;

                trial_seconds = ((clustered_spikeTimes - TTL_spk_idx_Start(move_i))/AO_spike_fs) - offset_seconds; % incorporate this in All_SpikesPerMove_Tbl
                SpkMoveTbl.(['C', num2str(cii), '_ts']){move_i} = trial_seconds;
            end
        end
    end

    % Add spike_ID, move_ID, and TTL_spk_idx_Start columns
    SpkMoveTbl.spike_trial_ID = spike_trial_ID;
    SpkMoveTbl.move_trial_ID = move_trial_ID;
    SpkMoveTbl.TTL_spk_idx_Start = TTL_spk_idx_Start;
    SpkMoveTbl.TTL_spk_idx_End = TTL_spk_idx_End;


    All_spikesPermove{spk_mat_name} = SpkMoveTbl;

end


% Define the standard column order
standard_col_order = {'MoveN', 'MoveType', 'BeginF', 'EndF', 'TTL_spk_idx_Start', 'TTL_spk_idx_End', 'spike_trial_ID', 'move_trial_ID'}; % add 'trial_seconds'

% Get the unique cluster IDs from all tables
cluster_ids = {};
clust_ts = {};
for tbl_i = 1:length(All_spikesPermove)
    tbl_1 = All_spikesPermove{tbl_i};

    % Check if tbl_1 is empty or not a table
    if isempty(tbl_1)
        continue
    end
    if ~istable(tbl_1)
        error('Entry %d in All_spikesPermove is not a table. Type: %s', tbl_i, class(tbl_1));
    end

    % Access table properties
    col_names = tbl_1.Properties.VariableNames;
    cluster_ids = [cluster_ids, col_names(startsWith(col_names, 'C'))];
    %clust_ts = []
end
cluster_ids = unique(cluster_ids);


% Combine standard columns with cluster IDs to create the full standard order
standard_col_order = [standard_col_order, cluster_ids];

% Initialize All_data as a cell array
All_data = {};
% Loop through each table in All_moveTbl_array
for tbl_i = 1:length(All_spikesPermove)
    tbl_1 = All_spikesPermove{tbl_i};
    if isempty(tbl_1)
        continue
    end

    % Initialize temp_data with NaNs or empty cells
    temp_data = cell(height(tbl_1), length(standard_col_order));

    % Align each table's data to the standard column order
    for col_idx = 1:length(standard_col_order)
        col_name = standard_col_order{col_idx};
        if ismember(col_name, tbl_1.Properties.VariableNames)
            temp_data(:, col_idx) = table2cell(tbl_1(:, col_name));
        end
    end

    % Append the aligned data to All_data
    All_data = [All_data; temp_data];
end

% Convert All_data to a table with the standard column order
All_SpikesPerMove_Tbl = cell2table(All_data, 'VariableNames', standard_col_order);

cd(SpikesPerMove_Dir)
save('All_SpikesPerMove_offset.mat',"All_SpikesPerMove_Tbl");


%% next steps

% movement indices from raw case vids
% 3/09/2023     - incomplete case contexts, revisit (Ruby)
% 3/23/2023     - complete, awesome case
% 4/05/2023     - complete, awesome case
% 4/13/2023_L   - complete, odd results, revisit
% 5/11/2023     - incomplete case contexts, revisit (Ruby)
% ...


% firing rate cut off for STN spike clusters
% bin size
% raster plots
% PSTHs
% baseline norm

% use this as template for LFP analysis

%% initial exploratory code

lastSpikeTS = spikeClInfo.SpikeTSindex(end);
numSpikes = spikeClInfo.SpikeTSindex(end) - spikeClInfo.SpikeTSindex(1);
num_seconds = numSpikes/AO_spike_fs;

% waveforms of spikes only necessary for figures
waves = spikeClInfo.waveForms.waves;
figure;
plot(waves)

