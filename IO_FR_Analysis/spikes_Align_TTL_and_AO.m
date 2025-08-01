%% Align IO spike data segments with corresponding movement data segments

clear; clc;
addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis'

%% Environment / Directory Set-up

% specify directory where case-specific data files are located
curPCname = getenv('COMPUTERNAME');

switch curPCname
    case 'DESKTOP-I5CPDO7'  % PC_1
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'DSKTP-JTLAB-EMR'  % Lab Desktop
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'NSG-M-H8J3X34'    % PC_2
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
end

cd(IO_DataDir)
Subject_AO = readtable('Subject_AO.xlsx');

%% Hardcode Case-specific Data directories

% Isolate specific CaseDate / studyID (StudyNum in Subject_AO csv)
% CaseDate = '03_09_2023'; % studyID = 1, ptID 1

CaseDate = '03_23_2023'; % studyID = 2, ptID 2    * % Use for INS 2025
% CaseDate = '04_05_2023'; % studyID = 3, ptID 2    * % Use for INS 2025

% CaseDate = '04_13_2023_bilateral'; % studyID = 4(L), 5(R), ptID 3

% CaseDate = '05_11_2023'; % studyID = 6, ptID 4
%  CaseDate = '05_18_2023_a'; % studyID = 7, ptID 4

% CaseDate = '05_18_2023_b_bilateral'; % LSTN: studyID = 8, ptID = 5    % Use for INS 2025
% RSTN: studyID = 9, ptID = 5

% CaseDate = '05_31_2023';  % studyID = 10, ptID 6

% CaseDate = '06_08_2023_bilateral'; % studyID = 11(L), 12(R), ptID = 7

% CaseDate = '07_13_2023_bilateral'; % studyID = 15(L), 16(R), ptID = 9

% define case-specific data directory
Case_DataDir = [IO_DataDir, filesep, CaseDate];

% Define directories where case-specific IO ephys data are located
RawDataDir = [Case_DataDir, filesep, 'Raw Electrophysiology MATLAB'];       % directory where raw MATLAB data files are located (case-specific)
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology'];       % directory where processed MATLAB data should be saved (case-specific)


%% Handle bilateral cases and hemisphere selection

isBilateral = contains(CaseDate, 'bilateral', 'IgnoreCase', true);

if isBilateral
    fprintf('[INFO] Bilateral case detected: %s\n', CaseDate);
    
    % Prompt user for hemisphere choice (LSTN or RSTN)
    CaseDate_hem = input('Enter hemisphere (LSTN or RSTN): ', 's');
    
    % Validate input
    validHems = {'LSTN','RSTN'};
    if ~ismember(CaseDate_hem, validHems)
        error('Invalid input. Please enter LSTN or RSTN.');
    end
else
    CaseDate_hem = ''; % No hemisphere for unilateral cases
end

% Append hemisphere folder if needed
if ~isempty(CaseDate_hem)
    ProcDataDir = fullfile(ProcDataDir, CaseDate_hem);
    fprintf('[INFO] Hemisphere-specific directory set: %s\n', ProcDataDir);
else
    fprintf('[INFO] Using base processed directory: %s\n', ProcDataDir);
end

% Define Clustered Spike Times directory
ClustSpkTimesDir = fullfile(ProcDataDir, 'ClusteredSpikeTimes'); % directory where clustered spike times should be saved (case-specific)
if ~isfolder(ClustSpkTimesDir)
    error('[ERROR] ClustSpkTimesDir does not exist: %s', ClustSpkTimesDir);
end
cd(ClustSpkTimesDir);


%% Define case-specific kinematic data directory

MoveDataDir = [IO_DataDir, filesep, 'Processed DLC'];
cd(MoveDataDir)

% % Specify case ID
% Move_CaseID = 'IO_03_09_2023_RSTN'; % studyID = 1, ptID 1 (processed, incomplete case)

 Move_CaseID = 'IO_03_23_2023_LSTN'; % studyID = 2, ptID 2 (processed, complete case) *
% Move_CaseID = 'IO_04_05_2023_RSTN'; % studyID = 3, ptID 2 (processed, complete case) *

% Move_CaseID = 'IO_04_13_2023_LSTN'; % studyID = 4, ptID 3 (processed, complete case) 
% Move_CaseID = 'IO_04_13_2023_RSTN'; % studyID = 5, ptID 3

% Move_CaseID = 'IO_05_11_2023_LSTN'; % studyID = 6, ptID 4 (processed, incomplete case)
% Move_CaseID = 'IO_05_18_2023_a_RSTN'; % studyID = 7, ptID 4

% Move_CaseID = 'IO_05_18_2023_b_LSTN'; % studyID = 8, ptID 5 (processed, complete case) *
% Move_CaseID = 'IO_05_18_2023_b_RSTN'; % studyID = 9, ptID 5

% Move_CaseID = 'IO_05_31_2023_LSTN'; % studyID = 10, ptID 6

% Move_CaseID ='IO_06_08_2023_LSTN'; % studyID = 11, ptID = 7 (processed, complete case) 
% Move_CaseID ='IO_06_08_2023_RSTN'; % studyID = 12, ptID = 7 (processed, incomplete case)

% Move_CaseID ='IO_07_13_2023_LSTN'; % studyID = 15, ptID = 9  
% Move_CaseID ='IO_07_13_2023_RSTN'; % studyID = 16, ptID = 9 


%% Directory for movement indices

% isolate case-specific kinematic data directory
Move_CaseDir = [MoveDataDir, filesep, Move_CaseID];

% data subfolders:
Move_CaseMats = [Move_CaseDir, filesep, 'mat folder'];
Move_CaseVideos = [Move_CaseDir, filesep, 'video folder'];
cd(Move_CaseVideos)

% % for kinematic analyses
% cd(Move_CaseMats)
% moveMat = dir('*.mat');
% moveMat_names = {moveMat.name};


%% define datastream sampling rates (Alpha Omega and Video sampling fs)

TTL_fs = 44000; % Hz
AO_spike_fs = 44000; % Hz
AO_LFP_fs = 1375; % Hz
DLC_fs = 100; % fps


%% define offset duration

offset_ms = 50; % milliseconds
offset_seconds = offset_ms / 1000; % seconds

% Calculate the number of TTL samples
offset_TTLs = round(TTL_fs * offset_seconds); % ensure value is integer

% for future function input: useOffset; when = 1,  offset = offset; when = 0, offset = 0

%% go to ClustSpkTimesDir
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

    % Find row of ao_MAT_file that corresponds with trial
    SubjectAO_row = Subject_AO(contains(Subject_AO.ao_MAT_file,ProcName),:);

    switch SubjectAO_row.stn_loc{1}(1)
        case 'd'
            depthName = 't'
        case 'c'
            depthName = 'c'
        case 'v'
            depthName = 'b'
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
    trial_seconds = nan(height(SpkMoveTbl), 1);

    % loop through trials and pull out BeginF and EndF indices in MoveIndex csv per corresponding trial
    for move_i = 1:height(SpkMoveTbl)

        % define frame indices based on task recording context (frame indices in TTL_Down)                    % get these from movement index csv
        frame_startTime = SpkMoveTbl.BeginF(move_i); % index of recording initiation
        frame_endTime = SpkMoveTbl.EndF(move_i); % index of recording termination

        % %  if frame_endTime exceeds length of TTL down, trim SpkMoveTble
        % Check if frame indices are within bounds
        if frame_startTime > length(TTL_Down) || frame_endTime > length(TTL_Down)
            warning('Frame index out of bounds. Skipping trial %d.', move_i);
            continue; % Skip this trial if indices are invalid
        end

        % extract TTL signals of interest based on task context
        TTL_samp_taskStart = TTL_Down(frame_startTime); % number of samples wrt TTL clock
        TTL_samp_taskEnd = TTL_Down(frame_endTime);

        % define AO recording times
        AO_startTime = spikeClInfo.AOstartTS;

        % calculate difference in clock startTimes
        time_offset = TTL_clockStart - AO_startTime; % time (seconds) wrt TTL clock

        % convert time offset to sample offset
        sample_offset = round(time_offset*AO_spike_fs); % number of samples

        % convert TTL sample indices by the sample offset to transform into AO_spike clock / spike sample domain
        TTL_spk_idx_Start(move_i) = TTL_samp_taskStart + sample_offset; % number of samples wrt AO clock
        TTL_spk_idx_End = TTL_samp_taskEnd + sample_offset;

        % incorporate offset
        TTL_spk_idx_Start(move_i) =  TTL_spk_idx_Start(move_i) - offset_TTLs; % start [50 ms] before

        % check # of clusters in spike file
        if numel(unique(spikeClInfo.clusterIDS)) == 1
            % find spike cluster indices
            clust_time = spikeClInfo.SpikeTSindex; % samples wrt AO clock

            % determine spikes in cluster within movement block
            spikes_in_move1 = clust_time > TTL_spk_idx_Start(move_i) & clust_time < TTL_spk_idx_End; % logical indicating spikes w/in moveblock in AO_time

            % get spike times in cluster within movement block
            clustered_spikeTimes = clust_time(spikes_in_move1);

            SpkMoveTbl.('C1'){move_i} = clustered_spikeTimes;
            %((temp_spks - temp_row.TTL_spk_idx_Start)/AO_spike_fs) - 0.05)

            trial_seconds = (clustered_spikeTimes - TTL_spk_idx_Start(move_i)/AO_spike_fs) - offset_seconds; % incorporate this in All_SpikesPerMove_Tbl
            SpkMoveTbl.('C1_ts'){move_i} = trial_seconds;

        else
            % determine cluster IDs
            clust_IDs = unique(spikeClInfo.clusterIDS);

            for cii = 1:length(clust_IDs)
                % find unique cluster ID indices
                clust_index = ismember(spikeClInfo.clusterIDS,clust_IDs(cii));
                clust_time = spikeClInfo.SpikeTSindex(clust_index); % samples wrt AO clock

                % determine spikes in distinct clusters within movement block
                spikes_in_move1 = clust_time > TTL_spk_idx_Start(move_i) & clust_time < TTL_spk_idx_End; % logical indicating spikes w/in moveblock in AO_time

                % get spike times in distinct clusters within movement block
                clustered_spikeTimes = clust_time(spikes_in_move1);

                SpkMoveTbl.(['C', num2str(cii)]){move_i} = clustered_spikeTimes;

                trial_seconds = (clustered_spikeTimes - TTL_spk_idx_Start(move_i)/AO_spike_fs) - offset_seconds; % incorporate this in All_SpikesPerMove_Tbl
                SpkMoveTbl.(['C', num2str(cii), '_ts']){move_i} = trial_seconds;
            end
        end
    end

    % Add spike_ID, move_ID, and TTL_spk_idx_Start columns
    SpkMoveTbl.spike_trial_ID = spike_trial_ID;
    SpkMoveTbl.move_trial_ID = move_trial_ID;
    SpkMoveTbl.TTL_spk_idx_Start = TTL_spk_idx_Start;

    All_spikesPermove{spk_mat_name} = SpkMoveTbl;

end


% Define the standard column order
standard_col_order = {'MoveN', 'MoveType', 'BeginF', 'EndF', 'TTL_spk_idx_Start', 'spike_trial_ID', 'move_trial_ID'}; % add 'trial_seconds'

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

% save All_SpikesPerMove_Tbl to a file
SpikesPerMove_Dir = [Case_DataDir, filesep, 'DLC_Ephys'];

% % For Bialteral cases:
 SpikesPerMove_Dir = [Case_DataDir, filesep, 'DLC_Ephys', filesep, CaseDate_hem]; % comment out when N/A

cd(SpikesPerMove_Dir)

save('All_SpikesPerMove_offset.mat',"All_SpikesPerMove_Tbl");


%% next steps

% movement indices from raw case vids
% 3/09/2023     - incomplete case contexts, revisit (Ruby)
% 3/23/2023     - complete, awesome case :)
% 4/05/2023     - partially complete
% 4/13/2023_L   - complete, odd results, revisit
% 4/13/2023_R   - incomplete case contexts
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

