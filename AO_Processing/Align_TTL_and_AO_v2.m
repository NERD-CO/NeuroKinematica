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

% isolate a specific CaseDate / studyID (StudyNum in Subject_AO csv)
CaseDate = '03_09_2023'; % studyID = 1
% CaseDate = '03_23_2023'; % studyID = 2

% define case-specific data directory
Case_DataDir = [IO_DataDir, filesep, CaseDate];

% directories where case-specific IO ephys data are located
RawDataDir = [Case_DataDir, filesep, 'Raw Electrophysiology MATLAB'];       % directory where raw MATLAB data files are located (case-specific)
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology'];       % directory where processed MATLAB data should be saved (case-specific)
ClustSpkTimesDir = [ProcDataDir, filesep, 'ClusteredSpikeTimes'];           % directory where clustered spike times should be saved (case-specific)


%% directory for movement indices

% define kinematic data directory
MoveDataDir = [IO_DataDir, filesep, 'Kinematic Analyses'];

% specify case ID
Move_CaseID = 'IO_03_09_2023_RSTN'; % studyID = 1
% Move_CaseID = 'IO_03_23_2023_LSTN'; % studyID = 2

% isolate case-specific kinematic data directory
Move_CaseDir = [MoveDataDir, filesep, Move_CaseID];

% data subfolders:
Move_CaseMats = [Move_CaseDir, filesep, 'mat folder'];
Move_CaseVideos = [Move_CaseDir, filesep, 'video folder'];

% isolate frames of movement type
% moveType = 'Hand O/C'


cd(Move_CaseMats)
moveMat = dir('*.mat');
moveMat_names = {moveMat.name};




%% define datastream sampling rates (Alpha Omega and Video sampling fs)

TTL_fs = 44000; % Hz
AO_spike_fs = 44000; % Hz
AO_LFP_fs = 1375; % Hz
DLC_fs = 100; % fps


%% go to ClustSpkTimesDir
cd(ClustSpkTimesDir)

% list of filename
SPKmatfiles = dir('*.mat');
SPKmatnames = {SPKmatfiles.name};

% All_moveTbl = table(); % table containing spike times per motor trial
All_moveTbl_array = {length(SPKmatnames)};

for spk_mat_names = 1:length(SPKmatnames)
    cd(ClustSpkTimesDir)

    load(SPKmatnames{spk_mat_names},'spikeClInfo')

    fileparts = split(SPKmatnames{spk_mat_names},'_');
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
        case 'v'
            depthName = 'b'
        otherwise
            depthName = 'c'
    end

    move_trial_ID = [depthName, num2str(SubjectAO_row.trialNum),'_', SubjectAO_row.depth{1}];

    % Generate list of Motor Index CSVs
    cd(Move_CaseVideos)
    moveCSV = dir('*.csv');
    moveCSV_names = {moveCSV.name};

    % Generate list of Motor Index CSVs (filters for CSVs that contain 'Move' string)
    moveCSV = moveCSV_names(contains(moveCSV_names,'Move'));

    if ~any(contains(moveCSV, move_trial_ID))  % checking if logical is false
        continue
    end

    moveTbl_name = moveCSV{contains(moveCSV, move_trial_ID)};

    moveTbl = readtable(moveTbl_name);
    moveTbl = moveTbl(cellfun(@(X) ~isempty(X), moveTbl.MoveType, 'UniformOutput', true),:); % clean moveTbl (remove zeros)

    spike_ID = repmat({ProcFile}, height(moveTbl), 1);
    move_ID = repmat({move_trial_ID}, height(moveTbl), 1);
   

    % loop through trials and pull out BeginF and EndF indices in MoveIndex csv per corresponding trial
    for move_i = 1:height(moveTbl)

        % define frame indices based on task recording context (frame indices in TTL_Down)                    % get these from movement index csv
        frame_startTime = moveTbl.BeginF(move_i); % index of recording initiation
        frame_endTime = moveTbl.EndF(move_i); % index of recording termination

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
        TTL_spk_idx_Start = TTL_samp_taskStart + sample_offset; % number of samples wrt AO clock
        TTL_spk_idx_End = TTL_samp_taskEnd + sample_offset;

        % check # of clusters in spike file
        if numel(unique(spikeClInfo.clusterIDS)) == 1

            % find spike cluster indices
            clust_time = spikeClInfo.SpikeTSindex; % samples wrt AO clock

            % determine spikes in cluster within movement block
            spikes_in_move1 = clust_time > TTL_spk_idx_Start & clust_time < TTL_spk_idx_End; % logical indicating spikes w/in moveblock in AO_time

            % get spike times in cluster within movement block
            clustered_spikeTimes = clust_time(spikes_in_move1);

            moveTbl.('C1'){move_i} = clustered_spikeTimes;

        else

            % determine clust IDs
            clust_IDs = unique(spikeClInfo.clusterIDS);

            for cii = 1:length(clust_IDs)
                % find unique cluster ID indices
                clust_index = ismember(spikeClInfo.clusterIDS,clust_IDs(cii));
                clust_time = spikeClInfo.SpikeTSindex(clust_index); % samples wrt AO clock

                % determine spikes in distinct clusters within movement block
                spikes_in_move1 = clust_time > TTL_spk_idx_Start & clust_time < TTL_spk_idx_End; % logical indicating spikes w/in moveblock in AO_time

                % get spike times in distinct clusters within movement block
                clustered_spikeTimes = clust_time(spikes_in_move1);

                moveTbl.(['C', num2str(cii)]){move_i} = clustered_spikeTimes;
            end

        end

    end

    All_moveTbl_array{spk_mat_names} = moveTbl;

    % moveTbl.spike_ID = spike_ID;
    % moveTbl.move_ID = move_ID;
    % All_moveTbl = [All_moveTbl; moveTbl];

end

% create a table containing =# rows in 
% spike_trail_name
% move_trial_name

%% next steps
% understand this code
% use this as template for LFP analysis
% movement indices from raw case vids
% 3/09/2023
% 3/23/2023

% code to add in movement indices

%% movement block structure
% define start frame of event
% define end frame of event
% define offset duration before & after event


% rec_offset = 500; % duration in ms, example
% rec_endTime = rec_startTime + rec_offset; % frame index for AO recording start-time

%% initial code

lastSpikeTS = spikeClInfo.SpikeTSindex(end);
numSpikes = spikeClInfo.SpikeTSindex(end) - spikeClInfo.SpikeTSindex(1);
num_seconds = numSpikes/AO_fs;

% waveforms of spikes only necessary for figures
waves = spikeClInfo.waveForms.waves;


