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
% CaseDate = '03_09_2023'; % studyID = 1
CaseDate = '03_23_2023'; % studyID = 2

% CaseDate = '06_08_2023_bilateral';
    % CaseDate_hem = 'LSTN'; % comment out when N/A

% define case-specific data directory
Case_DataDir = [IO_DataDir, filesep, CaseDate];

% directories where case-specific IO ephys data are located
RawDataDir = [Case_DataDir, filesep, 'Raw Electrophysiology MATLAB'];       % directory where raw MATLAB data files are located (case-specific)
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology'];       % directory where processed MATLAB data should be saved (case-specific)
    % ProcDataDir = [ProcDataDir, filesep, CaseDate_hem];                     % comment out when N/A

%% directory for movement indices

% define kinematic data directory
MoveDataDir = [IO_DataDir, filesep, 'Kinematic Analyses'];

% specify case ID
% Move_CaseID = 'IO_03_09_2023_RSTN'; % studyID = 1
Move_CaseID = 'IO_03_23_2023_LSTN'; % studyID = 2

% Move_CaseID ='IO_06_08_2023_LSTN';

% isolate case-specific kinematic data directory
Move_CaseDir = [MoveDataDir, filesep, Move_CaseID];

% data subfolders:
Move_CaseMats = [Move_CaseDir, filesep, 'mat folder'];
Move_CaseVideos = [Move_CaseDir, filesep, 'video folder'];

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

% Calculate number of TTL samples
offset_TTLs = round(TTL_fs * offset_seconds); % ensure value is integer

% Calculate number TTL samples in AO_LFP sample domain by downsampling TTL_fs
offset_TTLs_LFP = round((offset_TTLs/TTL_fs)*AO_LFP_fs); % ensure value is integer

% for future function input: useOffset; when = 1,  offset = offset; when = 0, offset = 0

%% go to ProcDataDir
cd(ProcDataDir)

% list of filename
LFPmatfiles = dir('*.mat');
LFPmatnames = {LFPmatfiles.name};

All_LFPsPermove = {length(LFPmatnames)};

for LFP_mat_name = 1:length(LFPmatnames)

    cd(ProcDataDir)
    load(LFPmatnames{LFP_mat_name},'ProcEphys')

    fileparts = split(LFPmatnames{LFP_mat_name},'_');
    ProcName = fileparts{2};

    ProcFile = LFPmatnames{contains(LFPmatnames, ProcName)};
    load(ProcFile, 'ProcEphys')

    % account for mult. electrode channels
    electrod_names = fieldnames(ProcEphys.LFP); % get num
    for e_names = 1:length(electrod_names)

        LFP_raw = ProcEphys.LFP.(electrod_names{e_names}).rawData; % dynamically index within a struct

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
        LFPMoveTbl = moveTbl(cellfun(@(X) ~isempty(X), moveTbl.MoveType, 'UniformOutput', true),:); % clean moveTbl (remove zeros)

        electrode_LFP_name = [ProcName,'_', electrod_names{e_names}]; %temp var

        % initialize arrays
        LFP_trial_ID = repmat({electrode_LFP_name}, height(LFPMoveTbl), 1);
        move_trial_ID = repmat({motor_trial_ID}, height(LFPMoveTbl), 1);
        TTL_LFP_idx_Start = nan(height(LFPMoveTbl), 1);

        % loop through trials and pull out BeginF and EndF indices in MoveIndex csv per corresponding trial
        for move_i = 1:height(LFPMoveTbl)

            % define frame indices based on task recording context (frame indices in TTL_Down)                    % get these from movement index csv
            frame_startTime = LFPMoveTbl.BeginF(move_i); % index of recording initiation
            frame_endTime = LFPMoveTbl.EndF(move_i); % index of recording termination

            % extract TTL signals of interest based on task context
            TTL_samp_taskStart = TTL_Down(frame_startTime); % number of samples wrt TTL clock
            TTL_samp_taskEnd = TTL_Down(frame_endTime);

            % downsample TTL_fs - convert ^ to time; multiply AO_LFP_fs
            TTL_samp_taskStart = round((TTL_samp_taskStart/TTL_fs)*AO_LFP_fs);
            TTL_samp_taskEnd = round((TTL_samp_taskEnd/TTL_fs)*AO_LFP_fs);

            % define AO recording times
            AO_startTime = ProcEphys.LFP.(electrod_names{e_names}).startTime;

            % calculate difference in clock startTimes
            time_offset = TTL_clockStart - AO_startTime; % time (seconds) wrt TTL clock

            % convert time offset to sample offset
            sample_offset = round(time_offset*AO_LFP_fs); % number of samples

            % convert TTL sample indices by the sample offset to transform into AO_LFP clock / LFP sample domain
            TTL_LFP_idx_Start(move_i) = TTL_samp_taskStart + sample_offset; % number of samples wrt AO clock
            TTL_LFP_idx_End = TTL_samp_taskEnd + sample_offset;

            % incorporate offset
            TTL_LFP_idx_Start(move_i) =  TTL_LFP_idx_Start(move_i) - offset_TTLs_LFP; % start [50 ms] before

            % determine LFPs within movement block
            LFPs_in_move1 = LFP_raw(TTL_LFP_idx_Start(move_i):TTL_LFP_idx_End); %  LFPs w/in moveblock in AO_time

            % populate LFPMoveTable row with LFPs within movement block
            LFPMoveTbl.('LFPs'){move_i} = LFPs_in_move1;

            % Add LFP_ID, move_ID, and TTL_LFP_idx_Start columns
            LFPMoveTbl.LFP_trial_ID = LFP_trial_ID;
            LFPMoveTbl.move_trial_ID = move_trial_ID;
            LFPMoveTbl.TTL_LFP_idx_Start = TTL_LFP_idx_Start;
            
        end

        All_LFPsPermove{LFP_mat_name} = LFPMoveTbl;

    end

end

% Define the standard column order
standard_col_order = {'MoveN', 'MoveType', 'BeginF', 'EndF', 'TTL_LFP_idx_Start', 'LFP_trial_ID', 'move_trial_ID', 'LFPs'};

% Initialize All_data as a cell array
All_data = {};
% Loop through each table in All_moveTbl_array
for tbl_i = 1:length(All_LFPsPermove)
    tbl_1 = All_LFPsPermove{tbl_i};
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
All_LFPsPerMove_Tbl = cell2table(All_data, 'VariableNames', standard_col_order);


% save All_SpikesPerMove_Tbl to a file
LFPsPerMove_Dir = [Case_DataDir, filesep, 'DLC_Ephys'];
    % LFPsPerMove_Dir = [Case_DataDir, filesep, 'DLC_Ephys', filesep, CaseDate_hem]; % comment out when N/A
cd(LFPsPerMove_Dir)
%save('All_LFPsPerMove_offset.mat',"All_LFPsPerMove_Tbl");


