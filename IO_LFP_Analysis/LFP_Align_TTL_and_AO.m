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


%% Case-specific Ephys Dir Input

% isolate a specific CaseDate / studyID (StudyNum in Subject_AO csv)
CaseDate = '03_23_2023';

% '03_09_2023'; % studyID = 1, ptID 1

% '03_23_2023'; % studyID = 2, ptID 2    * % Use for INS 2026
% '04_05_2023'; % studyID = 3, ptID 2    * % Use for INS 2026

% '04_13_2023_bilateral'; ptID 3
% studyID = 4(L), 5(R),

% '05_11_2023'; % studyID = 6, ptID 4
% '05_18_2023_a'; % studyID = 7, ptID 4

% '05_18_2023_b_bilateral';
% LSTN: studyID = 8, ptID = 5    % Use for INS 2026
% RSTN: studyID = 9, ptID = 5

% '05_31_2023';  % studyID = 10, ptID 6

% '06_08_2023_bilateral'; ptID = 7
% LSTN: studyID = 11,
% RSTN: studyID = 12(R),

% '07_13_2023_bilateral';
% studyID = 15(L), 16(R), ptID = 9


%% Case-specific Movement Dir Input

% Specify case ID in Processed Movement Data Dir
Move_CaseID = 'IO_03_23_2023_LSTN';

% 'IO_03_09_2023_RSTN'; % studyID = 1, ptID 1 (processed, incomplete case)

% 'IO_03_23_2023_LSTN'; % studyID = 2, ptID 2 (processed, complete case) *
% 'IO_04_05_2023_RSTN'; % studyID = 3, ptID 2 (processed, complete case) *

% 'IO_04_13_2023_LSTN'; % studyID = 4, ptID 3 (processed, complete case)
% 'IO_04_13_2023_RSTN'; % studyID = 5, ptID 3

% 'IO_05_11_2023_LSTN'; % studyID = 6, ptID 4 (processed, incomplete case)
% 'IO_05_18_2023_a_RSTN'; % studyID = 7, ptID 4

% 'IO_05_18_2023_b_LSTN'; % studyID = 8, ptID 5 (processed, complete case) *
% 'IO_05_18_2023_b_RSTN'; % studyID = 9, ptID 5

% 'IO_05_31_2023_LSTN'; % studyID = 10, ptID 6

% 'IO_06_08_2023_LSTN'; % studyID = 11, ptID = 7 (processed, complete case)
% 'IO_06_08_2023_RSTN'; % studyID = 12, ptID = 7 (processed, incomplete case)

% 'IO_07_13_2023_LSTN'; % studyID = 15, ptID = 9
% 'IO_07_13_2023_RSTN'; % studyID = 16, ptID = 9


%% Define Input and Output data directories

Case_DataDir = [IO_DataDir, filesep, CaseDate];
MoveDataDir = [IO_DataDir, filesep, 'Processed DLC'];

% Case-specific Input dirs
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology'];       % directory for processed ephys data and spike clusters
Move_CaseDir = [MoveDataDir, filesep, Move_CaseID];                         % directory for processed DLC data and Movement Indices

% Case-specific Output dir
ephysTbl_Dir = [Case_DataDir, filesep, 'DLC_Ephys'];                        % directory where all ephys per move-rep tables are located


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
    fprintf('[INFO] Hemisphere-specific input ephys directory set: %s\n', ProcDataDir);
    ephysTbl_Dir = fullfile(ephysTbl_Dir, CaseDate_hem);
    fprintf('[INFO] Hemisphere-specific output directory set: %s\n', ephysTbl_Dir);

else
    fprintf('[INFO] Using base ephys input directory: %s\n', ProcDataDir);
    fprintf('[INFO] Using base output directory: %s\n', ephysTbl_Dir);
end

%% Define datastream sampling rates (Alpha Omega and Video sampling fs)

TTL_fs = 44000; % Hz
AO_spike_fs = 44000; % Hz
AO_LFP_fs = 1375; % Hz
DLC_fs = 100; % fps


%% Define offset duration

offset_ms = 50; % milliseconds
offset_seconds = offset_ms / 1000; % seconds

% Calculate number of TTL samples
offset_TTLs = round(TTL_fs * offset_seconds); % ensure value is integer

% Calculate number TTL samples in AO_LFP sample domain by downsampling TTL_fs
offset_TTLs_LFP = round((offset_TTLs/TTL_fs)*AO_LFP_fs); % ensure value is integer

% for future function input: useOffset; when = 1,  offset = offset; when = 0, offset = 0


%% Define case-specific directory for movement indices per trial

cd(Move_CaseDir)

% Move_CaseDir data subfolders:
Move_CaseMats = [Move_CaseDir, filesep, 'mat folder'];      % contains processed DLC timeseries data (csv-to-mat)
Move_CaseVideos = [Move_CaseDir, filesep, 'video folder'];  % contains DLC-labeled videos and Movement Index CSVs
cd(Move_CaseVideos)

% % for kinematic analyses
% cd(Move_CaseMats)
% moveMat = dir('*.mat');
% moveMat_names = {moveMat.name};



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


%% save All_SpikesPerMove_Tbl to a file

cd(ephysTbl_Dir)
save('All_LFPsPerMove_offset.mat',"All_LFPsPerMove_Tbl");


