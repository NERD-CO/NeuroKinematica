%% IO LFP 

clear; clc;
addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_LFP_Analysis'

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


%% Inputs and Outputs

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


%% Define case-specific Ephys data directories

Case_DataDir = [IO_DataDir, filesep, CaseDate];

% Define directories where case-specific IO ephys data are located (inputs)
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology'];       % directory where processed ephys .matdata should be saved (inputs)
ephysTbl_Dir = [Case_DataDir, filesep, 'DLC_Ephys'];                        % directory where all ephys per move-rep tables are located (outputs)


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
    fprintf('[INFO] Using base input directory: %s\n', ProcDataDir);
    fprintf('[INFO] Using base output directory: %s\n', ephysTbl_Dir);
end


%% Define case-specific Kinematic data directory

MoveDataDir = [IO_DataDir, filesep, 'Processed DLC'];
cd(MoveDataDir)

% Specify case ID
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

%%  Isolate case-specific kinematic data directory

Move_CaseDir = [MoveDataDir, filesep, Move_CaseID];

%% Run align_LFPsPerMove_TTL function

% align_LFPsPerMove_TTL(Subject_AO, ProcDataDir, ClustSpkTimesDir, Move_CaseDir, ephysTbl_Dir)



%% Define case-specific directory for movement indices per trial

% Move_CaseDir data subfolders:
Move_CaseMats = [Move_CaseDir, filesep, 'mat folder'];      % contains processed DLC timeseries data (csv-to-mat)
Move_CaseVideos = [Move_CaseDir, filesep, 'video folder'];  % contains DLC-labeled videos and Movement Index CSVs
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

% Calculate number of TTL samples
offset_TTLs = round(TTL_fs * offset_seconds); % ensure value is integer

% Calculate number TTL samples in AO_LFP sample domain by downsampling TTL_fs
offset_TTLs_LFP = round((offset_TTLs/TTL_fs)*AO_LFP_fs); % ensure value is integer

% for future function input: useOffset; when = 1, use offset; when = 0, don't

%% Load all processed LFPs %%% fix this block later

cd(ProcDataDir)

% list of filename
LFPmatfiles = dir('*.mat');
LFPmatnames = {LFPmatfiles.name};

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

        % Find row of ao_MAT_file that corresponds with trial
        SubjectAO_row = Subject_AO(contains(Subject_AO.ao_MAT_file,ProcName),:);

        switch SubjectAO_row.stn_loc{1}(1)
            case 'd'
                depthName = 't'; % dorsal
            case 'v'
                depthName = 'b'; % ventral
            otherwise
                depthName = 'c'; % central
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

        electrode_LFP_name = [ProcName,'_', electrod_names{e_names}]; %temp var
        LFP_raw_depthName = [ProcName,'_', motor_trial_ID]; % temp var

    end

end


%% Load all LFPs per Movement Tbl %%% refine this block 

cd(ephysTbl_Dir)

% list of filename
Tblmatfiles = dir('*.mat');
Tblmatnames = {Tblmatfiles.name};

use_offset = 1;
for Tblmat_name = 1:length(Tblmatnames)
    fileparts = split(Tblmatnames{Tblmat_name},'_');
    if use_offset == 1
        Tbl_name = fileparts{2:3};
    else
        Tbl_name = fileparts{2};
    end
end


LFPsPerMove_Tbl = Tblmatnames{contains(Tblmatnames, 'All_LFPsPerMove_offset')};
load(LFPsPerMove_Tbl, 'All_LFPsPerMove_Tbl')

%% Quastions:

% MLFP vs LFP - does LFP = CLFP?

%% spectral interp

% loop on harmonics of 60 Hz

%% filters

% notch
% low-pass at 250 Hz
% high-pass at 0.5 or 1 Hz
% 1-3 kHz
% downsample to 500 Hz

%% PSDs on different context

% rest, H O/C, A Pro/Sup, E Flex/Exten

%% FOOOFs on different contexts

%% Hilbert transformation spec. freq. (inst. power)

% theta band
% low-beta band
% high-beta band
% gamma band
% HFO

% at specific movement timepoints

%% Spike-field coherence

% inst. beta power at specific time/context (move index)

%% LFP bursting analysis / continuious wavelet transformation




