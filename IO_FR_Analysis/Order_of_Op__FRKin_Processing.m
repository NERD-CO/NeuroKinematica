%% Order_of_Op__FRKin_Processing

clear; close all; clc;

%% functions:

% > align_SpikesPerMove_TTL

%% Environment / Directory Set-up

% specify directory where case-specific data files are located
curPCname = getenv('COMPUTERNAME');

switch curPCname
    case 'DESKTOP-I5CPDO7'  % PC_1
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'DSKTP-JTLAB-EMR'  % Lab Desktop
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        GitHub_RepoDir = 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis';
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis'
    case 'NSG-M-H8J3X34'    % PC_2
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        GitHub_RepoDir = 'C:\GitHub\NeuroKinematica\IO_FR_Analysis';
        addpath 'C:\GitHub\NeuroKinematica\IO_FR_Analysis'
end

cd(IO_DataDir)
Subject_AO = readtable('Subject_AO.xlsx');

%% Config - Define datastream sampling rates

TTL_fs = 44000;         % Hz, Alpha Omega TTL clock
AO_spike_fs = 44000;    % Hz, Alpha Omega spike sampling rate
AO_LFP_fs = 1375;       % Hz, Alpha Omega LFP sampling rate
DLC_fs = 100;           % fps, Video/DLC frame rate

%% define offset duration

pre_offset_ms = 50; % milliseconds
offset_seconds = pre_offset_ms / 1000; % seconds

% Calculate the number of TTL samples
offset_TTLs = round(TTL_fs * offset_seconds); % ensure value is integer

useOffset = true;
% If useOffset == true or pre_offset_ms>0, useOffset_spikes function returns the #
% of samples to pre-pad in spike sample domain based on a set pre-trial offset time
% (and a meta struct for reference).
% If useOffset == false or pre_offset_ms<=0, useOffset_spike function returns 0.

% for future function input: useOffset; when = 1,  offset = offset; when = 0, offset = 0

%% Run useOffset_spikes function

[offset_spike_samples, meta_Offset_spk] = useOffset_spikes(TTL_fs, AO_spike_fs, pre_offset_ms, useOffset);


%% Inputs

% Hardcode Case-specific Data directories

CaseDate = '08_23_2023';

% '03_23_2023';             % NER 2025, NANS 2026, INS 2026      % 1
% '04_05_2023';             % NER 2025, NANS 2026, INS 2026      % 2
% '05_18_2023_b_bilateral'; % NER 2025, NANS 2026, INS 2026      % 3 (LSTN)
% '05_31_2023';                                  % INS 2026      % 4
% '06_08_2023_bilateral';   % NER 2025, NANS 2026, INS 2026      % 5 (LSTN)
% '07_06_2023_bilateral';                        % INS 2026      % 6 (LSTN)
% '07_13_2023_bilateral';                        % INS 2026      % 7 (LSTN)
% '08_23_2023';                       % NANS 2026, INS 2026      % 8


%% Notes

% '03_09_2023'; % studyID = 1, ptID 1

% '03_23_2023'; % studyID = 2, ptID 2    * % Use for INS 2026
% '04_05_2023'; % studyID = 3, ptID 2    * % Use for INS 2026

% '04_13_2023_bilateral';
% LSTN: studyID = 4(L), ptID 3
% RSTN: studyID = 5(R), ptID 3

% '05_11_2023'; % studyID = 6, ptID 4

% '05_18_2023_a'; % studyID = 7, ptID 4

% '05_18_2023_b_bilateral';
% LSTN: studyID = 8, ptID = 5            * % Use for INS 2026
% RSTN: studyID = 9, ptID = 5

% '05_31_2023'; % studyID = 10, ptID 6   * % Use for INS 2026

% '06_08_2023_bilateral';
% LSTN: studyID = 11, ptID = 7           * % Use for INS 2026
% RSTN: studyID = 12, ptID = 7

% '07_06_2023_bilateral';
% LSTN: studyID = 14, ptID = 8           * % Use for INS 2026

% '07_13_2023_bilateral';
% LSTN: studyID = 15, ptID = 9           * % Use for INS 2026
% RSTN: studyID = 16, ptID = 9

% '08_23_2023';                          * % Use for INS 2026


%% Define case-specific Kinematic data directory

MoveDataDir = [IO_DataDir, filesep, 'Processed DLC'];
cd(MoveDataDir)

% Specify case ID
MoveDir_CaseID = 'IO_08_23_2023_RSTN'; % Adjust as needed

% 'IO_03_23_2023_LSTN'; % NER 2025, NANS 2026
% 'IO_04_05_2023_RSTN'; % NER 2025, NANS 2026
% 'IO_05_18_2023_b_LSTN'; % NER 2025, NANS 2026         %% errors using new pipeline
% 'IO_05_31_2023_LSTN';
% 'IO_06_08_2023_LSTN'; % NER 2025, NANS 2026
% 'IO_07_06_2023_LSTN';
% 'IO_07_13_2023_LSTN';
% 'IO_08_23_2023_RSTN'; % NANS 2026



%% Notes

% 'IO_03_09_2023_RSTN'; % studyID = 1, ptID 1 (processed, incomplete case)

% 'IO_03_23_2023_LSTN'; % studyID = 2, ptID 2 (processed, complete case) *
% 'IO_04_05_2023_RSTN'; % studyID = 3, ptID 2 (processed, complete case) *

% 'IO_04_13_2023_LSTN'; % studyID = 4, ptID 3 (processed, complete case)

% 'IO_05_18_2023_a_RSTN'; % studyID = 7, ptID 4

% 'IO_05_18_2023_b_LSTN'; % studyID = 8, ptID 5 (processed, complete case)*
% 'IO_05_18_2023_b_RSTN'; % studyID = 9, ptID 5

% 'IO_05_31_2023_LSTN'; % studyID = 10, ptID 6 (processed, complete case) *

% 'IO_06_08_2023_LSTN'; % studyID = 11, ptID 7 (processed, complete case) *

% 'IO_07_13_2023_LSTN'; % studyID = 15, ptID 9 (processed, complete case) *
% 'IO_07_13_2023_RSTN'; % studyID = 16, ptID 9 

% 'IO_08_23_2023_RSTN';                        (processed, complete case)*


%% Define case-specific data input and outputs directories

Ephys_CaseDir = [IO_DataDir, filesep, CaseDate];       % case-specific ephys data input directory 
Move_CaseDir = [MoveDataDir, filesep, MoveDir_CaseID]; % case-specific kinematic data directory

% Define directories where case-specific IO ephys data are located (inputs)
RawDataDir = [Ephys_CaseDir, filesep, 'Raw Electrophysiology MATLAB'];       % directory where raw ephys data mat files are located (case-specific)
ProcDataDir = [Ephys_CaseDir, filesep, 'Processed Electrophysiology'];       % directory where processed ephys data should be saved (case-specific)

% save All_SpikesPerMove_Tbl in this directory (outputs)
SpikesPerMove_Dir = [Ephys_CaseDir, filesep, 'DLC_Ephys'];
if ~exist(SpikesPerMove_Dir,'dir')
    mkdir(SpikesPerMove_Dir);
end


%% Define new output directory:

% Ephys_Kinematics
FR_Kin_Dir = [IO_DataDir, filesep, 'Ephys_Kinematics', filesep, 'FR_Kinematic_Analyses'];
Case_FRKin_Dir = [FR_Kin_Dir, filesep, CaseDate];
if ~exist(Case_FRKin_Dir,'dir')
    mkdir(Case_FRKin_Dir);
end


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
    SpikesPerMove_Dir = fullfile(SpikesPerMove_Dir, CaseDate_hem);
    Case_FRKin_Dir = fullfile(Case_FRKin_Dir, CaseDate_hem);
    if ~exist(Case_FRKin_Dir,'dir')
        mkdir(Case_FRKin_Dir);
    end
    fprintf('[INFO] Hemisphere-specific output directory set: %s\n', Case_FRKin_Dir);

else
    fprintf('[INFO] Using base input directory: %s\n', ProcDataDir);
    fprintf('[INFO] Using base output directory: %s\n', Case_FRKin_Dir);
end


%% Define Clustered Spike Times directory

ClustSpkTimesDir = fullfile(ProcDataDir, 'ClusteredSpikeTimes'); % directory where clustered spike times should be saved (case-specific)
if ~isfolder(ClustSpkTimesDir)
    error('[ERROR] ClustSpkTimesDir does not exist: %s', ClustSpkTimesDir);
end
cd(ClustSpkTimesDir);


%% Run align_SpikesPerMove_TTL function

All_SpikesPerMove_Tbl = align_SpikesPerMove_TTL(Subject_AO, AO_spike_fs, TTL_fs, ProcDataDir, ClustSpkTimesDir, Move_CaseDir, pre_offset_ms, useOffset, Case_FRKin_Dir);

%%
