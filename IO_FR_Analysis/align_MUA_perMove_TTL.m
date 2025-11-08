

% Align IO spike data segments with corresponding movement data segments per trial


%% define offset duration

offset_seconds = pre_offset_ms / 1000; % seconds

% Calculate the number of TTL samples
offset_TTLs = round(TTL_fs * offset_seconds); % ensure value is integer

% for future function input: useOffset; when = 1,  offset = offset; when = 0, offset = 0

%% Run useOffset_spikes function

[offset_spike_samples, meta_Offset_spk] = useOffset_spikes(TTL_fs, AO_spike_fs, pre_offset_ms, useOffset);

%% Directory Setup

curPCname = getenv('COMPUTERNAME');
switch curPCname
    case 'DESKTOP-I5CPDO7'  % PC_1
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'DSKTP-JTLAB-EMR'  % Lab Desktop
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\DLC_Processing'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\zetatest'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\Kinematic Analyses'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_LFP_Analysis'
    case 'NSG-M-H8J3X34'    % PC_2
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        addpath 'C:\GitHub\NeuroKinematica\DLC_Processing'
        addpath 'C:\GitHub\NeuroKinematica\IO_FR_Analysis'
        addpath 'C:\GitHub\NeuroKinematica\zetatest'
        addpath 'C:\GitHub\NeuroKinematica\Kinematic Analyses'
        addpath 'C:\GitHub\NeuroKinematica\IO_LFP_Analysis'
end
cd(IO_DataDir)
Subject_Hem_CaseMap = readtable('Subject_Hem_MetaSummary.xlsx'); % Subject_Hem_MetaSummary.xlsx


%% Config - Define datastream sampling rates (Alpha Omega and Video sampling fs)

TTL_fs = 44000;         % Hz, Alpha Omega TTL clock
AO_spike_fs = 44000;    % Hz, Alpha Omega spike sampling rate
AO_LFP_fs = 1375;       % Hz, Alpha Omega LFP sampling rate
DLC_fs = 100;           % fps, Video/DLC frame rate


%% Config - Define pre-trial offset duration for useOffset_spikes function

pre_offset_ms = 50; % milliseconds
offset_seconds = pre_offset_ms / 1000; % seconds

% Calculate number of TTL samples
offset_TTLs = round(TTL_fs * offset_seconds); % ensure value is integer

% Calculate number TTL samples in AO_LFP sample domain by downsampling TTL_fs
offset_spikes = round((offset_TTLs/TTL_fs)*AO_spike_fs); % ensure value is integer

useOffset = true;
% If useOffset == true or pre_offset_ms>0, useOffset_spikes function returns the #
% of samples to pre-pad in spike sample domain based on a set pre-trial offset time
% (and a meta struct for reference).
% If useOffset == false or pre_offset_ms<=0, useOffset_spike function returns 0.


%% Inputs:

% Ephys data folder:
CaseDate = '04_05_2023'; % Adjust as needed

% 1: '%% Directory Setup

curPCname = getenv('COMPUTERNAME');
switch curPCname
    case 'DESKTOP-I5CPDO7'  % PC_1
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'DSKTP-JTLAB-EMR'  % Lab Desktop
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\DLC_Processing'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\zetatest'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\Kinematic Analyses'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_LFP_Analysis'
    case 'NSG-M-H8J3X34'    % PC_2
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        addpath 'C:\GitHub\NeuroKinematica\DLC_Processing'
        addpath 'C:\GitHub\NeuroKinematica\IO_FR_Analysis'
        addpath 'C:\GitHub\NeuroKinematica\zetatest'
        addpath 'C:\GitHub\NeuroKinematica\Kinematic Analyses'
        addpath 'C:\GitHub\NeuroKinematica\IO_LFP_Analysis'
end
cd(IO_DataDir)
Subject_Hem_CaseMap = readtable('Subject_Hem_MetaSummary.xlsx'); % Subject_Hem_MetaSummary.xlsx


%% Config - Define datastream sampling rates (Alpha Omega and Video sampling fs)

TTL_fs = 44000;         % Hz, Alpha Omega TTL clock
AO_spike_fs = 44000;    % Hz, Alpha Omega spike sampling rate
AO_LFP_fs = 1375;       % Hz, Alpha Omega LFP sampling rate
DLC_fs = 100;           % fps, Video/DLC frame rate


%% Config - Define pre-trial offset duration for useOffset_spikes function

pre_offset_ms = 50; % milliseconds
offset_seconds = pre_offset_ms / 1000; % seconds

% Calculate number of TTL samples
offset_TTLs = round(TTL_fs * offset_seconds); % ensure value is integer

% Calculate number TTL samples in AO_LFP sample domain by downsampling TTL_fs
offset_spikes = round((offset_TTLs/TTL_fs)*AO_spike_fs); % ensure value is integer

useOffset = true;
% If useOffset == true or pre_offset_ms>0, useOffset_spikes function returns the #
% of samples to pre-pad in spike sample domain based on a set pre-trial offset time
% (and a meta struct for reference).
% If useOffset == false or pre_offset_ms<=0, useOffset_spike function returns 0.


%% Inputs:

% Ephys data folder:
CaseDate = '03_23_2023'; % Adjust as needed

% 1: '03_23_2023';             % NER 2025
% 1: '04_05_2023';             % NER 2025

% 2: '04_13_2023_bilateral';   % GRC 2026      %% errors - Unrecognized field name "dblZetaT".

% 3: '05_18_2023_b_bilateral'; % NER 2025      %% errors using new pipeline --> check MIs;
%   %   % Other error: Unrecognized field name "dblZetaT".

% 4: '05_31_2023';             % INS 2026

% 5: '06_08_2023_bilateral';   % NER 2025      %% errors using new pipeline --> check MIs;
%   %   % Other error: Error using findpeaks>parse_inputs (line 225): Data set must contain at least 3 samples.

% 6: '07_06_2023_bilateral';   % INS 2026
% 7: '07_13_2023_bilateral';   % INS 2026
% 10: '08_23_2023';             % INS 2026
% 12: '11_30_2023_bilateral';   % GRC 2026


% Kinematic data folder:
MoveDir_CaseID = 'IO_03_23_2023_LSTN'; % Adjust as needed

% 'IO_03_23_2023_LSTN';   % NER 2025
% 'IO_04_05_2023_RSTN';   % NER 2025
% 'IO_04_13_2023_LSTN';   % GRC 2026
% 'IO_05_18_2023_b_LSTN'; % NER 2025
% 'IO_05_31_2023_LSTN';   % INS 2026
% 'IO_06_08_2023_LSTN';   % NER 2025
% 'IO_07_06_2023_LSTN';   % INS 2026
% 'IO_07_13_2023_LSTN';   % INS 2026
% 'IO_07_13_2023_RSTN';   % GRC 2026
% 'IO_08_23_2023_RSTN';   % INS 2026
% 'IO_11_30_2023_LSTN';   % GRC 2026


%% Data folder paths

Case_DataDir = fullfile(IO_DataDir, CaseDate);
% ephysTbl_Dir = fullfile(Case_DataDir,'DLC_Ephys');
MoveDataDir = fullfile(IO_DataDir, 'Processed DLC');
KinematicsDir = fullfile(IO_DataDir, 'Kinematic Analyses');
Ephys_Kin_Dir = fullfile(IO_DataDir, 'Ephys_Kinematics');
FR_Kin_Dir = fullfile(Ephys_Kin_Dir, 'FR_Kinematic_Analyses');
Case_FRKin_Dir = fullfile(FR_Kin_Dir, CaseDate); % input dir (new ephysTbleDir) and output (results) dir


%% Handle Bilateral Cases

% Specify hemisphere in command window if bilateral
isBilateral = contains(CaseDate,'bilateral','IgnoreCase',true);
if isBilateral
    fprintf('[INFO] Bilateral case detected: %s\n',CaseDate);
    CaseDate_hem = input('Enter hemisphere (LSTN or RSTN): ','s');
    if ~ismember(CaseDate_hem,{'LSTN','RSTN'})
        error('Invalid hemisphere');
    end
    % Append hemisphere-specific folder
    Case_FRKin_Dir = fullfile(Case_FRKin_Dir, CaseDate_hem);
else
    CaseDate_hem = '';
end
fprintf('[INFO] Loading spike data from: %s\n', Case_FRKin_Dir);


%% Load All_SpikesPerMove_Tbl

cd(Case_FRKin_Dir)
Tbl_list = dir('*Spikes*.mat');
Tbl_names = {Tbl_list.name};
spk_case = Tbl_names{contains(Tbl_names, 'offset')}; % offset version preferred
load(spk_case, 'All_SpikesPerMove_Tbl');

%% Data folder paths

Case_DataDir = fullfile(IO_DataDir, CaseDate);
% ephysTbl_Dir = fullfile(Case_DataDir,'DLC_Ephys');
MoveDataDir = fullfile(IO_DataDir, 'Processed DLC');
KinematicsDir = fullfile(IO_DataDir, 'Kinematic Analyses');
Ephys_Kin_Dir = fullfile(IO_DataDir, 'Ephys_Kinematics');
FR_Kin_Dir = fullfile(Ephys_Kin_Dir, 'FR_Kinematic_Analyses');
Case_FRKin_Dir = fullfile(FR_Kin_Dir, CaseDate); % input dir (new ephysTbleDir) and output (results) dir


%% Handle Bilateral Cases

% Specify hemisphere in command window if bilateral
isBilateral = contains(CaseDate,'bilateral','IgnoreCase',true);
if isBilateral
    fprintf('[INFO] Bilateral case detected: %s\n',CaseDate);
    CaseDate_hem = input('Enter hemisphere (LSTN or RSTN): ','s');
    if ~ismember(CaseDate_hem,{'LSTN','RSTN'})
        error('Invalid hemisphere');
    end
    % Append hemisphere-specific folder
    Case_FRKin_Dir = fullfile(Case_FRKin_Dir, CaseDate_hem);
else
    CaseDate_hem = '';
end
fprintf('[INFO] Loading spike data from: %s\n', Case_FRKin_Dir);


%% Load All_SpikesPerMove_Tbl

cd(Case_FRKin_Dir)
Tbl_list = dir('*Spikes*.mat');
Tbl_names = {Tbl_list.name};
spk_case = Tbl_names{contains(Tbl_names, 'offset')}; % offset version preferred
load(spk_case, 'All_SpikesPerMove_Tbl');

%%
cd(Case_DataDir)
ProcEphys_dir = [Case_DataDir, filesep, 'Processed Electrophysiology'];
cd(ProcEphys_dir)

% Loop through rows in AllSpikesPerMove table
for t_i = 1:height(All_SpikesPerMove_Tbl)
    temp_mat = All_SpikesPerMove_Tbl.spike_trial_ID{t_i};
    load(temp_mat);
    allSPIKEdata = ProcEphys.Spike;

    % loop through electrodes (E1 ... En)
    num_E = numel(fieldnames(allSPIKEdata));
    for e_i = 1:num_E

        spikeDataRaw = allSPIKEdata.(['E', num2str(e_i)]).rawData;
        spkFS = allSPIKEdata.(['E', num2str(e_i)]).Hz*1000;

        sampLen = round(spkFS/1000);

        thresh = (round(std(double(spikeDataRaw)) * 4) + mean(double(spikeDataRaw)));
        noise = (round(std(double(spikeDataRaw)) * 10) + mean(double(spikeDataRaw)));
        minDist = round(spkFS/1000) * 1.5; % was 1.5 ms

        % MUA (unsorted spikes that exeeded thresholds
        [ waveData ] = extractWaveforms_Clz_v07(double(spikeDataRaw), thresh, noise, minDist, spkFS, 1); %%%%%%%%%%%%%%%%%%%%% PNFLAG = 1
        waveForms = waveData;

        peakIndex = waveData.allWavesInfo.alllocs; % MUA spike times

        
        % allocs within trial_start to trail_end --> MUA

        % time math --> MUA_ts

        % just store E1 (don't worry about storing MUA and MUA_ts for more than)



    end


    % MUA

    % MUA_ts






end

