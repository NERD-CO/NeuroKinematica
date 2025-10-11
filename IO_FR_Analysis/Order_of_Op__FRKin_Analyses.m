%% Order of Operations - FR-Kinematic Correlation Script

clear; close all; clc;


%% Functions in FR Analysis, and FR-Kinematic Correlation Pipeline

% compute_FRperMove_perSTNdepth(CaseDate, ephysTbl_Dir, All_SpikesPerMove_Tbl)
% [FR_SummaryTbl] = run_IO_FR_Analysis_and_Plotting(CaseDate, CaseDate_hem, ephysTbl_Dir, ephys_offset, FR_Kin_Dir);
% [kinTbl, kinSummaryTbl] = run_MovementFeatureAnalysis_IO_v2(IO_DataDir, MoveDataDir, MoveDir_CaseID);
% merge_FRKin_SummaryTbls(IO_DataDir, ephysTbl_Dir, ephys_offset, MoveDir_CaseID, FR_SummaryTbl, kinSummaryTbl)
% aggregate_FRKinematic_Correlations(KinematicsDir, Ephys_Kin_Dir, FR_Kin_Dir)


%% Directory Setup

curPCname = getenv('COMPUTERNAME');
switch curPCname
    case 'DESKTOP-I5CPDO7'  % PC_1
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'DSKTP-JTLAB-EMR'  % Lab Desktop
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\DLC_Processing'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\Kinematic Analyses'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_LFP_Analysis'
    case 'NSG-M-H8J3X34'    % PC_2
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        addpath 'C:\GitHub\NeuroKinematica\DLC_Processing'
        addpath 'C:\GitHub\NeuroKinematica\IO_FR_Analysis'
        addpath 'C:\GitHub\NeuroKinematica\Kinematic Analyses'
        addpath 'C:\GitHub\NeuroKinematica\IO_LFP_Analysis'
end


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

CaseDate = '07_13_2023_bilateral'; % Adjust as needed

% '03_23_2023'; % NER 2025
% '04_05_2023'; % NER 2025
% '05_18_2023_b_bilateral'; % NER 2025      %% errors using new pipeline
% '05_31_2023'; % INS 2026
% '06_08_2023_bilateral'; % NER 2025        %% errors using new pipeline
% '07_06_2023_bilateral'; % INS 2026
% '07_13_2023_bilateral'; % INS 2026
% '08_23_2023'; % NANS 2026

MoveDir_CaseID = 'IO_07_13_2023_LSTN'; % Adjust as needed

% 'IO_03_23_2023_LSTN'; % NER 2025
% 'IO_04_05_2023_RSTN'; % NER 2025
% 'IO_05_18_2023_b_LSTN'; % NER 2025
% 'IO_05_31_2023_LSTN';
% 'IO_06_08_2023_LSTN'; % NER 2025
% 'IO_07_06_2023_LSTN';
% 'IO_07_13_2023_LSTN';
% 'IO_08_23_2023_RSTN'; % NANS 2026


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


%% OPTIONAL: Case-specific cleaning (remove duplicates if needed)

% for 3_23_2023 case - remove duplicates / only plot for primary electrode
if strcmp(char(CaseDate), '03_23_2023')
    % All_SpikesPerMove_Tbl = All_SpikesPerMove_Tbl(158:end,1:13); % Comment or adjust as needed
    All_SpikesPerMove_Tbl = All_SpikesPerMove_Tbl(168:end,1:13);
elseif strcmp(char(CaseDate), '04_25_2023')
    All_SpikesPerMove_Tbl = [All_SpikesPerMove_Tbl(1:68,1:11); All_SpikesPerMove_Tbl(133:197,1:11); All_SpikesPerMove_Tbl(266:329,1:11)];
end

%% ===== Functions =====

%% Run MaxSpkDuration_Raster_PSTH

% outputs max spige segment duration and a raster + PSTH for all trials
[Max_SpikeDuration_samples, spikesMatrix] = MaxSpkDuration_Raster_PSTH(CaseDate, All_SpikesPerMove_Tbl, AO_spike_fs, Case_FRKin_Dir);

Max_SpkDur_seconds = Max_SpikeDuration_samples/AO_spike_fs; % seconds
Max_SpkDus_ms = Max_SpkDur_seconds * 1000; % milliseconds

%% plot_Raster_PSTH per move rep

plot_Raster_PSTH(CaseDate, All_SpikesPerMove_Tbl, AO_spike_fs, Case_FRKin_Dir)

%% Run compute_FRperMove_perSTNdepth_v3

% Outputs FR per moveRep and a summary per moveType
[FR_perTrialRep_All, FR_perMoveType_perDepth_Summary] = compute_FRperMove_perSTNdepth_v3(CaseDate, All_SpikesPerMove_Tbl, AO_spike_fs, Case_FRKin_Dir);


%% Run run_IO_FR_Analysis_and_Plotting

% Outputs FR summary per move trial per depth
[FR_SummaryTbl] = run_IO_FR_Analysis_and_Plotting(CaseDate, CaseDate_hem, All_SpikesPerMove_Tbl, Case_FRKin_Dir, useOffset, FR_Kin_Dir);

%%% update this %%%

% Example:
% IO_plotting_vCombined('03_23_2023', '', 1);                 % unilateral
% IO_plotting_vCombined('05_18_2023_b_bilateral', 'LSTN', 1); % bilateral

%% Run run_MovementFeatureAnalysis_IO_v2

fprintf('[INFO] Loading movement data from: %s\n', MoveDataDir);

% Outputs kinematics summary per move trial per depth
[kinTbl, kinSummaryTbl] = run_MovementFeatureAnalysis_IO_v2(IO_DataDir, MoveDataDir, MoveDir_CaseID);


%% Run merge_FRKin_SummaryTbls

merge_FRKin_SummaryTbls(IO_DataDir, Case_FRKin_Dir, useOffset, MoveDir_CaseID, FR_SummaryTbl, kinSummaryTbl)


%% Run aggregate_FRKinematic_Correlations

aggregate_FRKinematic_Correlations(KinematicsDir, Ephys_Kin_Dir, Case_FRKin_Dir)


%%