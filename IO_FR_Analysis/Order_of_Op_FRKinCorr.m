%% Order of Operations - FR-Kinematic Correlation Script

clear; clc;
addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis'
addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\Kinematic Analyses'


%% Functions in FR Analysis, and FR-Kinematic Correlation Pipeline

% > spikes_Align_TTL_and_AO
% > compute_meanFR_perMove_perSTNdepth_v2
% IO_plotting_vCombined(CaseDate, CaseDate_hem, ephys_offset)
% [kinTbl, kinSummaryTbl] = run_MovementFeatureAnalysis_IO_v2(IO_DataDir, MoveDataDir, MoveDir_CaseID);


%% Directory Setup

curPCname = getenv('COMPUTERNAME');
switch curPCname
    case 'DESKTOP-I5CPDO7'
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    otherwise
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
end

%% Inputs:

CaseDate = '04_05_2023'; % Adjust as needed
% '03_23_2023';
% '04_05_2023';
% '05_18_2023_b_bilateral';
% '06_08_2023_bilateral';

MoveDir_CaseID = 'IO_04_05_2023_RSTN'; % Adjust as needed
% 'IO_03_23_2023_LSTN';
% 'IO_04_05_2023_RSTN';
% 'IO_05_18_2023_b_LSTN';
% 'IO_06_08_2023_LSTN'

ephys_offset = 1;

Case_DataDir = fullfile(IO_DataDir, CaseDate);
ephysTbl_Dir = fullfile(Case_DataDir,'DLC_Ephys');
MoveDataDir = fullfile(IO_DataDir, 'Processed DLC');
KinematicsDir = fullfile(IO_DataDir, 'Kinematic Analyses');
Ephys_Kin_Dir = fullfile(IO_DataDir, 'Ephys_Kinematics');
FR_Kin_Dir = fullfile(Ephys_Kin_Dir, 'FR_Kinematic_Analyses');
Case_FR_Kin = fullfile(FR_Kin_Dir, CaseDate); % case-specific results from run_FR_KinematicCorr saved here


isBilateral = contains(CaseDate,'bilateral','IgnoreCase',true);
if isBilateral
    fprintf('[INFO] Bilateral case detected: %s\n',CaseDate);
    CaseDate_hem = input('Enter hemisphere (LSTN or RSTN): ','s');
    if ~ismember(CaseDate_hem,{'LSTN','RSTN'}), error('Invalid hemisphere'); end
    ephysTbl_Dir = fullfile(ephysTbl_Dir,CaseDate_hem);
else
    CaseDate_hem = '';
end
fprintf('[INFO] Loading spike data from: %s\n', ephysTbl_Dir);


%% Run run_IO_FR_Analysis_and_Plotting

[FR_SummaryTbl] = run_IO_FR_Analysis_and_Plotting(CaseDate, CaseDate_hem, ephysTbl_Dir, ephys_offset, FR_Kin_Dir);

% Example:
% IO_plotting_vCombined('03_23_2023', '', 1);                 % unilateral
% IO_plotting_vCombined('05_18_2023_b_bilateral', 'LSTN', 1); % bilateral

%% Run run_MovementFeatureAnalysis_IO_v2

fprintf('[INFO] Loading movement data from: %s\n', MoveDataDir);

[kinTbl, kinSummaryTbl] = run_MovementFeatureAnalysis_IO_v2(IO_DataDir, MoveDataDir, MoveDir_CaseID);


%% Run run_FR_KinematicCorr

merge_FRKin_SummaryTbls(IO_DataDir, ephysTbl_Dir, ephys_offset, MoveDir_CaseID, FR_SummaryTbl, kinSummaryTbl)


%% Run aggregate_FRKinematic_Correlations

aggregate_FRKinematic_Correlations(KinematicsDir, Ephys_Kin_Dir, FR_Kin_Dir)

