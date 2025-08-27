%% Order of Operations - FR-Kinematic Correlation Script

clear; clc;
addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\DLC_Processing'
addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis'
addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\Kinematic Analyses'


%% Functions in FR Analysis, and FR-Kinematic Correlation Pipeline

% compute_FRperMove_perSTNdepth(CaseDate, ephysTbl_Dir, All_SpikesPerMove_Tbl)
% [FR_SummaryTbl] = run_IO_FR_Analysis_and_Plotting(CaseDate, CaseDate_hem, ephysTbl_Dir, ephys_offset, FR_Kin_Dir);
% [kinTbl, kinSummaryTbl] = run_MovementFeatureAnalysis_IO_v2(IO_DataDir, MoveDataDir, MoveDir_CaseID);
% merge_FRKin_SummaryTbls(IO_DataDir, ephysTbl_Dir, ephys_offset, MoveDir_CaseID, FR_SummaryTbl, kinSummaryTbl)
% aggregate_FRKinematic_Correlations(KinematicsDir, Ephys_Kin_Dir, FR_Kin_Dir)


%% Directory Setup

curPCname = getenv('COMPUTERNAME');
switch curPCname
    case 'DESKTOP-I5CPDO7'
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    otherwise
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
end


%% Inputs:

CaseDate = '05_18_2023_b_bilateral'; % Adjust as needed
% '03_23_2023';
% '04_05_2023';
% '05_18_2023_b_bilateral';
% '06_08_2023_bilateral';

MoveDir_CaseID = 'IO_05_18_2023_b_LSTN'; % Adjust as needed
% 'IO_03_23_2023_LSTN';
% 'IO_04_05_2023_RSTN';
% 'IO_05_18_2023_b_LSTN';
% 'IO_06_08_2023_LSTN'


ephys_offset = 1;

% Data folder paths
Case_DataDir = fullfile(IO_DataDir, CaseDate);
ephysTbl_Dir = fullfile(Case_DataDir,'DLC_Ephys');
MoveDataDir = fullfile(IO_DataDir, 'Processed DLC');
KinematicsDir = fullfile(IO_DataDir, 'Kinematic Analyses');
Ephys_Kin_Dir = fullfile(IO_DataDir, 'Ephys_Kinematics');
FR_Kin_Dir = fullfile(Ephys_Kin_Dir, 'FR_Kinematic_Analyses');
Case_FR_Kin = fullfile(FR_Kin_Dir, CaseDate); % case-specific results from run_FR_KinematicCorr saved here

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
    ephysTbl_Dir = fullfile(ephysTbl_Dir, CaseDate_hem); 
else
    CaseDate_hem = '';
end
fprintf('[INFO] Loading spike data from: %s\n', ephysTbl_Dir);


%% Load All_SpikesPerMove_Tbl

cd(ephysTbl_Dir)
Tbl_list = dir('*Spikes*.mat');
Tbl_names = {Tbl_list.name};
spk_case = Tbl_names{contains(Tbl_names, 'offset')}; % offset version preferred
load(spk_case, 'All_SpikesPerMove_Tbl');

%% OPTIONAL: Case-specific cleaning (remove duplicates if needed)

% for 3_23_2023 case - remove duplicates / only plot for primary electrode
if CaseDate == '03_23_2023'
    All_SpikesPerMove_Tbl = All_SpikesPerMove_Tbl(158:end,1:13); % Comment or adjust as needed
else
    All_SpikesPerMove_Tbl = All_SpikesPerMove_Tbl;
end


%% ===== Functions =====

%% Run compute_FRperMove_perSTNdepth

compute_FRperMove_perSTNdepth(CaseDate, ephysTbl_Dir, All_SpikesPerMove_Tbl)

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

