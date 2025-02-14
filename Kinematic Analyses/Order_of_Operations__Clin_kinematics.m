%% Order of Operations - Clinical Kinematic Analyses

clear all; close all; clc;

%% I. Processing

addpath C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\DLC_Processing

% 1) Run CSV-to-MAT processing function: 
% run 'run_dlc_processCSV2.m'

% 2) Generate movement indices 
% C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\DLC_VideoCheck_GUI
% run 'ER_DLC_MoceCheck_Dual_v4_mlapp

%% II. Analysis 

addpath 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\Kinematic Analyses'

%% 1) Define function inputs

% Define mainDir

mainDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\Kinematic Analyses';

% Define caseID (date or subjectID) and hemisphere

% case_hem = '09_12_2023_LSTN';
% case_hem = '09_12_2023_RSTN';

case_hem = 'MDTII_MDT5_LSTN';
% case_hem = 'MDTII_MDT5_RSTN';


%% 2) Process and visualize movement timeseries data

run_MovementProcessing_Clin_v2(mainDir, case_hem)

%% 3) by movement type

%run_MovementProcessing_Clin_v3(mainDir, casedate_hem)


%% Optional: Process and visualize movement timeseries data cleaned by artifact rejection function

run_MovementProcessing_Clin_artifactRejection(mainDir, case_hem)


%% 3) stats

% Define switch case inputs (for run_MovementStats_Clin_v2 function)

hemisphere = 'L';
% hemisphere = 'R';

switch hemisphere
    case 'L'

        mainDir2 = [mainDir , filesep , 'MDTII_MDT5_LSTN'];

    case 'R'

        mainDir2 = [mainDir , filesep , 'MDTII_MDT5_RSTN'];
end


run_MovementStats_Clin_v2(mainDir, case_hem, hemisphere)

% fix lines 87-91 in ttest_plot_2States_RandTrim s.t. output files save
% within specified outDataDir folder (movementStats)





