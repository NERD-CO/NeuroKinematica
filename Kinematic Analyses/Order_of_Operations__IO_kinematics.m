%% Order of Operations - IO Kinematic Analyses

% cd 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\Kinematic Analyses'
addpath 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\Kinematic Analyses'

%% 1) Define function inputs

% Define mainDir

mainDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\Kinematic Analyses';

% Define casedate and hemisphere

% casedate_hem = '03_09_2023_RSTN';     % does not work
% casedate_hem = 'IO_05_11_2023_LSTN';
 casedate_hem = 'IO_06_08_2023_LSTN';
% casedate_hem = 'IO_06_08_2023_RSTN';


%% 2a) Process and visualize movement timeseries data - v1

% run_MovementProcessing_IO_v1(mainDir, casedate_hem)

%% 2b) Process and visualize movement timeseries data - v2

run_MovementProcessing_IO_v2(mainDir, casedate_hem)

%% 2c) Process and visualize movement timeseries data - v2

% run_MovementProcessing_IO_jatv1(mainDir, casedate_hem)


%% Optional: Process and visualize movement timeseries data cleaned by artifact rejection function

% run_MovementProcessing_IO_artifactRejection(mainDir, casedate_hem)


%% 3) stats

% run_MovementStats_IO_v1(mainDir, casedate_hem)






