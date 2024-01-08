%% Order of Operations - Clinical Kinematic Analyses

% cd 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\Kinematic Analyses'
addpath 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\Kinematic Analyses'

%% 1) Define function inputs

% Define mainDir

mainDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\Kinematic Analyses';

% Define casedate and hemisphere

casedate_hem = '09_12_2023_LSTN';
% casedate_hem = '09_12_2023_RSTN';


% Define switch case inputs (for run_MovementStats_Clin_v2 function)

hemisphere = 'L';

switch hemisphere
    case 'L'

        mainDir2 = [mainDir , filesep , '09_12_2023_LSTN'];

    case 'R'

        mainDir2 = [mainDir , filesep , '09_12_2023_RSTN'];
end


%% 2) Process and visualize movement timeseries data

run_MovementProcessing_Clin_v2(mainDir, casedate_hem)


%% Optional: Process and visualize movement timeseries data cleaned by artifact rejection function

% run_MovementProcessing_Clin_artifactRejection(mainDir, casedate_hem)


%% 3) stats

run_MovementStats_Clin_v2(mainDir, casedate_hem, hemisphere)

% fix lines 87-91 in ttest_plot_2States_RandTrim s.t. output files save
% within specified outDataDir folder (movementStats)





