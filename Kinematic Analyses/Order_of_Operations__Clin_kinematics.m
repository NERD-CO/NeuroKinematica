%% Order of Operations - Clinical Kinematic Analyses

% cd 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\Kinematic Analyses'
addpath 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\Kinematic Analyses'

%% 1) Define function inputs

casedate = '09_12_2023';
hemisphere = 'L'; 

% Completed cases:
% casedate = '09_12_2023';
    % hemisphere = 'L';
    % hemisphere = 'R'

    

%% 2) Process and visualize movement timeseries data

run_MovementProcessing_Clin_v1(casedate, hemisphere)


%% Optional: Process and visualize movement timeseries data cleaned by artifact rejection function

run_MovementProcessing_Clin_artifactRejection(casedate, hemisphere)


%% 3) stats

run_MovementStats_Clin_v1(casedate, hemisphere)

% fix lines 87-91 in ttest_plot_2States_RandTrim s.t. output files save
% within specified outDataDir folder (movementStats)





