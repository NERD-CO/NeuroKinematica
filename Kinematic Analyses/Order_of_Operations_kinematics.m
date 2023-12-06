%% Order of Operations

%% 1) Define function inputs


casedate = '09_12_2023';
hemisphere = 'L';


%% 2) Process and visualize movement timeseries data

run_MovementProcessing_v1(casedate, hemisphere)


%% Optional: Process and visualize movement timeseries data cleaned by artifact rejection function

run_MovementProcessing_artifactRejection(casedate, hemisphere)


%% 3) stats

run_MovementStats_v1(casedate, hemisphere)

% fix lines 87-91 in ttest_plot_2States_RandTrim s.t. output files save
% within specified outDataDir folder (movementStats)





