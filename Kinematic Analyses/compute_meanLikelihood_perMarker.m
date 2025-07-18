function [] = compute_meanLikelihood_perMarker(mainDir, case_hem)

%% Compute mean likelihood per marker 

% Notes-to-self: Data Dirs (for Karen)
% GDrive (mainDLC, MDTII)
% Local (MDTII, MDT 5) - 1 Representative Subject
% % Raw: Z:\RadcliffeE\Transfer_KarenG\KG_DLC Models_new\MDT5
% % Processing: Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc\MDTII_5_kg
% % Analysis: Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\Kinematic Analyses

% specify directory where case-specific data files are located
curPCname = getenv('COMPUTERNAME');

switch curPCname
    case 'DESKTOP-I5CPDO7'  % PC_1
        Clin_Kin_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\Kinematic Analyses';
    case 'DSKTP-JTLAB-EMR'  % Lab Desktop
        Clin_Kin_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\Kinematic Analyses';
    case 'NSG-M-H8J3X34'    % PC_2
        Clin_Kin_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\Kinematic Analyses';
end

cd(Clin_Kin_DataDir)
 
%%% Per Project: Create a metaData tracking sheet for subject data in clinical context %%%
% Subject_AO = readtable('Subject_AO.xlsx');


%% Analyze participant data isolated by caseID and hemisphere

% Define caseID (date or subjectID) and hemisphere

% case_hem = '09_12_2023_LSTN';
% case_hem = '09_12_2023_RSTN';

  case_hem = 'MDTII_MDT5_LSTN';
% case_hem = 'MDTII_MDT5_RSTN';

caseDir = [Clin_Kin_DataDir , filesep , case_hem];
cd(caseDir)


%% Isolate dlc outputs of interest

% Generate list of dlc-video-labeled CSV files
mainCSV = dir('*.csv');
mainCSV2 = {mainCSV.name};

% Generate list of dlc-video-labeled MAT files
mainMat = dir('*.mat');
mainMAT2 = {mainMat.name};

% Generate list of Motor Index CSVs (filters for CSVs that contain 'Move' string)
moveCSV = mainCSV2(contains(mainCSV2,'Move'));


%% Main function loop - Compute and Save Average Likelihood per Marker

% create an outputs directory
outputDir = [caseDir filesep 'processedMovement' filesep 'likelihood_perMarker'];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Initialize overall results table for condition-wise statistics
Likelihood_PerCondition = table();

% Loop through CSV files - Raw Data Processing
for csv_i = 1:length(moveCSV)

    tmpCSV = moveCSV{csv_i};

    % Split file names to extract relevant parts (dateID, sessID, and hemID)
    nameParts = split(tmpCSV,'_');
    dateID = nameParts{1};
    sessID = nameParts{3};
    hemID = nameParts{8};
    matName_title = [dateID , '-' , sessID, '-', hemID]

    % Find and load corresponding dlcDAT MAT file
    matTempfind = [dateID , '_' , sessID];
    matInd = contains(mainMAT2 , matTempfind);
    matName = mainMAT2{matInd};
    load(matName, 'outDATA'); % Load only the required matData

    % Perform artifact rejection for accepting/rejecting frames (default CT = 0.5) 
    CT = 0.5; % confidence threshold (adjustable)
    [outDATA_NaN] = artifactRejection(outDATA, CT);
    % Use outDATA_NaN instead of outDATA for further processing

    % Process dlcDAT MAT file (all points, all frames) per vid first (Split column names of outDATA)
    colNames = outDATA_NaN.Properties.VariableNames; % outDATA_NaN should be a table containing labeled coordinate data from DeepLabCut w/ x,y, values replaced by NaN for marker likelihood values < CT (processed)
    % colNames = outDATA.Properties.VariableNames;
    colNames2 = cellfun(@(x) split(x,'_'), colNames, 'UniformOutput',false);
    uniqueMarkers = unique(cellfun(@(x) x{1}, colNames2, 'UniformOutput', false));
    uniqueMarkers = uniqueMarkers(~matches(uniqueMarkers, 'frames')); % Remove 'frames' column if present

    % Load MoveIndex CSV file to extract specific portions of kinematic timeseries data from dlcDAT MAT file
    moveINDtab = readtable(tmpCSV);
    moveINDtab = moveINDtab(moveINDtab.BeginF ~= 0 & moveINDtab.EndF ~= 0, :); % Remove invalid rows

    % Align start and stop frames / Define frame range for analysis
    firstBegin = moveINDtab.BeginF(1) - 1; % assigned to value of the first element in the BeginF column of moveINDtab - 1 (because MATLAB uses 1-based indexing)
    lastEnd = moveINDtab.EndF(height(moveINDtab)) - 1; % assigned to value of the last element in the EndF column of moveINDtab - 1
    % lastEnd = moveINDtab.EndF(end) - 1;

    % Compute the average likelihood value per anatomical marker 
    % within the specific portions of kinematic timeseries data (b/t the firstBegin and lastEnd frames)
    % Save the average likelihood value per anatomical marker as a table and in CSV format within outputDir

    % Initialize results table
    Likelihood_resultsTable = table();
    
   % Compute likelihood statistics for each marker
    for i = 1:length(uniqueMarkers)
        markerName = uniqueMarkers{i};
        likelihood_col = strcat(markerName, '_likelihood');

        if ismember(likelihood_col, colNames) % Ensure column exists
            likelihood_values = outDATA_NaN.(likelihood_col)(firstBegin:lastEnd);

            % Compute mean, standard error (SE), and standard deviation (STDev)
            mean_likelihood = mean(likelihood_values, 'omitnan');
            se_likelihood = std(likelihood_values, 'omitnan') / sqrt(length(likelihood_values));
            std_likelihood = std(likelihood_values, 'omitnan');

            % Store in results table
            Likelihood_resultsTable = [Likelihood_resultsTable; table({markerName}, mean_likelihood, se_likelihood, std_likelihood, ...
                          'VariableNames', {'Marker', 'Mean_Likelihood', 'SE_Likelihood', 'STD_Likelihood'})];
        end
    end

    % Save results to CSV
    outputFile = fullfile(outputDir, [matName_title, '_Likelihood_Stats.csv']);
    writetable(Likelihood_resultsTable, outputFile);
    
    % Compute averages across all markers per Condition (tmpCSV)
    avg_mean_likelihood = mean(Likelihood_resultsTable.Mean_Likelihood, 'omitnan');
    avg_se_likelihood = mean(Likelihood_resultsTable.SE_Likelihood, 'omitnan');
    avg_std_likelihood = mean(Likelihood_resultsTable.STD_Likelihood, 'omitnan');

    % Store in overall results table
    Likelihood_PerCondition = [Likelihood_PerCondition; table({tmpCSV}, avg_mean_likelihood, avg_se_likelihood, avg_std_likelihood, ...
                   'VariableNames', {'Condition', 'Avg_Mean_Likelihood', 'Avg_SE_Likelihood', 'Avg_STD_Likelihood'})];

    % Display values in command window
    fprintf('\nCondition: %s\n', tmpCSV);
    fprintf('Average Mean Likelihood: %.4f\n', avg_mean_likelihood);
    fprintf('Average SE Likelihood: %.4f\n', avg_se_likelihood);
    fprintf('Average STD Likelihood: %.4f\n', avg_std_likelihood);
end

% Save overall summary results as CSV
overallResultsFile = fullfile(outputDir, 'Likelihood_PerCondition_Summary.csv');
writetable(Likelihood_PerCondition, overallResultsFile);

cd(outputDir);

end

%% Artifact Rejection

function [outDATA_NaN] = artifactRejection(outDATA, CT)
% Inputs:
% - outDATA: labeled timeseries data from DLC (original)
% - CT: Confidence Threshold for accepting/rejecting frames (default = 0.5)
% Outputs:
% - outDATA_NaN: x,y, values replaced by NaN for marker likelihood values < CT (processed)

% Initialize variables
numMarkers = round(width(outDATA) / 3); % Assuming each marker has x, y, and likelihood columnm
outDATA_NaN = outDATA;

% Loop through markers (13) and evaluate each frame for each marker
for markerIdx = 1:numMarkers
    tempMarker_confidence = table2array(outDATA(:, markerIdx*3));
    confidence_logical = tempMarker_confidence <  CT; % logical index
    % replace x,y, values w/ NaN for marker likelihood values < CT
    outDATA_x2y = (markerIdx*3)-2: (markerIdx*3)-1;
    outDATA_NaN(confidence_logical, outDATA_x2y) = repmat({NaN},sum(confidence_logical),2); % replicates dimensions by given row and collumn size
end

end

