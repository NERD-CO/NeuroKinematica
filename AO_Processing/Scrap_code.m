%% Scrap code

% navigate to raw MAT file data
cd(RawDataDir)
Case_Matfile = dir('*.mat');             % get list of all .mat files in the current directory
Case_MatfileNames = {Case_Matfile.name}; % create cell array containing names of the .mat files found in previous step

% query rows in summaryXLSX by specific pt hemisphere/studyID (single integer correlated with StudyNum)
% filteredXLSX = summaryXLSX(summaryXLSX.StudyNum == studyID, :);

% identify rows with non-empty/non-NAN cells in the trialNum column 
% trial_rows = ~isnan(filteredXLSX.trialNum); % logical operation
% trial_rows = ~cellfun(@isempty, filteredXLSX.trialNum) & ... ~cellfun(@(x) any(isnan(x)), filtered_table.trialNum);

% extract relevant .mat filenames in the ao_MAT_file column
% mat_filelist = filteredXLSX.ao_MAT_file(trial_rows); % correspond with non-NAN/non-empty cells in the trialNum column
% output = cell array of relevant .mat filenames

% create struct containing only the .mat files in IO_Matfile (list of all .mat files in the studyID directory) that match the filenames in mat_filelist

% initialize an empty struct to store matched files
MATmatchedFiles = struct;

% iterate over each file in mat_filelist
for i = 1:height(mat_filelist)

    % get current filename from the table
    currentFile = mat_filelist.MAT_filenames{i};

    % Check if the current file exists in IO_MatfileNames
    % if ismember(currentFile, Case_MatfileNames)
        % If there is a match, load the mat file and save it in the struct
        % Replace '.' with '_' to create a valid field name
        fieldName = strrep(currentFile, '.', '_');
        MATmatchedFiles.(fieldName) = load(currentFile);
    % end
end