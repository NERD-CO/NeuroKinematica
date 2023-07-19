%% Function outline

% Goal: 
    % query by pt/study ID
    % extract fields of interest in mat file
    % restructure relevant data into matrix

%% Required resources

% Summary XLSX file 
    % trial IDs and mat filenames
% Directories for raw data locs
% Directories for saving post-processed data

%% Scrap (?):

% load XLSX file location
xlsxLoc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
cd(xlsxLoc)

%% pseudocode

% isolate a specific studyID
studyID = 1;

% hardcode directories
IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
RawDataDir = [IO_DataDir, filesep, 'Raw Electrophysiology MATLAB']; % filesep saves based OS syntax
ProcDataDir = [IO_DataDir, filesep, 'Processed Electrophysiology'];

% load summaryXLSX table (save in GitHub repo)
summaryXLSX = readtable("Subject_AO.xlsx");

% query rows in summaryXLSX by specific pt hemisphere/studyID (single integer correlated with StudyNum)
filteredXLSX = summaryXLSX(summaryXLSX.StudyNum == studyID, :);

% identify rows with non-empty/non-NAN cells in the trialNum column 
trial_rows = ~isnan(filteredXLSX.trialNum); % logical operation
% trial_rows = ~cellfun(@isempty, filteredXLSX.trialNum) & ... ~cellfun(@(x) any(isnan(x)), filtered_table.trialNum);

% extract relevant .mat filenames in the ao_MAT_file column
mat_files = filteredXLSX.ao_MAT_file(trial_rows); % correspond with non-NAN/non-empty cells in the trialNum column

% output cell array of relevant .mat filenames
if isnumeric(mat_files) % check the type of mat_files, if it's a cell array no need to convert
    mat_files = num2cell(mat_files); % if it's a numeric vector, convert it to cell array
end


% loop through relevant mat files in cell array
    % combine RawDataDir with filename to creat load location
    % load() % matfile
for i = 1:height(mat_files)
    % use fuction(s) to isolate fields of interest 
        % look at code in GitHub repo: save_DLCprocFiles_er
        % 1. Find all spike files 'CSPK'
        % 2. Find all LFP 'CLFP'
        % 3. Find all mLFP 'CMacro_LFP'
        % 4. Find all TTL 'CDIG'
        % Optional Find EMG when used

    % create new struct containing fields of interest

    % save into one struct

    % save into new directory with new name

    % save(filename, new struct var) % filename ~ ProcDataDir combined with mat fileame
end

% save it: 

%% functions

function mat_files = save_IO_mat_files(studyID)

    % hardcode directories
    IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    RawDataDir = [IO_DataDir, filesep, 'Raw Electrophysiology MATLAB']; % filesep saves based OS syntax
    ProcDataDir = [IO_DataDir, filesep, 'Processed Electrophysiology'];

    % load summaryXLSX table (save in GitHub repo)
    summaryXLSX = readtable("Subject_AO.xlsx");

    % query rows in summaryXLSX by specific pt hemisphere/studyID (single integer correlated with StudyNum)
    filteredXLSX = summaryXLSX(summaryXLSX.StudyNum == studyID, :);

    % identify rows with non-empty/non-NAN cells in the trialNum column 
    trial_rows = ~isnan(filteredXLSX.trialNum); % logical operation

    % extract relevant .mat filenames in the ao_MAT_file column
    mat_files = filteredXLSX.ao_MAT_file(trial_rows); % correspond with non-NAN/non-empty cells in the trialNum column

    % output cell array of relevant .mat filenames
    if isnumeric(mat_files) % check the type of mat_files, if it's a cell array no need to convert
        mat_files = num2cell(mat_files); % if it's a numeric vector, convert it to cell array
    end

end