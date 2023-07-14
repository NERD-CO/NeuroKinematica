%% Function outline

% Goal: 
    % query by pt/study ID
    % extract fields of interest in mat file
    % restructure relevant data into matrix

% Required resources

% Summary XLSX file 
    % mat filenames and trial IDs
% Directories for raw data locs
% Directories for saving post-processed data

%% pseudocode
% load XLSX file location
% 1: create >= 1 input arg. for pt ID (single integer correlated with StudyNum)
    % extract rows in mat file 

% summaryXLSX saved in GitHub repo
summaryXLSX = readtable("Subject_AO.xlsx");

% hardcode directories
SubjectDataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
RawDataDir = [SubjectDataDir, filesep, 'Raw Electrophysiology MATLAB']; % filesep saves based OS syntax
ProcDataDir = [SubjectDataDir, filesep, 'Processed Electrophysiology'];

% ID non NAN rows in col. 6 - create logical vec
% with relevant rows (1s), extract through relevant mat file names in col. 7 (output: cell array of mat files)
% loop through relevant mat files

% combine RawDataDir with filename to creat load location
% load() % matfile
% *isolate fields of interest - look at code in GitHub repo
% create new struct containing fields of interest
% save it: save(filename, new struct var) % filename ~ ProcDataDir combined with mat fileame


