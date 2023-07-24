%% Function outline

% Goal: 
    % query rows in summaryXLSX by specific pt hemisphere/studyID (single integer correlated with StudyNum)
    % identify rows with non-empty/non-NAN cells in the trialNum column 
    % extract relevant .mat filenames in the ao_MAT_file column (% output cell array of relevant .mat filenames)
    % extract Ephys fields of interest in relevant .mat files
    % restructure relevant data into output matrix

% Required resources
    % Summary XLSX file 
        % trial IDs and mat filenames
    % Directories for raw data locs
    % Directories for saving post-processed data

%% pseudocode
% function [] = save_IO_EphysProcFiles(IO_DataDir,summaryXLSX,studyID,Case_DataDir,RawDataDir,ProcDataDir)

% hardcode directories
IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';              % directory where all IO data is located
Case_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\03_09_2023'; % directory where case-specific data files are located 
RawDataDir = [Case_DataDir, filesep, 'Raw Electrophysiology MATLAB'];                                % directory where raw MATLAB data files are located (case-specific)
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology'];                                % directory where processed MATLAB data should be saved (case-specific)

% load XLSX file location
xlsxLoc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
cd(xlsxLoc)

% load summaryXLSX table (save in GitHub repo)
summaryXLSX = readtable("Subject_AO.xlsx");

% isolate a specific studyID
studyID = 1;

%% call functions

% extract relevant .mat filenames in the ao_MAT_file column
mat_filelist = save_IO_mat_filenames(studyID);

% extract relevant info from relevant .mat files in mat_filelist
save_IO_mat_ProcFiles(mat_filelist, Case_DataDir);
 

%% main functions


