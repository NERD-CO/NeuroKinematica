%% Script / function outline

% function [] = save_IO_EphysProcFiles(studyID, Case_DataDir)

% Goal: 
    % query rows in summaryXLSX by specific pt hemisphere/studyID (single integer correlated with StudyNum)
    % identify rows with non-empty/non-NAN cells in the trialNum column 
    % extract relevant .mat filenames in the ao_MAT_file column (output cell array of relevant .mat filenames)
    % extract Ephys fields of interest in relevant .mat files
    % restructure relevant data into output matrix
    % save processed ephys data matrix in new folder for each case

% Required resources
    % Summary XLSX file 
        % trial IDs and mat filenames
    % Directories for raw data locs
    % Directories for saving post-processed data


%% CCC

clc
close 
clear

%% Machine-specific data directory Input

% specify directory where case-specific data files are located 
curPCname = getenv('COMPUTERNAME'); % for windows

switch curPCname
    case 'DESKTOP-I5CPDO7' 
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';  
    case 'NSG-M-H8J3X34'   % EMR PC (Precision 5490)
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative'; 
    case 'DSKTP-JTLAB-EMR' % EMR Lab Desktop (DigitalStorm Lumos)
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';  
end

%% Subject-specific data directory Inputs

% isolate a specific studyID and CaseDate
studyID = 29;
CaseDate = '12_06_2023_bilateral';

Case_DataDir = [IO_DataDir, filesep, CaseDate];


%% Completed cases:

% 1: '03_09_2023'
% 2: '03_23_2023'
% 3: '04_05_2023'
% 4: '04_13_2023'
% 5: '04_13_2023'
% 6: '05_11_2023'; *ACC
% 7: '05_18_2023_a'; *ACC
% 8: '05_18_2023_b'
% 9: '05_18_2023_b'
% 10: '05_31_2023'; *ACC
% 11: '06_08_2023'; *ACC
% 12: '06_08_2023'; *ACC
% 13: '07_06_2023_bilateral';
% 14: '07_06_2023_bilateral'; *ACC
% 15: '07_13_2023_bilateral'; *ACC
% 16: '07_13_2023_bilateral'; *ACC
% 17: '07_20_2023'; *ACC
% 18: '07_26_2023'; *ACC
% 19: '08_10_2023_bilateral'; 
% 20: '08_10_2023_bilateral'; 
% 21: '08_23_2023'
% 22: '08_24_2023
% 23: '10_04_2023'
% 24: '11_01_2023' *ACC
% 25: '11_08_2023' *ACC
% 26: '11_09_2023' *ACC
% 27: '11_30_2023_bilateral' *ACC
% 28: '11_30_2023_bilateral' *ACC
% 29: '12_06_2023_bilateral' *ACC
% 30: '12_06_2023_bilateral' *ACC


%% Hardcode directories

% directory where all IO data is located
RawDataDir = [Case_DataDir, filesep, 'Raw Electrophysiology MATLAB'];       % directory where raw MATLAB data files are located (case-specific)
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology'];       % directory where processed MATLAB data should be saved (case-specific)
if ~exist(ProcDataDir,'dir')
    mkdir(ProcDataDir);
end

% load XLSX file location
cd(IO_DataDir)

% load summaryXLSX table (save in GitHub repo)
summaryXLSX = readtable("Subject_AO.xlsx");


%% call main function 1

% extract relevant .mat filenames in the ao_MAT_file column
[mat_filelist, ACC_check] = save_IO_mat_filenames(studyID, IO_DataDir);


%% call main function 2

% extract relevant info from relevant .mat files in mat_filelist
save_IO_mat_ProcFiles(mat_filelist, Case_DataDir, ACC_check);

cd(IO_DataDir);

%% Ephys Processing done :)
