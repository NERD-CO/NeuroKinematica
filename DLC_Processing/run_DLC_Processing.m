function [] = run_DLC_Processing(Case_DataDir)

% addpath C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\DLC_Processing
addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\DLC_Processing'

DLC_Processing_Dir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\Processed DLC';
cd(DLC_Processing_Dir)

%% Input - define case

casID = 'IO_2023-08-23_RSTN';
Case_DataDir = [DLC_Processing_Dir, filesep, casID];


%% Define case-specific subfolders (csv = input; mat = output)

csv_Loc = [Case_DataDir, filesep, 'csv folder'];
mat_Loc = [Case_DataDir, filesep, 'mat folder'];

% change dir to case-specific csv folder
cd(csv_Loc)

% create cell array of all avail csv files in csv folder
csv_array = dir('*.csv'); % *lists all possible .csv filetypes in dir (as struct)
csv_array2 = {csv_array.name}; % outputs list of the csv filenames (as cell array)


%% Run CSV-to-MAT processing function

% loop through all avail csv files in csv folder and run 'dlcIO_processCSV2' funtion
for csv_i = 1:length(csv_array2)

    tmpCSV = [csv_Loc, filesep, csv_array2{csv_i}];
    dlc_processCSV2('userLOC', tmpCSV, 'saveLOC', mat_Loc)
    
end


%% completed IO (Thesis Aim 1) cases

% Batch 1:
% 03_09_2023 - LSTN (P1, StudyID 1)
% 05_11_2023 - LSTN (P4, StudyID 6)
% 06_08_2023 - LSTN
% 06_08_2023 - RSTN

% Batch 2
% 03_23_2023 - LSTN (P2, StudyID 2)
% 04_05_2023 - RSTN (P2, StudyID 3)
% 04_13_2023 - LSTN (P3, StudyID 4)
% 05_18_2023_b - LSTN


%% completed Clin (Thesis Aim 2, v1) cases 

% 09_12_2023 - LSTN (P7)
    % vid_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc\09_12_2023_LSTN_v2\video folder'; % vid_dir 
    % csv_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc\09_12_2023_LSTN_v2\csv folder'; % csv_dir (input)
    % mat_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc\09_12_2023_LSTN_v2\mat folder'; % mat_dir (output)

% 09_12_2023 - RSTN (P7)
    % vid_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc\09_12_2023_RSTN_v1\video folder'; % vid_dir 
    % csv_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc\09_12_2023_RSTN_v1\csv folder'; % csv_dir (input)
    % mat_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc\09_12_2023_RSTN_v1\mat folder'; % mat_dir (output)

%% completed Clin MDTII cases

% MDTII_MDT5_LSTN
    % vid_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc\MDTII_5_kg\MDTII_MDT5_LSTN\video folder';
    % csv_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc\MDTII_5_kg\MDTII_MDT5_LSTN\csv folder';
    % mat_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc\MDTII_5_kg\MDTII_MDT5_LSTN\mat folder';

% MDTII_MDT5_RSTN
    % vid_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc\MDTII_5_kg\MDTII_MDT5_RSTN\video folder';
    % csv_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc\MDTII_5_kg\MDTII_MDT5_RSTN\csv folder';
    % mat_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc\MDTII_5_kg\MDTII_MDT5_RSTN\mat folder';


end

