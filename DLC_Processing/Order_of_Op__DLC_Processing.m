%% Order of Operations - IO DLC Data Processing

clear; clc;

% addpath C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\DLC_Processing
addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\DLC_Processing'

% addpath 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\Kinematic Analyses'
addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\Kinematic Analyses'

% Folder for all processed DLC cases 
IO_procDLC = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\Processed DLC';
Clin_procDLC = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc';

% Folder for all kinematic analyses  
IO_kinematicData = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\Kinematic Analyses'; 
Clin_kinematicData = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\Kinematic Analyses';

cd(IO_procDLC);
% cd(Clin_procDLC);

%% 1) Define DLC processing function inputs

% Define casedate and hemisphere:
casedate_hem = 'IO_05_31_2023_LSTN';

% 'IO_03_09_2023_RSTN';  
% 'IO_03_23_2023_LSTN';
% 'IO_04_05_2023_RSTN';
% 'IO_04_13_2023_LSTN';
% 'IO_05_18_2023_a_RSTN';
% 'IO_05_31_2023_LSTN';
% 'IO_06_08_2023_LSTN';
% 'IO_2023_07_06_LSTN';
% 'IO_2023_07_13_LSTN';
% 'IO_2023_07_13_RSTN';
% 'IO_2023_08_23_RSTN';

% Define case-specific kinematic data dir
DLCProc_caseID = [IO_procDLC, filesep, casedate_hem];


%% 2) Run DLC Processing function

DLCProc_dir = 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\DLC_Processing';
cd(DLCProc_dir)

% Convert DLC timeseries data per video from .csv to .mat
run_DLC_Processing(IO_procDLC, DLCProc_caseID)


%% 3) Fill dropped frames to align front and side cam frames in video folder

cd(DLCProc_caseID)
path2videos = [DLCProc_caseID, filesep, 'video folder'];
quality = 100; 

% run JAT video processing function
fillDroppedFrames_JT_v2(path2videos, quality)


%% 4) Generate Movement Indices for each video

MovIndexGUI_dir = 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\DLC_VideoIndexing_GUI';
cd(MovIndexGUI_dir)

% ER_DLC_MoveCheck_Dual_v4.mlapp


%% Next Steps

% Gor to Order of Operation Scripts in IO_FR_Analysis
FRKin_dir = 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis';
cd(FRKin_dir)

% > Order_of_Op__FRKin_Processing
% > Order_of_Op__FRKin_Analyses

% run everything up to:
% [kinTbl, kinSummaryTbl] = run_MovementFeatureAnalysis_IO_v2(IO_DataDir, MoveDataDir, MoveDir_CaseID);


%% Stats













