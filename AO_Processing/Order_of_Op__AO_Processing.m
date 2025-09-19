%% Order_of_Op__AO_Processing

clear; close all; clc;

%% functions:

% > trial_ID_determination
% > save_IO_EphysProcFiles
    % contains functions (2)
    % > [mat_filelist, ACC_check] = save_IO_mat_filenames(studyID, IO_DataDir); (1)
    % > save_IO_mat_ProcFiles(mat_filelist, Case_DataDir, ACC_check); (2)

%%% consider combining the 2 scripts ^ into 1 function %%%


%% Environment / Directory Set-up

% specify directory where case-specific data files are located
curPCname = getenv('COMPUTERNAME');

switch curPCname
    case 'DESKTOP-I5CPDO7'  % PC_1
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'DSKTP-JTLAB-EMR'  % Lab Desktop
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\AO_Processing'
    case 'NSG-M-H8J3X34'    % PC_2
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        addpath 'C:\GitHub\NeuroKinematica\AO_Processing'
end

cd(IO_DataDir)
Subject_AO = readtable('Subject_AO.xlsx');