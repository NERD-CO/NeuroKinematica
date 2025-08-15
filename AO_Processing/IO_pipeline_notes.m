%% IO pipeline notes - chronologically order tasks and functions

%% Metadata Org. (itemized steps in sheet 2: Data Tracking)

% Box --> Synology --> Data Architecture (in Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative)

% 1) IO_DataSummary_Radcliffe.xlsx (OneDrive)
% sheet 1: Summary
% sheet 2: Data Tracking
% sheet 3: Subject_AO
% sheet 4: Subject_DLC
% sheet 5: Demographics
% sheet 6: BCI 2025

% 2) Subject_AO.xlsx (in Synology ... \\som-nsg-r-ao1)
% sheet 3 in (1) is filled out


%% MATLAB workflow

% AO_Processing (GitHub\NeuroKinematica\AO_Processing)
%%% consider combining the following 2 scripts into 1 function %%%
% > trial_ID_determination
% > save_IO_EphysProcFiles
    % contains functions (2)
    % > [mat_filelist, ACC_check] = save_IO_mat_filenames(studyID, IO_DataDir); (1)
    % > save_IO_mat_ProcFiles(mat_filelist, Case_DataDir, ACC_check); (2)


% Spike Clustering (GitHub\NeuroKinematica\SpikeAssessor_ER)
% > CluterReviewer.mlapp


% FR Analysis (GitHub\NeuroKinematica\IO_FR_Analysis)
% > spikes_Align_TTL_and_AO
% > compute_meanFR_perMove_perSTNdepth_v2


% FR Raster Plotting


% LFP Analysis






%% Unfamiliar stuff in AO_Processing repo - review with JAT

% > updateTrialInformation (similar to trial_ID_determination - choose/use  1)
% > script2assess_updateTrialInfo

% > matfileTTLcheck
% > accelupdate_fromJAT
% > accelCheck_script
% > accelAOcheck

% New stuff
% > save_IO_mat_ProcFiles_genericLOC

% Older stuff (?)
% getLISTof_recDepths_er
% save_DLCprocFiles_er

