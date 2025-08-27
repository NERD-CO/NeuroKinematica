%% IO pipeline notes - chronologically order tasks and functions

%% Metadata Org. (itemized steps in sheet 2: Data Tracking)

% Box --> Synology --> Data Architecture (in Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative)

% IO_DataSummary_Radcliffe.xlsx (OneDrive)
    % sheet 1: Summary
    % sheet 2: Data Tracking
    % sheet 3: Subject_AO
    % sheet 4: Subject_DLC
    % sheet 5: Demographics
    % sheet 6: BCI 2025

% Subject_AO.xlsx (in Synology ... \\som-nsg-r-ao1)
    % sheet 3 in IO_DataSummary_Radcliffe.xlsx ^


%% MATLAB workflow

% 1) Ephys Processing (GitHub\NeuroKinematica\AO_Processing)                % done!
% > trial_ID_determination
% > save_IO_EphysProcFiles
    % contains functions (2)
    % > [mat_filelist, ACC_check] = save_IO_mat_filenames(studyID, IO_DataDir); (1)
    % > save_IO_mat_ProcFiles(mat_filelist, Case_DataDir, ACC_check); (2)
%%% consider combining the 2 scripts ^ into 1 function %%%


% 2) Spike Clustering (GitHub\NeuroKinematica\SpikeAssessor_ER)             % in progress
% > CluterReviewer.mlapp


% 3) DLC Processing (GitHub\NeuroKinematica\DLC_Processing)                 % in progress
% > Order_of_Op__DLC_Processing
    % > run_DLC_Processing
    % > fillDroppedFrames_JT_v2

% 4) Movement Indexing (GitHub\NeuroKinematica\DLC_VideoIndexing_GUI)       % in progress
% > ER_DLC_MoveCheck_Dual_v4.mlapp


% FR Processing (GitHub\NeuroKinematica\IO_FR_Analysis)                     % in progress
% > Order_of_Op__FRKin_Processing
    % > align_SpikesPerMove_TTL


% FR Analysis & Plotting (GitHub\NeuroKinematica\IO_FR_Analysis)            % in progress
% > Order_of_Op__FRKin_Analyses
    % > compute_FRperMove_perSTNdepth
    % > [FR_SummaryTbl] = run_IO_FR_Analysis_and_Plotting
    % > [kinTbl, kinSummaryTbl] = run_MovementFeatureAnalysis_IO_v2
    % > merge_FRKin_SummaryTbls
    % > aggregate_FRKinematic_Correlations


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

