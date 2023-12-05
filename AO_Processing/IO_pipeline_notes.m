% %% IO pipeline notes - chronologically order tasks and functions

%% AO_Processing

% Current workflow
% > trial_ID_determination
% > save_IO_EphysProcFiles
    % contains functions (2)
    % > [mat_filelist, ACC_check] = save_IO_mat_filenames(studyID, IO_DataDir); (1)
    % > save_IO_mat_ProcFiles(mat_filelist, Case_DataDir, ACC_check); (2)


% Less familiar stuff - review with JAT
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

