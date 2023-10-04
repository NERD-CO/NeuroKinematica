%% Script for running 'dlcIO_processCSV2' funtion

% folder for all processed dlc cases 
outer_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\post_DLCproc'; 

% define case-specific file locations
vid_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\post_DLCproc\03_09_2023\video folder'; % vid_dir 
csv_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\post_DLCproc\03_09_2023\csv folder'; % csv_dir (input)
mat_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\post_DLCproc\03_09_2023\mat folder'; % mat_dir (output)

% change dir to case-specific csv folder
cd(csv_Loc)

% create cell array of all avail csv files in csv folder
csv_array = dir('*.csv'); % *lists all possible .csv filetypes in dir (as struct)
csv_array2 = {csv_array.name}; % outputs list of the csv filenames (as cell array)


% loop through all avail csv files in csv folder and run 'dlcIO_processCSV2' funtion
for csv_i = 1:length(csv_array2)

    tmpCSV = [csv_Loc, filesep, csv_array2{csv_i}];
    dlcIO_processCSV2('userLOC', tmpCSV, 'saveLOC', mat_Loc)
    
end


%% completed cases

% 03_09_2023
% vid_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\post_DLCproc\03_09_2023\video folder'; % vid_dir 
% csv_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\post_DLCproc\03_09_2023\csv folder'; % csv_dir (input)
% mat_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\post_DLCproc\03_09_2023\mat folder'; % mat_dir (output)

