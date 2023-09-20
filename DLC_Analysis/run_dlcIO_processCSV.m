%% Script for running 'dlcIO_processCSV' funtion

% folder for all processed dlc cases 
outer_Loc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\post_DLCproc'; 

% define case-specicfc file locations
vid_dir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\post_DLCproc\03_09_2023\video folder'; % vid_dir
csv_dir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\post_DLCproc\03_09_2023\csv folder'; % csv_dir 
save_matLoc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\post_DLCproc\03_09_2023\mat folder'; % mat_dir

% change dir to case-specific csv folder
% cd csv_dir
cd 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\post_DLCproc\03_09_2023\csv folder'; % csv_dir

% create cell array of all avail csv files
csv_array = dir('*.csv'); % *list all possible .csv filetypes in dir (as struct)
csv_array2 = {csv_array.name}; % outputs list of cell array of the csv filenames


% loop through all avail csv files in csv folder and run 'dlcIO_processCSV' funtion
for i = 1:length(csv_array2)
    csv_vidLoc = [csv_dir, filesep, csv_array2{i}];
    
    
    
    dlcIO_processCSV('dirLOC',1,'userLOC',csv_vidLoc,'saveLOC',save_matLoc)
end

