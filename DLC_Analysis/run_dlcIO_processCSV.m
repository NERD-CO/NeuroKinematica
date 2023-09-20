
csv_vidLoc = 'C:\Users\radclier\Downloads\csv_t1';
save_matLoc = 'C:\Users\radclier\Downloads\csv_t1';


% change dir to case-specific csv folder
% create cell array of all avail csv files
% create dir stem (csv_dir = '';)
% loop through all avail csv files in csv folder and run _processCSV

csv_dir = '';
save_matLoc = ''; % mat_dir
cd csv_dir
csv_array = dir('*.csv'); % *list all possible .csv filetypes in dir (as struct)
csv_array2 = {csv_array.name}; % outputs list of cell array of the csv filenames

for i = 1:length(csv_array2)
    csv_vidLoc = [csv_dir, filesep, csv_array2{i}];
   





    dlcIO_processCSV('dirLOC',1,'userLOC',csv_vidLoc,'saveLOC',save_matLoc)

end

