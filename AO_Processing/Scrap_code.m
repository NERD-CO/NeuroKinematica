%% Scrap code

% navigate to raw MAT file data
cd(RawDataDir)
Case_Matfile = dir('*.mat');             % get list of all .mat files in the current directory
Case_MatfileNames = {Case_Matfile.name}; % create cell array containing names of the .mat files found in previous step

% create struct containing only the .mat files in IO_Matfile (list of all .mat files in the studyID directory) that match the filenames in mat_filelist

% initialize an empty struct to store matched files
MATmatchedFiles = struct;

% iterate over each file in mat_filelist
for i = 1:height(mat_filelist)

    % get current filename from the table
    currentFile = mat_filelist.MAT_filenames{i};

    % Check if the current file exists in IO_MatfileNames
    % if ismember(currentFile, Case_MatfileNames)
        % If there is a match, load the mat file and save it in the struct
        % Replace '.' with '_' to create a valid field name
        fieldName = strrep(currentFile, '.', '_');
        MATmatchedFiles.(fieldName) = load(currentFile);
    % end
end