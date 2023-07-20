function mat_filelist = save_IO_mat_files(studyID)

    % hardcode directories
    IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    RawDataDir = [IO_DataDir, filesep, 'Raw Electrophysiology MATLAB']; % filesep saves based OS syntax
    ProcDataDir = [IO_DataDir, filesep, 'Processed Electrophysiology'];

    % load summaryXLSX table (save in GitHub repo)
    summaryXLSX = readtable("Subject_AO.xlsx");

    % query rows in summaryXLSX by specific pt hemisphere/studyID (single integer correlated with StudyNum)
    filteredXLSX = summaryXLSX(summaryXLSX.StudyNum == studyID, :);

    % identify rows with non-empty/non-NAN cells in the trialNum column 
    trial_rows = ~isnan(filteredXLSX.trialNum); % logical operation

    % extract relevant .mat filenames in the ao_MAT_file column
    mat_filelist = filteredXLSX.ao_MAT_file(trial_rows); % correspond with non-NAN/non-empty cells in the trialNum column

    % output cell array of relevant .mat filenames
    if isnumeric(mat_filelist) % check the type of mat_files, if it's a cell array no need to convert
        mat_filelist = num2cell(mat_filelist); % if it's a numeric vector, convert it to cell array
    end

end
