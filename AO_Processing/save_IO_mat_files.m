function mat_filelist = save_IO_mat_files(studyID)

    % load XLSX file location
    xlsxLoc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    cd(xlsxLoc) 

    % load summaryXLSX table (save in GitHub repo)
    summaryXLSX = readtable("Subject_AO.xlsx");

    % query rows in summaryXLSX by specific pt hemisphere/studyID (single integer correlated with StudyNum)
    filteredXLSX = summaryXLSX(summaryXLSX.StudyNum == studyID, :);

    % identify rows with non-empty/non-NAN cells in the trialNum column 
    trial_rows = ~isnan(filteredXLSX.trialNum); % logical operation

    % extract relevant .mat filenames in the ao_MAT_file column
    mat_filelist = filteredXLSX.ao_MAT_file(trial_rows); % correspond with non-NAN/non-empty cells in the trialNum column
    % output = cell array of relevant .mat filenames

    % convert cell array to table
    T_mat_filelist = cell2table(mat_filelist, 'VariableNames', {'MAT_filenames'});

    % specify output filename
    out_mat_filelist = fullfile(xlsxLoc, 'mat_filelist.xlsx');

    % write table to an Excel file
    writetable(T_mat_filelist, out_mat_filelist);
    
end
