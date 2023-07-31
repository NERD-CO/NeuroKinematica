function mat_filelist = save_IO_mat_filenames(studyID , xlsxLoc)

% load XLSX file location
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

end
