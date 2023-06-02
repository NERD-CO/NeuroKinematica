% read in csv file
% loop through depths of interest per stn region

csvLoc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
cd(csvLoc) 

summaryCSV = readtable("Subject_AO.csv");

% isolate a specific subject
studyID = 1;
studyDataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\03_09_2023\Raw Electrophysiology MATLAB';

studyTable = summaryCSV(ismember(summaryCSV.StudyNum,studyID),:);
studyTblIndex = ismember(summaryCSV.StudyNum,studyID); % logical index = loc of subject

% create and extract list of unique stn locations
stn_locs = unique(studyTable.stn_loc);

% loop through stn locations
for sti = 1:length(stn_locs)
    % create counter
    keepCount = 1;
    % hold current location
    temp_loc = stn_locs{sti};
    % find relevant rows of table
    stnlTable = studyTable(matches(studyTable.stn_loc,temp_loc),:);
    stnlTblIndex = matches(summaryCSV.stn_loc,temp_loc); % logical index of stn depth
    sumSTNIndex = studyTblIndex & stnlTblIndex; % links subject with stn depth
    nameParts = cellfun(@(x) split(x,'.'),stnlTable.ao_MAT_file,'UniformOutput', false);
    fileOrder = cellfun(@(x) str2double(x{2}(end)),nameParts, 'UniformOutput', true);
    % loop through files per stn location
    for stf = 1:height(stnlTable)
        temp_file = stnlTable.ao_MAT_file{stf};
        % find loc of temp file
        fileTblIndex = matches(summaryCSV.ao_MAT_file,temp_file);
        temp_dir = [studyDataDir,filesep,temp_file];
        load(temp_dir)
        % do we care about this depth?
        matTemp = temp_dir;
        matftemp1 = whos(matfile(matTemp));
        matVarList = {matftemp1.name}; % extract columns of cell array

        ttlCHECK = matches('CDIG_IN_1_KHz',matVarList); % logical - if 1, we care ..maybe.
        % How many ttls?
        if ttlCHECK
           ttl_num =  length(CDIG_IN_1_Up);
           sec_thresh = 60*30;
           if ttl_num < sec_thresh
              summaryCSV.trialNum(fileTblIndex) = NaN;
           else
               % populate row with ID
               summaryCSV.trialNum(fileTblIndex) = keepCount;
               keepCount = keepCount +1;
           end
        else
            summaryCSV.trialNum(fileTblIndex) = NaN;
        end

    end
end

% save new CSV with trial ID
% cd to CSV save loc
% writetable(summaryCSV,'Subject_AO.csv") % fill trial column

stopTest = 1;
