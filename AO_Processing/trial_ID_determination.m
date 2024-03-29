
% navigate to summaryXLSX file location (IO_DataDir)
xlsxLoc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
cd(xlsxLoc) 

% read in summaryXLSX table 
summaryXLSX = readtable("Subject_AO.xlsx");

%% Inputs: isolate a specific subject/case date
studyID = 18;
studyMatDataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\07_26_2023\Raw Electrophysiology MATLAB';

%% Completed subjects/cases:
% 1: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\03_09_2023\Raw Electrophysiology MATLAB'
% 2: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\03_23_2023\Raw Electrophysiology MATLAB'
% 3: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\04_05_2023\Raw Electrophysiology MATLAB'
% 4: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\04_13_2023\Raw Electrophysiology MATLAB\LH'
% 5: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\04_13_2023\Raw Electrophysiology MATLAB\RH'
% 6: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\05_11_2023\Raw Electrophysiology MATLAB'
% 7: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\05_18_2023_a\Raw Electrophysiology MATLAB'
% 8: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\05_18_2023_b\Raw Electrophysiology MATLAB\LH'
% 9: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\05_18_2023_b\Raw Electrophysiology MATLAB\RH'
% 10: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\05_31_2023\Raw Electrophysiology MATLAB'
% 11: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\06_08_2023\Raw Electrophysiology MATLAB'
% 12: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\06_08_2023\Raw Electrophysiology MATLAB'
% 13: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\07_06_2023_bilateral\Raw Electrophysiology MATLAB'
% 14: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\07_06_2023_bilateral\Raw Electrophysiology MATLAB'
% 15: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\07_13_2023_bilateral\Raw Electrophysiology MATLAB'
% 16: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\07_13_2023_bilateral\Raw Electrophysiology MATLAB'
% 17: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\07_20_2023\Raw Electrophysiology MATLAB'
% 18: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\07_26_2023\Raw Electrophysiology MATLAB'

%% function

% Define data table and indexing variables:
studyTable = summaryXLSX(ismember(summaryXLSX.StudyNum,studyID),:);
studyTblIndex = ismember(summaryXLSX.StudyNum,studyID); % logical index = loc of subject

% create and extract list of unique stn locations
stn_locs = unique(studyTable.stn_loc);

% loop through depths of interest per stn region
for sti = 1:length(stn_locs)
    % create counter
    keepCount = 1;
    % hold current location
    temp_loc = stn_locs{sti};
    % find relevant rows of table
    stnlTable = studyTable(matches(studyTable.stn_loc,temp_loc),:);

    % stnlTblIndex = matches(summaryCSV.stn_loc,temp_loc); % logical index of stn depth
    % sumSTNIndex = studyTblIndex & stnlTblIndex; % links subject with stn depth
    % nameParts = cellfun(@(x) split(x,'.'),stnlTable.ao_MAT_file,'UniformOutput', false);
    % fileOrder = cellfun(@(x) str2double(x{2}(end)),nameParts, 'UniformOutput', true);

    % loop through files per stn location
    for stf = 1:height(stnlTable)
        temp_file = stnlTable.ao_MAT_file{stf};
        % find loc of temp file
        fileTblIndex = matches(summaryXLSX.ao_MAT_file,temp_file); % notes row to save relvant experimental rec. ID 
        temp_dir = [studyMatDataDir,filesep,temp_file];
        % load(temp_dir)

        % do we care about this depth?
        matftemp = whos(matfile(temp_dir)); % look at filenames without loading file content
        matVarList = {matftemp.name}; % extract columns of cell array (filenames)
        ttlCHECK = matches('CDIG_IN_1_KHz',matVarList); % logical - if 1, we care ..maybe.
        tfilenum = str2double(extractBefore(extractAfter(temp_file,'F'), '.'));
        if height(stnlTable) == 4 && tfilenum == 1
            firstfilecheck = 1;
        else
            firstfilecheck = 0;
        end

        % How many ttls?
        if firstfilecheck
            summaryXLSX.trialNum(fileTblIndex) = NaN;
        else

            if ttlCHECK
                load(temp_dir,"CDIG_IN_1_Down")
                ttl_num =  length(CDIG_IN_1_Down) % 60 frames per sec.
                ttl_thresh = 60*25; % ((30 sec, 1800 ttls; 28 sec, 1680 ttls; 25 sec, 1500 ttls) 
                if ttl_num < ttl_thresh
                    summaryXLSX.trialNum(fileTblIndex) = NaN;
                else
                    % populate row with ID
                    summaryXLSX.trialNum(fileTblIndex) = keepCount;
                    keepCount = keepCount +1;
                end
            else
                summaryXLSX.trialNum(fileTblIndex) = NaN; % why is this condition here?
            end
        end

    end
end

% save new CSV with trial ID
cd(xlsxLoc) 
writetable(summaryXLSX,'Subject_AO.xlsx') % fill trial column

% output: csv file with trialNum column filled with relevant experimental iteration (per STN location)

