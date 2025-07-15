%% Goal: Update Subject_AO excel file with trialNum iteration per STN depth based on defined criteria

% specify directory where case-specific data files are located
curPCname = getenv('COMPUTERNAME');

switch curPCname
    case 'DESKTOP-I5CPDO7'  % PC_1
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'DSKTP-JTLAB-EMR'  % Lab Desktop
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'NSG-M-H8J3X34'    % PC_2
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
end

% navigate to summaryXLSX file location (IO_DataDir)
cd(IO_DataDir)

% read in summaryXLSX table (Subject_AO)
summaryXLSX = readtable("Subject_AO.xlsx");


%% Inputs: isolate a specific subject/case date

studyID = 20;
studyMatDataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\08_10_2023_bilateral\Raw Electrophysiology MATLAB';


%% Define Case-specific Data Directories %%% next step for streamlining script into function %%%

% Case Date = 
% Case_DataDir = [IO_DataDir, filesep, CaseDate];
% cd(Case_DataDir)

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
cd(IO_DataDir) 
writetable(summaryXLSX,'Subject_AO.xlsx') % fill trial column

%%

% output: updated Subject_AO.xlsx file with trialNum column filled with relevant/qualifying experimental iteration per STN location

