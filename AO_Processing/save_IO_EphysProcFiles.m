%% Function outline

% Goal: 
    % query by pt/study ID
    % extract fields of interest in mat file
    % restructure relevant data into matrix

% Required resources
    % Summary XLSX file 
        % trial IDs and mat filenames
    % Directories for raw data locs
    % Directories for saving post-processed data

%% pseudocode
% function [] = save_IO_EphysProcFiles(IO_DataDir,summaryXLSX,studyID,Case_DataDir,RawDataDir,ProcDataDir)

% hardcode directories
IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';              % directory where all IO data is located
Case_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\03_09_2023'; % directory where case-specific data files are located 
RawDataDir = [Case_DataDir, filesep, 'Raw Electrophysiology MATLAB'];                                % directory where raw MATLAB data files are located (case-specific)
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology'];                                % directory where processed MATLAB data should be saved (case-specific)

% load XLSX file location
xlsxLoc = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
cd(xlsxLoc)

% load summaryXLSX table (save in GitHub repo)
summaryXLSX = readtable("Subject_AO.xlsx");

% isolate a specific studyID
studyID = 1;

%% call functions

% extract relevant .mat filenames in the ao_MAT_file column
mat_filelist = save_IO_mat_files(studyID);

% extract relevant info from relevant .mat files in mat_filelist
save_IO_mat_ProcFiles(mat_filelist, Case_DataDir);
 

%% main functions

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

end

function [] = save_IO_mat_ProcFiles(mat_filelist, Case_DataDir)

% hardcode directories
RawDataDir = [Case_DataDir, filesep, 'Raw Electrophysiology MATLAB']; % directory where raw MATLAB data files are located (case-specific)
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology']; % directory where processed MATLAB data should be saved (case-specific)

% navigate to directory where raw MATLAB data files are located
cd(RawDataDir)

% extract relevant info from relevant .mat files in mat_filelist
for i = 1:height(mat_filelist)

    tmpFilename = mat_filelist{i};          % retrieve i-th .mat filename from mat_filelist
    matFileInfo = matfile(tmpFilename);     % create matfile object representing the .mat file specified by tmpFilename
    matFileVars1 = whos(matFileInfo);       % use 'whos' function to get info about variables stored in .mat file represented by matFileInfo
    matFileVars2 = {matFileVars1.name};     % extract names of all variables in the .mat file and store them in cell array matFileVars2.

    % list Ephys filetypes of interest (conditions)
    ftypes = {'CSPK', 'CLFP', 'CMacro_LFP', 'CDIG'};
    
    % isolate fields of interest per filetype
    for f = 1:4
        Ftype = ftypes{f};

        % use getFILEinfo function to process the data based on the ftype - look at code in GitHub repo: save_DLCprocFiles_er
        outStruct = getFILEinfo(Ftype, matFileVars2, tmpFilename); % create new struct containing fields of interest

        switch Ftype
            case 'CSPK'
                % 1. Find all spike files and extract relevant fields via getFILEinfo
                % save to outStruct
                ProcEphys.Spike = outStruct;

            case 'CLFP'
                % 2. Find all LFP files and extract relevant fields via getFILEinfo
                % save to outStruct
                ProcEphys.LFP = outStruct;

            case 'CMacro_LFP'
                % 3. Find all mLFP files and extract relevant fields via getFILEinfo
                % save to outStruct
                ProcEphys.MLFP = outStruct;

            case 'CDIG'
                % 4. Find all TTL files and extract relevant fields via getFILEinfo
                % save to outStruct
                ProcEphys.TTL = outStruct;
        end

    end

    % 5. Find EMG when used
    % 6. Find ACC when used

    % save into new directory with new name
    saveName = ['Processed_',tmpFilename];
    cd(ProcDataDir)
    save(saveName,'ProcEphys');

end
end


%% helper functions

function [outStruct] = getFILEinfo(fTYPE,varITEMS,mfname)

switch fTYPE
    case {'CSPK','CLFP','CMacro_LFP'} % add EMG and ACC ftypes to this case later
        indiCES = contains(varITEMS,fTYPE); % find relevant fields of ftype
        % Get list
        varLIST = varITEMS(indiCES);
        % how many micro electrodes
        % Get info between '-'
        eleContent = extract(varLIST,digitsPattern); % digitsPatterns - look for all pos. dig. in string
        eleIDs = extractAfter(unique(eleContent),1);
        wholeEleID = unique(eleContent);

        for ei = 1:numel(wholeEleID) % numel - # of elements in array
            % Hz
            [freqItem] = getVARid(varLIST , wholeEleID{ei} , fTYPE, '_KHz');
            outStruct.(['E',num2str(eleIDs{ei})]).Hz = load(mfname,freqItem);

            % Raw data
            [dataItem] = getVARid(varLIST , wholeEleID{ei} , fTYPE, '');
            outStruct.(['E',num2str(eleIDs{ei})]).rawData = load(mfname,dataItem);

            % Start time
            [startTitem] = getVARid(varLIST , wholeEleID{ei} , fTYPE, '_TimeBegin');
            outStruct.(['E',num2str(eleIDs{ei})]).startTime = load(mfname,startTitem);

            % End time
            [endTitem] = getVARid(varLIST , wholeEleID{ei} , fTYPE, '_TimeEnd');
            outStruct.(['E',num2str(eleIDs{ei})]).endTime = load(mfname,endTitem);
        end

    case 'CDIG'
        indiCES = contains(varITEMS,fTYPE);
        % Get list
        varLIST = varITEMS(indiCES);

        % Hz
        [freqItem] = getVARid(varLIST , 'IN_1' , fTYPE, '_KHz');
        outStruct.TTL.Hz = load(mfname,freqItem);

        % Down
        [freqItem] = getVARid(varLIST , 'IN_1' , fTYPE, '_Down');
        outStruct.TTL.Down = load(mfname,freqItem);

        % Start time
        [freqItem] = getVARid(varLIST , 'IN_1' , fTYPE, '_TimeBegin');
        outStruct.TTL.startTime = load(mfname,freqItem);

        % End time
        [freqItem] = getVARid(varLIST , 'IN_1' , fTYPE, '_TimeEnd');
        outStruct.TTL.endTime = load(mfname,freqItem);
end

end



function [varNAME] = getVARid(vLIST , wholeEleID , fTYPE1, fTYPE2)

freqItem = vLIST(matches(vLIST,[fTYPE1,'_', wholeEleID ,fTYPE2]));
varNAME = freqItem{1};

end

