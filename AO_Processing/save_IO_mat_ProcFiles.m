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
    ftypes = {'CSPK', 'CLFP', 'CMacro_LFP', 'CDIG', 'CEMG', 'CACC'};

    % isolate fields of interest per filetype
    for f = 1:6
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
                % Debugging 'otherwise' case
            
            case 'CEMG'
                % 5. Find all EMG files and extract relevant fields via getFILEinfo
                % save to outStruct
                ProcEphys.EMG = outStruct;

            case 'CACC'
                % 6. Find all accelerometry files and extract relevant fields via getFILEinfo
                % save to outStruct
                ProcEphys.ACC = outStruct;
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

function [outStruct] = getFILEinfo(fTYPE, varITEMS, mfname)

switch fTYPE

    case {'CSPK','CLFP','CMacro_LFP', 'CEMG', 'CACC' } % add EMG and ACC ftypes to this case later
        indiCES = contains(varITEMS, fTYPE); % find relevant fields of ftype
        % Get list
        varLIST = varITEMS(indiCES);
        % how many micro electrodes
        % Get info between '-'
        eleContent = extract(varLIST,digitsPattern); % digitsPatterns - look for all pos. dig. in string
        eleIDs = extractAfter(unique(eleContent),1);
        wholeEleID = unique(eleContent);

        % Debugging - print eleIDs
        disp(['eleIDs: ', strjoin(eleIDs, ', ')]);
        % return an empty struct if eleIDs is empty
        if isempty(eleIDs)
            outStruct = struct;
            return
        end

        for ei = 1:numel(wholeEleID) % numel - # of elements in array
            % Hz
            [freqItem] = getVARid(varLIST, wholeEleID{ei}, fTYPE, '_KHz');
            % outStruct.(['E',num2str(eleIDs{ei})]).Hz = load(mfname,freqItem);
            if ~isempty(freqItem)
                outStruct.(['E',num2str(eleIDs{ei})]).Hz = load(mfname,freqItem);
            else
                warning('freqItem is empty. Skipping load operation for E%d.Hz.', ei);
            end

            % Raw data
            [dataItem] = getVARid(varLIST, wholeEleID{ei}, fTYPE, '');
            % outStruct.(['E',num2str(eleIDs{ei})]).rawData = load(mfname,dataItem);
            if ~isempty(dataItem)
                outStruct.(['E',num2str(eleIDs{ei})]).rawData = load(mfname,dataItem);
            else
                warning('dataItem is empty. Skipping load operation for E%d.rawData.', ei);
            end

            % Start time
            [startTitem] = getVARid(varLIST, wholeEleID{ei}, fTYPE, '_TimeBegin');
            % outStruct.(['E',num2str(eleIDs{ei})]).startTime = load(mfname,startTitem);
            if ~isempty(startTitem)
                outStruct.(['E',num2str(eleIDs{ei})]).startTime = load(mfname,startTitem);
            else
                warning('startTitem is empty. Skipping load operation for E%d.startTime.', ei);
            end

            % End time
            [endTitem] = getVARid(varLIST, wholeEleID{ei}, fTYPE, '_TimeEnd');
            % outStruct.(['E',num2str(eleIDs{ei})]).endTime = load(mfname,endTitem);
            if ~isempty(endTitem)
                outStruct.(['E',num2str(eleIDs{ei})]).endTime = load(mfname,endTitem);
            else
                warning('endTitem is empty. Skipping load operation for E%d.endTime.', ei);
            end
        end

    case 'CDIG'
        indiCES = contains(varITEMS, fTYPE);
        % Get list
        varLIST = varITEMS(indiCES);

        % Debugging
        if isempty(varLIST)
            disp('Warning: Variable list (vLIST) is empty for CDIG file type. Skipping operations...');
            outStruct = struct;
            return
        end

        % Hz
        [freqItem] = getVARid(varLIST, 'IN_1', fTYPE, '_KHz');
        % outStruct.TTL.Hz = load(mfname,freqItem);
        [freqItem] = getVARid(varLIST , 'IN_1' , fTYPE, '_KHz');
        if ~isempty(freqItem)
            outStruct.TTL.Hz = load(mfname,freqItem);
        else
            warning('freqItem is empty. Skipping load operation for TTL.Hz.');
        end

        % Down
        [freqItem] = getVARid(varLIST, 'IN_1', fTYPE, '_Down');
        % outStruct.TTL.Down = load(mfname,freqItem);
        if ~isempty(freqItem)
            outStruct.TTL.Down = load(mfname,freqItem);
        else
            warning('freqItem is empty. Skipping load operation for TTL.Down.');
        end

        % Start time
        [freqItem] = getVARid(varLIST, 'IN_1', fTYPE, '_TimeBegin');
        % outStruct.TTL.startTime = load(mfname,freqItem);
        if ~isempty(freqItem)
            outStruct.TTL.startTime = load(mfname,freqItem);
        else
            warning('freqItem is empty. Skipping load operation for TTL.startTime.');
        end

        % End time
        [freqItem] = getVARid(varLIST, 'IN_1', fTYPE, '_TimeEnd');
        % outStruct.TTL.endTime = load(mfname,freqItem);
        if ~isempty(freqItem)
            outStruct.TTL.endTime = load(mfname,freqItem);
        else
            warning('freqItem is empty. Skipping load operation for TTL.endTime.');
        end

end
end


function [varNAME] = getVARid(vLIST, wholeEleID, fTYPE1, fTYPE2)

% Debugging
disp(['vLIST: ', strjoin(vLIST, ', ')]);
disp(['wholeEleID: ', wholeEleID]);
disp(['fTYPE1: ', fTYPE1]);
disp(['fTYPE2: ', fTYPE2]);

freqItem = vLIST(matches(vLIST,[fTYPE1, '_', wholeEleID, fTYPE2]));

if isempty(freqItem)
    warning('Variable list is empty. Returning empty varNAME');
    varNAME = [];
else
    varNAME = freqItem{1};
end

end
