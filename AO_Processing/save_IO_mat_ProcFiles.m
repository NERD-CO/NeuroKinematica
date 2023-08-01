function [] = save_IO_mat_ProcFiles(mat_filelist, Case_DataDir, ACC_IO)

% hardcode directories
RawDataDir = [Case_DataDir, filesep, 'Raw Electrophysiology MATLAB']; % directory where raw MATLAB data files are located (case-specific)
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology']; % directory where processed MATLAB data should be saved (case-specific)

% navigate to directory where raw MATLAB data files are located
cd(RawDataDir)

% extract relevant info from relevant .mat files in mat_filelist
for i = 1:height(mat_filelist)

    % Added 7/31/2023
    cd(RawDataDir)

    tmpFilename = mat_filelist{i};          % retrieve i-th .mat filename from mat_filelist
    matFileInfo = matfile(tmpFilename);     % create matfile object representing the .mat file specified by tmpFilename
    matFileVars1 = whos(matFileInfo);       % use 'whos' function to get info about variables stored in .mat file represented by matFileInfo
    matFileVars2 = {matFileVars1.name};     % extract names of all variables in the .mat file and store them in cell array matFileVars2.

    % list Ephys filetypes of interest (conditions)
    if ACC_IO(i) == 0
        ftypes = {'CSPK', 'CLFP', 'CMacro_LFP', 'CDIG'}; % 'CEMG'
    else
        ftypes = {'CSPK', 'CLFP', 'CMacro_LFP', 'CDIG', 'CACC'}; % 'CEMG'
    end
    % isolate fields of interest per filetype
    for f = 1:length(ftypes) % 1:6
        Ftype = ftypes{f};

        % use getFILEinfo function to process the data based on the ftype - look at code in GitHub repo: save_DLCprocFiles_er
        outStruct = getFILEinfo(Ftype, matFileVars2, tmpFilename, ACC_IO(i)); % create new struct containing fields of interest

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

                % case 'CEMG'
                %     % 5. Find all EMG files and extract relevant fields via getFILEinfo
                %     % save to outStruct
                %     ProcEphys.EMG = outStruct;

                %
            case 'CACC'
                % 6. Find all accelerometry files and extract relevant fields via getFILEinfo
                % save to outStruct
                ProcEphys.ACC = outStruct;
        end

    end

    % save into new directory with new name
    saveName = ['Processed_',tmpFilename];
    cd(ProcDataDir)
    save(saveName,'ProcEphys');

end
end


%% helper functions

function [outStruct] = getFILEinfo(fTYPE, varITEMS, mfname, accelCount)

switch fTYPE

    case {'CSPK','CLFP','CMacro_LFP', 'CEMG'} % add EMG and ACC ftypes to this case later
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
            [varFreq] = extractVarName(mfname,freqItem);
            outStruct.(['E',num2str(eleIDs{ei})]).Hz = varFreq;

            % Raw data
            [dataItem] = getVARid(varLIST, wholeEleID{ei}, fTYPE, '');
            [varData] = extractVarName(mfname,dataItem);
            outStruct.(['E',num2str(eleIDs{ei})]).rawData = varData;

            % Start time
            [startTitem] = getVARid(varLIST, wholeEleID{ei}, fTYPE, '_TimeBegin');
            [varStime] = extractVarName(mfname,startTitem);
            outStruct.(['E',num2str(eleIDs{ei})]).startTime = varStime;
            % End time

            [endTitem] = getVARid(varLIST, wholeEleID{ei}, fTYPE, '_TimeEnd');
            [varEtime] = extractVarName(mfname,endTitem);
            outStruct.(['E',num2str(eleIDs{ei})]).endTime = varEtime;
        end

    case 'CDIG'
        indiCES = contains(varITEMS, fTYPE);
        % Get list
        varLIST = varITEMS(indiCES);
        % Hz
        [freqItem] = getVARid(varLIST, 'IN_1', fTYPE, '_KHz');
        [varFreq] = extractVarName(mfname,freqItem);
        outStruct.Hz = varFreq;

        % Down
        [downItem] = getVARid(varLIST, 'IN_1', fTYPE, '_Down');
        [varDown] = extractVarName(mfname,downItem);
        outStruct.Down = varDown;

        % Start time
        [startTitem] = getVARid(varLIST, 'IN_1', fTYPE, '_TimeBegin');
        [varStime] = extractVarName(mfname,startTitem);
        outStruct.startTime = varStime;

        % End time
        [endTitem] = getVARid(varLIST, 'IN_1', fTYPE, '_TimeEnd');
        [varEtime] = extractVarName(mfname,endTitem);
        outStruct.endTime = varEtime;

    case 'CACC'
        indiCES = contains(varITEMS, fTYPE); % find relevant fields of ftype
        % Get list
        varLIST = varITEMS(indiCES);
        % how many micro electrodes
        % Get info between '-'
        elePresent1 = cellfun(@(x) x(22:end), varLIST, 'UniformOutput', false);
        elePresent2 = cellfun(@(x) str2double(x(1)), elePresent1, 'UniformOutput', true);
        AccelAvailable = 1:accelCount;
        eleIDs = elePresent2(ismember(elePresent2, AccelAvailable));
        wholeEleID = unique(eleIDs);

        outStruct = struct;

        for ei = 1:numel(wholeEleID) % numel - # of elements in array
            % Hz
            [freqItem] = getVARid(varLIST, wholeEleID(ei), fTYPE, '_KHz');
            [varFreq] = extractVarName(mfname,freqItem{1});
            outStruct.(['ACC',num2str(wholeEleID(ei))]).Hz = varFreq;

            % Raw data
            [dataItem] = getVARid(varLIST, wholeEleID(ei), fTYPE, 'DATA');
            % outStruct.(['E',num2str(eleIDs{ei})]).rawData = load(mfname,dataItem);

            % loop through 3 axes of each accel
            for di = 1:length(dataItem)
                tmpLoadF = load(mfname,dataItem{di});
                tmpLoadFns = fieldnames(tmpLoadF);
                outStruct.(['ACC',num2str(wholeEleID(ei))]).rawData(di,:) = tmpLoadF.(tmpLoadFns{1});
            end

            % Start time
            [startTitem] = getVARid(varLIST, wholeEleID(ei), fTYPE, '_TimeBegin');
            [varStime] = extractVarName(mfname,startTitem{1});
            outStruct.(['ACC',num2str(wholeEleID(ei))]).startTime = varStime;

            % End time
            [endTitem] = getVARid(varLIST, wholeEleID(ei), fTYPE, '_TimeEnd');
            [varEtime] = extractVarName(mfname,endTitem{1});
            outStruct.(['ACC',num2str(wholeEleID(ei))]).endTime = varEtime;

        end
end
end


function [varNAME] = getVARid(vLIST, wholeEleID, fTYPE1, fTYPE2)

% Debugging
% disp(['vLIST: ', strjoin(vLIST, ', ')]);
% disp(['wholeEleID: ', wholeEleID]);
% disp(['fTYPE1: ', fTYPE1]);
% disp(['fTYPE2: ', fTYPE2]);

if matches(fTYPE1,'CACC')

    if matches(fTYPE2,'DATA')

        % Get the Sensor
        vLISTa = extractBefore(extractAfter(vLIST,21),2);
        sensorList = cellfun(@(x) str2double(x), vLISTa, 'UniformOutput',true);
        sensorLog = sensorList == wholeEleID;
        xyzList = extractAfter(vLIST,25);

        xyzIDs = {'X','Y','Z'};
        freqItem = cell(1,3);
        for xi = 1:length(xyzIDs)
            xyz = contains(xyzList,xyzIDs{xi});
            fieldLIST = vLIST(sensorLog & xyz);
            fList_acceli = fieldLIST{matches(extractAfter(fieldLIST,25),xyzIDs{xi})};
            freqItem{xi} = fList_acceli;
        end

    else

        % Get X
        if wholeEleID == 1
            sensorID = '01';
            eleID = '1';
        else
            sensorID = '04';
            eleID = '2';
        end
        tempFINDname = ['CACC_3___',sensorID,'___Sensor_',eleID , '___X',fTYPE2];

        freqItem = vLIST(matches(vLIST,tempFINDname));
    end

varNAME = freqItem;

else
    freqItem = vLIST(matches(vLIST,[fTYPE1, '_', wholeEleID, fTYPE2]));

    if isempty(freqItem)
        warning('Variable list is empty. Returning empty varNAME');
        varNAME = [];
    else
        varNAME = freqItem{1};
    end
end


end



function [varExtract] = extractVarName(inMfname,inMvar)

tmpLoadF = load(inMfname,inMvar);
tmpLoadFns = fieldnames(tmpLoadF);
varExtract = tmpLoadF.(tmpLoadFns{1});

end
