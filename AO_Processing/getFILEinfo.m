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
