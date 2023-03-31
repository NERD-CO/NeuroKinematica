function [] = save_DLCprocFiles_er(dlcLISTdir , rawMATdir , dlcDATAdir)
%% Navigate to folder with DLC file list table
cd(dlcLISTdir)
dlcMatfile = dir('*.mat');
dlcfileName = {dlcMatfile.name};
dlcfileName = dlcfileName{1};

load(dlcfileName,'outTableF')


%% Loop through file list - it is a table
for mi = 1:height(outTableF)

    cd(rawMATdir)

    tmpFname = outTableF.FullFile{mi};
    matFileInfo = matfile(tmpFname);
    matFileVars1 = whos(matFileInfo);
    matFileVars2 = {matFileVars1.name};

    %% Use matfile function to
    % 1. Find all spike files
    [outStructSPK] = getFILEinfo('CSPK',matFileVars2,tmpFname);
    % 2. Find all LFP
    [outStructLFP] = getFILEinfo('CLFP',matFileVars2,tmpFname);
    % 3. Find all mLFP
    [outStructMLFP] = getFILEinfo('CMacro_LFP',matFileVars2,tmpFname);
    % 4. Find all TTL
    [outStructTTL] = getFILEinfo('CDIG',matFileVars2,tmpFname);
    % Optional Find EMG when used

    % Save into one Struct
    dlcDepths.Spike = outStructSPK;
    dlcDepths.LFP = outStructLFP;
    dlcDepths.MLFP = outStructMLFP;
    dlcDepths.TTL = outStructTTL;

    % Save into new directory with new name
    sAVEname = ['DLCao_',tmpFname];
    cd(dlcDATAdir)
    save(sAVEname,'dlcDepths');


end

end



function [outStruct] = getFILEinfo(fTYPE,varITEMS,mfname)

switch fTYPE
    case {'CSPK','CLFP','CMacro_LFP'}
        indiCES = contains(varITEMS,fTYPE);
        % Get list
        varLIST = varITEMS(indiCES);
        % how many micro electrodes
        % Get info between '-'
        eleContent = extract(varLIST,digitsPattern);
        eleIDs = extractAfter(unique(eleContent),1);
        wholeEleID = unique(eleContent);
        
        for ei = 1:numel(wholeEleID)

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

        % cell array of stuff to save

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

switch fTYPE1
    case {'CSPK','CLFP','CMacro_LFP'}
        freqItem = vLIST(matches(vLIST,[fTYPE1,'_', wholeEleID ,fTYPE2]));
        varNAME = freqItem{1};


    case 'CDIG'
        freqItem = vLIST(matches(vLIST,[fTYPE1,'_', wholeEleID ,fTYPE2]));
        varNAME = freqItem{1};

end


end

