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