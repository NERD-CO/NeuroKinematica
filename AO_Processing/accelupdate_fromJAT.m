


% if data xyz

% if not data x


if matches(fTYPE1,'CACC')

    if matches(fTYPE2,'DATA')

        % Get the Sensor
        vLISTa = extractBefore(extractAfter(vLIST,21),2);
        sensorList = cellfun(@(x) str2double(x), vLISTa, 'UniformOutput',true);
        sensorLog = sensorList == wholeEleID;
        xyzList = extractAfter(vLIST,25);

        xyzIDs = {'X','Y','Z'};
        freqItemS = cell(1,3);
        for xi = 1:length(xyzIDs)
            xyz = contains(xList,xyzIDs{xi});
            fieldLIST = vLIST(sensorLog & xyz);
            fList_acceli = fieldLIST{matches(extractAfter(fieldLIST,25),xyzIDs{xi})};
            freqItemS{xi} = fList_acceli;
        end

    else
        % Get the Sensor
        vLISTa = extractBefore(extractAfter(vLIST,21),2);
        sensorList = cellfun(@(x) str2double(x), vLISTa, 'UniformOutput',true);
        sensorLog = sensorList == wholeEleID;
        xList = extractAfter(vLIST,25);
        xLog = contains(xList,'X');
        fieldLIST = vLIST(sensorLog & xLog);
      
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

end