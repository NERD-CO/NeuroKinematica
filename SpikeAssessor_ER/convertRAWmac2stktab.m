function [stkTabRaw, stkTabBP] = convertRAWmac2stktab(tsCheck,rawTS,Fs,timeFlag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Incoming matrix will be oriented by channel x sampling point
% Need to transpose

switch tsCheck
    case 'seeg'
        % check orientation
        if size(rawTS,2) > size(rawTS,1)
            rawTStrk = double(transpose(rawTS));
        else
            rawTStrk = double(rawTS);
        end

        tempTAB = array2table(rawTStrk);

        if timeFlag
            mTime = transpose(seconds((1:height(tempTAB))/Fs));
            tempTAB.TimeSecs = mTime;
            stkTabRaw = table2timetable(tempTAB,'RowTimes','TimeSecs');
        else
            stkTabRaw = tempTAB;
        end

        % Loop through
        % Bipolar mean of neighboring (1-2); mean - 1,

        bpTempTab = nan(height(rawTStrk),width(rawTStrk)-1);
        for ti = 1:width(rawTStrk)-1
            tBPi = double(rawTStrk(:,ti)) - mean(double(rawTStrk(:,ti:ti+1)),2);
            bpTempTab(:,ti) = tBPi;
        end

        tempTAB2 = array2table(bpTempTab);

        if timeFlag
            tempTAB2.TimeSecs = mTime;
            stkTabBP = table2timetable(tempTAB2,'RowTimes','TimeSecs');
        else
            stkTabBP = tempTAB2;
        end


    case 'psd'
        if size(rawTS,2) > size(rawTS,1)
            rawTStrk = transpose(rawTS);
        else
            rawTStrk = rawTS;
        end

        tempTAB = array2table(rawTStrk);

        stkTabRaw = tempTAB;
        stkTabBP = nan;

end



end