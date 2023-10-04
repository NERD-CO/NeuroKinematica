


pcname = getenv('COMPUTERNAME');
hemisphere = 'L';

switch pcname
    case 'DESKTOP-I5CPDO7'  %%% JAT

        mainDir = 'W:\RadcliffeE\INS_2024';

    case ''     %%% ER

        mainDir = '';

end


switch hemisphere
    case 'L'

        mainDir2 = [mainDir , filesep , 'LSTN'];

    case 'R'

        mainDir2 = [mainDir , filesep , 'LSTN'];

end


cd(mainDir2)

% List of CSV files
mainCSV = dir('*.csv');
mainCSV2 = {mainCSV.name};

% List of Motor Index CSVs
moveCSV = mainCSV2(contains(mainCSV2,'Move'));

% List of MAT motor indicies
mainMotorIND = dir('*.mat');
mainMotorIND2 = {mainMotorIND.name};



% Loop through CSV files
for ci = 1:length(moveCSV)

    tmpCSV = moveCSV{ci};

    nameParts = split(tmpCSV,'_');

    dateID = nameParts{1};
    sessID = nameParts{3};

    matTempfind = [dateID , '_' , sessID];

    matMotInd = contains(mainMotorIND2 , matTempfind);

    matMotName = mainMotorIND2{matMotInd};

    load(matMotName)
    % process main MAT first

    colNames = outDATA.Properties.VariableNames;
    colNames2 = cellfun(@(x) split(x,'_'), colNames,...
        'UniformOutput',false);
    colNames3 = unique(cellfun(@(x) x{1}, colNames2,...
        'UniformOutput',false));
    colNames4 = colNames3(~matches(colNames3,'frames'));

    euclidall = zeros(height(outDATA)-1,length(colNames4));

    for labI = 1:length(colNames4)

        tmpLabel_x = [colNames4{labI} , '_x'];
        tmpLabel_y = [colNames4{labI} , '_y'];

        tmpXdata = outDATA.(tmpLabel_x);
        tmpYdata = outDATA.(tmpLabel_y);

        labelData = [tmpXdata , tmpYdata];

        for fraI = 1:height(labelData)
            if fraI ~= height(labelData)

                point1 = labelData(fraI,:);
                point2 = labelData(fraI + 1,:);
                euclidall(fraI , labI) = pdist2(point1 , point2);

            end
        end

    end

    % 2. Reduce to Finger tips
    ftipInds = contains(colNames4,'fTip');
    ftipEuclid = euclidall(:,ftipInds);

    % 3. Average Finger tips
    ftipAverage = mean(ftipEuclid,2);

    % 4. Index using motor index
    moveINDtab = readtable(tmpCSV);
    % clean up
    moveINDtab = moveINDtab(~moveINDtab.BeginF == 0,:);

    % First Begin 
    firstBegin = moveINDtab.BeginF(1) - 1; % To align with Euclidean distance frames
    lastEnd = moveINDtab.EndF(height(moveINDtab)) - 1;


    % FingerTip average block

    ftipAveBlk = ftipAverage(firstBegin:lastEnd);

    % 1. Use smoothdata function 'guassian' smooth out the edges
    % 2. compute the half width -- peakfind [documention prominences]
    % 3. compute amplitude of peaks
    % 4. compute peak distance
    







end






















