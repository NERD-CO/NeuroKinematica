function [] = preView_MoveProcess_IO_jatv1(mainDir , csvChoose)

% Analyze data isolated by casedate and hemisphere

mainDirVID = [mainDir , filesep , 'video folder'];
mainDirMAT = [mainDir , filesep , 'mat folder'];
mainDirCSV = [mainDir , filesep , 'csv folder'];

% Isolate dlc outputs of interest
cd(mainDirCSV)
% Generate list of dlc-video-labeled CSV files
mainCSV = dir('*.csv');
rawCSV = {mainCSV.name};
cd(mainDirMAT)
% Generate list of dlc-video-labeled MAT files
mainMat = dir('*.mat');
mainMAT2 = {mainMat.name};
cd(mainDirVID)
mainCSVb = dir('*.csv');
mainCSVb2 = {mainCSVb.name};
% Generate list of Motor Index CSVs (filters for CSVs that contain 'Move' string)
moveCSV = mainCSVb2(contains(mainCSVb2,'Move'));

switch csvChoose
    case 'move'

        useCSV = moveCSV;

    case 'raw'

        useCSV = rawCSV;

end

% EUC indicies
cd(mainDir)


%% Main function
% create an outputs directory
outputDir = [mainDir , filesep , 'processedMovement'];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
% Define framerate of videos (time conversion factor)
% fps = 60; % frames per second

% Convert distance units to mm (distance conversion factor)
pixels_to_mm = 2.109; % 232 mm / 110 pxl = 2.1091 mm per pixel
% Anthropometry: vertical distance from the bottom of the chin (menton) to the top of the head: https://upload.wikimedia.org/wikipedia/commons/0/06/AvgHeadSizes.png
% US adult male, 50th percentile: Avg. = 23.2 cm, 9.1 inches


% Loop through CSV files - Raw Data Processing
for csv_i = 1:length(useCSV)

    tmpCSV = useCSV{csv_i};

    % Split file names to extract relevant parts (dateID, sessID, and hemID)
    nameParts = split(tmpCSV,'_');
    dateID = nameParts{1};
    sessID = nameParts{3};
    camID = nameParts{5};
    runID = nameParts{2};
    % hemID = nameParts{8};
    % matName_title = [dateID , '-' , sessID, '-', hemID];

    % Find and load corresponding dlcDAT MAT file
    matTempfind = [dateID , '_' , sessID,'_',runID,'_',camID];
    matInd = contains(mainMAT2 , matTempfind);
    matName = mainMAT2{matInd};
    cd(mainDirMAT)
    load(matName , 'outDATA')

    % Call artifact rejection function
    % CT = 0.5; % Set confidence threshold (CT) for accepting/rejecting frames (default CT = 0.5), adjust
    % [outDATA_NaN] = artifactRejection(outDATA, CT); % Use outDATA_NaN instead of outDATA for further processing

    % Process dlcDAT MAT file (all points, all frames) per vid first (Split column names of outDATA)
    colNames = outDATA.Properties.VariableNames; % outDATA should be a table containing labeled coordinate data from DeepLabCut
    colNames2 = cellfun(@(x) split(x,'_'), colNames,...
        'UniformOutput',false);
    colNames3 = unique(cellfun(@(x) x{1}, colNames2,...
        'UniformOutput',false));
    colNames4 = colNames3(~matches(colNames3,'frames'));

    % Initialize 'euclidall' to store Euclidean distances between successive points
    euclidall = zeros(height(outDATA)-1,length(colNames4));

    % Iterate over each label and compute Euclidean distance for each frame
    for label_i = 1:length(colNames4)

        tmpLabel_x = [colNames4{label_i} , '_x'];
        tmpLabel_y = [colNames4{label_i} , '_y'];

        tmpXdata = outDATA.(tmpLabel_x);
        tmpYdata = outDATA.(tmpLabel_y);

        labelData = [tmpXdata , tmpYdata];

        for frame_i = 1:height(labelData)
            if frame_i ~= height(labelData)
                point1 = labelData(frame_i,:);
                point2 = labelData(frame_i + 1,:);
                euclidall(frame_i , label_i) = pdist2(point1 , point2);
            end
        end
    end

    % Convert distance variables to mm usng conversion factor
    euclidall = euclidall * pixels_to_mm; % converting euclidean distances to mm

    % Filter the computed distances related to fingertip movements
    tIME = transpose(seconds(0:1/60:(length(euclidall)-1)/60));

    euclidTT = array2timetable(euclidall,'RowTimes',tIME);
    stackedplot(euclidTT)

    tmpCSVpng = extractBefore(tmpCSV,'.');
    % fileNAME = [moveTypeIDs{mmi}, ' ', sessID , ' ' , runID];
    title(tmpCSVpng)
    cd(mainDirCSV)


    saveas(gcf,tmpCSVpng,'png')



end




end


