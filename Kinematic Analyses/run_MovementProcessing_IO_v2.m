function [] = run_MovementProcessing_IO_v2(mainDir, casedate_hem)

% Goal: Process and visualize movement timeseries data based on videos that have been anatomically labeled (13pt per frame) and analyzed via a trained DeepLabCut model

%% Analyze data isolated by casedate and hemisphere

% Define casedate and hemisphere

% casedate_hem = '03_09_2023_RSTN';

mainDir2 = [mainDir , filesep , casedate_hem];
cd(mainDir2)

%%
mainDir_VID = [mainDir2 , filesep , 'video folder'];
mainDir_MAT = [mainDir2 , filesep , 'mat folder'];
mainDir_CSV = [mainDir2 , filesep , 'csv folder'];


%% Isolate dlc outputs of interest

% Generate list of dlc-video-labeled CSV files
cd(mainDir_CSV)
% mainCSV = dir('*.csv');
% mainCSV2 = {mainCSV.name};

% Generate list of dlc-video-labeled MAT files
cd(mainDir_MAT)
mainMat = dir('*.mat');
mainMAT2 = {mainMat.name};

cd(mainDir_VID)
mainCSVb = dir('*.csv');
mainCSVb2 = {mainCSVb.name};

% Generate list of Motor Index CSVs (filters for CSVs that contain 'Move' string)
moveCSV = mainCSVb2(contains(mainCSVb2,'Move'));

% EUC indicies
cd(mainDir2)
eucINDICIES = readtable("EUC_Indicies.xlsx");


%% Define time conversion factor and distance conversion factor

% Define framerate of videos (time conversion factor)
fps = 60; % frames per second

% Convert distance units to mm (distance conversion factor)
pixels_to_mm = 2.109; % 232 mm / 110 pxl = 2.1091 mm per pixel
% Anthropometry: vertical distance from the bottom of the chin (menton) to the top of the head: https://upload.wikimedia.org/wikipedia/commons/0/06/AvgHeadSizes.png
% US adult male, 50th percentile: Avg. = 23.2 cm, 9.1 inches
% Subject in video frames: Avg. = 110 pixels

% *** define distance conversion factor on standardized calibration measure
% moving forward


%% Define/initialize variables for recording conditions

% STN depths (ideally 3)
% 1) dorsal STN (t1, t2, t3):  3 sessions
% 2) central STN (c1, c2, c3): 3 sessions
% 3) ventral STN (b1, b2, b3): 3 sessions

%% Main function

% create an outputs directory
outputDir = [mainDir2 filesep 'processedMovement'];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Loop through CSV files - Raw Data Processing
for csv_i = 1:length(moveCSV)

    tmpCSV = moveCSV{csv_i};

    % Split file names to extract relevant parts (dateID, sessID, and hemID)
    nameParts = split(tmpCSV,'_');
    dateID = nameParts{1};
    sessID = nameParts{3};
    camID = nameParts{5};
    runID = nameParts{2};
    hemID = nameParts{8};
    matName_title = [dateID , '-' , sessID, '-', hemID];

    % Find and load corresponding dlcDAT MAT file
    matTempfind = [dateID , '_' , sessID,'_',runID,'_',camID];
    matInd = contains(mainMAT2 , matTempfind);
    matName = mainMAT2{matInd};
    cd(mainDir_MAT)
    load(matName , 'outDATA')

    % Call artifact rejection function
    CT = 0.5; % Set confidence threshold (CT) for accepting/rejecting frames (default CT = 0.5), adjust
    [outDATA_NaN] = artifactRejection(outDATA, CT); % Use outDATA_NaN instead of outDATA for further processing

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

    curVID = eucINDICIES(matches(eucINDICIES.videoID,tmpCSV),:);

    handOCindsT = curVID(matches(curVID.moveID,'HAND OC'),:);
    handOCinds = mean(euclidall(:,handOCindsT.eucID),2);
    handPSindsT = curVID(matches(curVID.moveID,'HAND PS'),:);
    handPSinds = mean(euclidall(:,handPSindsT.eucID),2);
    armEFindsT = curVID(matches(curVID.moveID,'ARM EF'),:);
    armEFinds = mean(euclidall(:,armEFindsT.eucID),2);

    % tIME = transpose(seconds(0:1/60:(length(euclidall)-1)/60));

    % euclidTT = array2timetable(euclidall,'RowTimes',tIME)
    %  stackedplot(euclidTT)

    % Average the computed distances related to fingertip movements
    % fTipAverage = mean(fTipEuclid,2);
    % elbowAverage = mean(elBlowEuclid,2);

    % Process dlcDAT MAT files using MoveIndex CSV files to select specific portions of the averaged fingertip distances
    cd(mainDir_VID)
    moveINDtab = readtable(tmpCSV);
    moveINDtab = moveINDtab(~moveINDtab.BeginF == 0,:); % clean up - filters out rows in moveINDtab where the BeginF field is zero.
    moveINDtab = moveINDtab(~moveINDtab.EndF == 0,:); % clean up - filters out rows in moveINDtab where the EndF field is zero.


    % LOOP THROUGH MOVE TYPES
    moveTypeIDs = unique(moveINDtab.MoveType);
    moveTypeNum = numel(moveTypeIDs);

    for mmi = 1:moveTypeNum

        moveTypeTab = moveINDtab(matches(moveINDtab.MoveType,moveTypeIDs{mmi}),:);

        % Align with Euclidean distance frames
        firstBegin = moveTypeTab.BeginF(1) - 30; % assigned to value of the first element in the BeginF column of moveINDtab - 1 (because MATLAB uses 1-based indexing)
        lastEnd = moveTypeTab.EndF(height(moveTypeTab)) + 30; % assigned to value of the last element in the EndF column of moveINDtab - 1

        switch moveTypeIDs{mmi}
            case 'HAND OC'
                moveEUC2use = handOCinds;
                window_Width = 15;
            case 'HAND PS'
                moveEUC2use = handPSinds;
                window_Width = 15;
            case 'ARM EF'
                moveEUC2use = armEFinds;
                window_Width = 25;
        end

        if lastEnd > height(moveEUC2use)
            lastEnd = height(moveEUC2use);
        end

        % Extract and store the average fingertip distance for specified frames (in fTip Average Block)
        moveAveBlk = moveEUC2use(firstBegin:lastEnd); % extracts subset of fTipAverage w/in frame range from firstBegin to lastEnd (represents specific portion of data where specified movement is detected, as indicated in MoveIndex CSV file)

        % Smooth out edges -- smoothdata function w/ 'guassian' method
        % window_Width = 5; % set windowWidth as needed
        smoothed_moveAveBlk = smoothdata(moveAveBlk, 'gaussian', window_Width); % read documentation, window overlap

        switch moveTypeIDs{mmi}
            case 'HAND OC'

                ampMEAN = mean(smoothed_moveAveBlk,'omitnan');
                ampSTD = std(smoothed_moveAveBlk,'omitnan');
                ampThresh = (ampMEAN + (ampSTD/2))*0.4; % 0.75

                [peaks, locs, widths, prominences] = findpeaks(smoothed_moveAveBlk,...
                    MinPeakHeight=mean(smoothed_moveAveBlk), MinPeakDistance=15,...
                    MinPeakProminence=ampThresh, Annotate ='extents');

                close all
                figure;

                findpeaks(smoothed_moveAveBlk,...
                    MinPeakHeight=mean(smoothed_moveAveBlk), MinPeakDistance=15,...
                    MinPeakProminence=ampThresh, Annotate ='extents')

            case 'HAND PS'


                ampMEAN = mean(smoothed_moveAveBlk,'omitnan');
                ampSTD = std(smoothed_moveAveBlk,'omitnan');
                ampThresh = (ampMEAN + (ampSTD/2))*0.3;

                [peaks, locs, widths, prominences] = findpeaks(smoothed_moveAveBlk,...
                    MinPeakHeight=mean(smoothed_moveAveBlk), MinPeakDistance=10,...
                    MinPeakProminence=ampThresh, Annotate ='extents');

                close all
                figure;

                findpeaks(smoothed_moveAveBlk,...
                    MinPeakHeight=mean(smoothed_moveAveBlk), MinPeakDistance=5,...
                    MinPeakProminence=ampThresh, Annotate ='extents')


            case 'ARM EF'
                [peaks, locs, widths, prominences] = findpeaks(smoothed_moveAveBlk,...
                    MinPeakHeight=mean(smoothed_moveAveBlk)*1.2, MinPeakDistance=20, Annotate ='extents');

                close all
                figure;

                findpeaks(smoothed_moveAveBlk,...
                    MinPeakHeight=mean(smoothed_moveAveBlk)*1.2, MinPeakDistance=20, Annotate ='extents')

        end

        xlabel('time (s)');
        ylabel('amplitude');
        hold off
        fileNAME = [moveTypeIDs{mmi}, ' ', sessID , ' ' , runID];
        title(fileNAME)
        cd(outputDir)
        legend('none')

        saveas(gcf,fileNAME,'png')
    end

        % set thresholds based on mean +/- 3xStd
        % MinPeakHeight = mean(smoothed_fTipAveBlk);
        % MinPeakProminence = mean(smoothed_fTipAveBlk)*2;

        % Rubric for movement type -OR- iterative algorthmic approach
        % MinPeakDistance_HandMov =
        % MinPeakDistance_PronSup =
        % MinPeakDistance_FlexExtend =
        % MinPeakDistance_FingerTap =

        % Find peak amplitudes and compute widths -- findpeaks function [review documentation]
        [peaks_fps, locs_fps, widths_fps, prominences_fps] = findpeaks(smoothed_moveAveBlk, fps, MinPeakHeight=mean(smoothed_moveAveBlk), MinPeakDistance=0.15, MinPeakProminence=mean(smoothed_moveAveBlk)*2, Annotate ='extents');


        % Convert distance variables to mm usng distance conversion factor
        amplitudes = peaks_fps * pixels_to_mm; % converting amplitudes to mm

        % Timepoints (in seconds) rather than frame numbers using video sampling rate (Fs) conversion factor
        timepoints_fps = locs_fps; 

        % Compute distances between consecutive peaks
        peakDists_fps = diff(timepoints_fps); % by timepoint (in seconds)

        halfWidths_fps = widths_fps / 2;
        % slope

        % Define unique name for the findpeaks results based on the current CSV name
        findpeaks_output = [outputDir , filesep , 'findpeaks_output_' ,...
            moveTypeIDs{mmi}, '_', tmpCSV(1:end-44) '.csv']; % assumes tmpCSV is a string ending in '.csv'

        % Create a table to store results based on computed variables
        T1 = table(timepoints_fps, locs_fps, peaks_fps, amplitudes, prominences_fps, widths_fps, halfWidths_fps, 'VariableNames', {'Timepoints', 'Locations', 'Peaks', 'Amplitudes', 'Prominences', 'Widths', 'HalfWidths'});

        % Write results table to a CSV file
        writetable(T1, findpeaks_output);

        % Define unique name for the outDATA_NaN output file based on the current CSV/MAT filename
        outDATA_NaN_filename = [outputDir filesep 'cleaned_dlcDAT_NaN_' tmpCSV(1:end-44) '.csv']; % Removes '.csv' from tmpCSV and appends 'outDATA_NaN_'
        % outDATA_NaN_filename = [outputDir filesep 'outDATA_NaN_' tmpCSV(1:end-44) '.csv']; % Removes '.csv' from tmpCSV and appends 'outDATA_NaN_'

        cd(outputDir)
        % Write outDATA_NaN table to a CSV file in the specified output directory
        writetable(outDATA_NaN, outDATA_NaN_filename);

    end

end


%% sub-functions

function [outDATA_NaN] = artifactRejection(outDATA, CT)
% Inputs:
% - outDATA: labeled timeseries data from DLC (original)
% - CT: Confidence Threshold for accepting/rejecting frames (default = 0.5)
% Outputs:
% - outDATA_NaN: x,y, values replaced by NaN for marker likelihood values < CT (processed)

% Initialize variables
numMarkers = round(width(outDATA) / 3); % Assuming each marker has x, y, and likelihood columnm
outDATA_NaN = outDATA;

% Loop through markers (13) and evaluate each frame for each marker
for markerIdx = 1:numMarkers
    tempMarker_confidence = table2array(outDATA(:, markerIdx*3));
    confidence_logical = tempMarker_confidence <  CT; % logical index
    % replace x,y, values w/ NaN for marker likelihood values < CT
    outDATA_x2y = (markerIdx*3)-2: (markerIdx*3)-1;
    outDATA_NaN(confidence_logical, outDATA_x2y) = repmat({NaN},sum(confidence_logical),2); % replicates dimensions by given row and collumn size
end

end
