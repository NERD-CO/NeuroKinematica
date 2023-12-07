function [] = run_MovementProcessing_Clin_v1(casedate, hemisphere)

% Goal: Process and visualize movement timeseries data based on videos that have been anatomically labeled (13pt per frame) and analyzed via a trained DeepLabCut model

%% Directory set-up - Navigate b/t machines
pcname = getenv('COMPUTERNAME');

switch pcname
    case 'DESKTOP-I5CPDO7'   %%% JAT Desktop

        % mainDir = '';

    case 'DSKTP-JTLAB-EMR'   %%% ER Desktop

        mainDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\Kinematic Analyses';

    case 'NSG-M-FQBPFK3'     %%% ER PC

        mainDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\Kinematic Analyses';
end


%% Analyze data isolated by casedate and hemisphere

% Define casedate and hemisphere
% casedate = '09_12_2023';
% hemisphere = 'L';


mainDir2 = [mainDir , filesep , casedate];

switch hemisphere
    case 'L'

        mainDir3 = [mainDir2 , filesep , 'LSTN'];

    case 'R'

        mainDir3 = [mainDir2 , filesep , 'RSTN'];
end

cd(mainDir3)


%% Isolate dlc outputs of interest

% Generate list of dlc-video-labeled CSV files
mainCSV = dir('*.csv');
mainCSV2 = {mainCSV.name};

% Generate list of dlc-video-labeled MAT files
mainMat = dir('*.mat');
mainMAT2 = {mainMat.name};

% Generate list of Motor Index CSVs (filters for CSVs that contain 'Move' string)
moveCSV = mainCSV2(contains(mainCSV2,'Move'));


%% Main function

% create an outputs directory
outputDir = [mainDir3 filesep 'processedMovement'];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Define framerate of videos (time conversion factor)
fps = 60; % frames per second

% Convert distance units to mm (distance conversion factor)
pixels_to_mm = 2.109; % 232 mm / 110 pxl = 2.1091 mm per pixel
% Anthropometry: vertical distance from the bottom of the chin (menton) to the top of the head: https://upload.wikimedia.org/wikipedia/commons/0/06/AvgHeadSizes.png
% US adult male, 50th percentile: Avg. = 23.2 cm, 9.1 inches
% Subject in video frames: Avg. = 110 pixels


% Loop through CSV files - Raw Data Processing
for csv_i = 1:length(moveCSV)

    tmpCSV = moveCSV{csv_i};

    % Split file names to extract relevant parts (dateID, sessID, and hemID)
    nameParts = split(tmpCSV,'_');
    dateID = nameParts{1};
    sessID = nameParts{3};
    hemID = nameParts{8};
    matName_title = [dateID , '-' , sessID, '-', hemID]

    % Find and load corresponding dlcDAT MAT file
    matTempfind = [dateID , '_' , sessID];
    matInd = contains(mainMAT2 , matTempfind);
    matName = mainMAT2{matInd};
    load(matName)

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
    fTipInds = contains(colNames4,'fTip');
    fTipEuclid = euclidall(:,fTipInds);

    % Average the computed distances related to fingertip movements
    fTipAverage = mean(fTipEuclid,2);

    % Process dlcDAT MAT files using MoveIndex CSV files to select specific portions of the averaged fingertip distances
    moveINDtab = readtable(tmpCSV);
    moveINDtab = moveINDtab(~moveINDtab.BeginF == 0,:); % clean up - filters out rows in moveINDtab where the BeginF field is zero.
    moveINDtab = moveINDtab(~moveINDtab.EndF == 0,:); % clean up - filters out rows in moveINDtab where the EndF field is zero.

    % Align with Euclidean distance frames
    firstBegin = moveINDtab.BeginF(1) - 1; % assigned to value of the first element in the BeginF column of moveINDtab - 1 (because MATLAB uses 1-based indexing)
    lastEnd = moveINDtab.EndF(height(moveINDtab)) - 1; % assigned to value of the last element in the EndF column of moveINDtab - 1

    % Extract and store the average fingertip distance for specified frames (in fTip Average Block)
    fTipAveBlk = fTipAverage(firstBegin:lastEnd); % extracts subset of fTipAverage w/in frame range from firstBegin to lastEnd (represents specific portion of data where specified movement is detected, as indicated in MoveIndex CSV file)

    % Smooth out edges -- smoothdata function w/ 'guassian' method
    window_Width = 5; % set windowWidth as needed
    smoothed_fTipAveBlk = smoothdata(fTipAveBlk, 'gaussian', window_Width); % read documentation, window overlap

    % set thresholds based on mean +/- 3xStd
    MinPeakHeight = mean(smoothed_fTipAveBlk);
    MinPeakProminence = mean(smoothed_fTipAveBlk)*2;

    % Rubric for movement type -OR- iterative algorthmic approach
    % MinPeakDistance_HandMov =
    % MinPeakDistance_PronSup =
    % MinPeakDistance_FlexExtend =
    % MinPeakDistance_FingerTap =

    % Find peak amplitudes and compute widths -- findpeaks function [review documentation]
    [peaks_tmp, locs_tmp, widths_tmp, prominences_tmp] = findpeaks(smoothed_fTipAveBlk, fps, MinPeakHeight=mean(smoothed_fTipAveBlk), MinPeakDistance=0.15, MinPeakProminence=mean(smoothed_fTipAveBlk)*2, Annotate ='extents');
    [peaks, locs, widths, prominences] = findpeaks(smoothed_fTipAveBlk, MinPeakHeight=mean(smoothed_fTipAveBlk), MinPeakDistance=15, MinPeakProminence=mean(smoothed_fTipAveBlk)*2, Annotate ='extents');

    % Convert distance variables to mm usng distance conversion factor
    amplitudes = peaks * pixels_to_mm; % converting amplitudes to mm

    % Compute timepoints from locs (vector of integer indices corresponding to video frame number)
    timepoints = locs / fps; % Convert frame numbers to time (in seconds) using video sampling rate (Fs) conversion factor

    % Compute distances between consecutive peaks
    peakDists_frames = diff(locs); % by frame indice
    peakDists_p = diff(timepoints); % by timepoint (in seconds)

    % Convert frame-relative variables to seconds using time conversion factor
    widths_fps = widths / fps; % converting widths to seconds
    halfWidths = widths_fps / 2;
    % slope

    % normalize within patient - relative to bl (Off, off)

    % Compute timepoints for all indices
    timepoints__fTipAveBlk = (1:length(smoothed_fTipAveBlk))/fps;

    % Plot smooth movement for each CSV iteration
    subplot(length(moveCSV), 1, csv_i);
    hold on
    % Adjust the parameters in findpeaks
    findpeaks(smoothed_fTipAveBlk, timepoints__fTipAveBlk, MinPeakHeight=mean(smoothed_fTipAveBlk), MinPeakDistance=0.15, MinPeakProminence=mean(smoothed_fTipAveBlk)*2, Annotate ='extents');
    % define axes labels and subplot titles
    xlabel('time (s)');
    ylabel('amplitude');
    hold off
    title(['Smooth Hand O/C Movement, ', num2str(matName_title)])

    % Define unique name for the findpeaks results based on the current CSV name
    findpeaks_output = [outputDir filesep 'findpeaks_output_' tmpCSV(1:end-44) '.csv']; % assumes tmpCSV is a string ending in '.csv'

    % Create a table to store results based on computed variables
    T1 = table(timepoints, locs, peaks, amplitudes, prominences, widths_fps, halfWidths, 'VariableNames', {'Timepoints', 'Locations', 'Peaks', 'Amplitudes', 'Prominences', 'Widths', 'HalfWidths'});

    % Write results table to a CSV file
    writetable(T1, findpeaks_output);

    % Define unique name for the outDATA_NaN output file based on the current CSV/MAT filename
    outDATA_NaN_filename = [outputDir filesep 'cleaned_dlcDAT_NaN_' tmpCSV(1:end-44) '.csv']; % Removes '.csv' from tmpCSV and appends 'outDATA_NaN_'
    % outDATA_NaN_filename = [outputDir filesep 'outDATA_NaN_' tmpCSV(1:end-44) '.csv']; % Removes '.csv' from tmpCSV and appends 'outDATA_NaN_'

    % Write outDATA_NaN table to a CSV file in the specified output directory
    writetable(outDATA_NaN, outDATA_NaN_filename);

end

cd(outputDir)
% save([outputDir filesep 'outDATA_NaN.csv']);


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



