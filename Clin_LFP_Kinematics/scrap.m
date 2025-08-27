% Old method
%% Determine duration (in seconds) of LEFT STN streaming sessions in OFF Med JSON Session Report

L_sessTimes_1_cell = stream_LEFT_1.FirstPacketDateTime; % cell
L_sessTimes_1_table =  cell2table(L_sessTimes_1_cell); % table
L_sessTimes_1 = table2array(L_sessTimes_1_table); % array

% Initialize array to store trimmed session times (only min, sec, millisec)
trimmed_L_sessTimes_1 = strings(size(L_sessTimes_1));

% Loop through each time string
for L_time_i = 1:length(L_sessTimes_1)
    % Convert string to datetime
    dt = datetime(L_sessTimes_1(L_time_i), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z');

    % Format to keep only minute, second, and millisecond
    trimmed_L_sessTimes_1(L_time_i) = datestr(dt, 'MM:SS.FFF');
end

% Convert the trimmed times to duration
trimmed_L_sessDurs_1 = duration(trimmed_L_sessTimes_1, 'InputFormat', 'mm:ss.SSS');

% Find the offset from the first time
L_sessTimeOffsets_1 = trimmed_L_sessDurs_1 - trimmed_L_sessDurs_1(1);

% Initialize a base time
baseTime = duration(0, 0, 0);

% Apply the offset to the base time
L_uniformTimes_1 = baseTime + L_sessTimeOffsets_1;

% Calculate differences between consecutive times
L_sessDurations_1 = diff(L_uniformTimes_1);

% Convert the duration array to seconds
L_sessDurations_1_seconds = seconds(L_sessDurations_1);

%% Determine duration (in seconds) of RIGHT STN streaming sessions in OFF Med JSON Session Report

R_sessTimes_1_cell = stream_RIGHT_1.FirstPacketDateTime; % cell
R_sessTimes_1_table =  cell2table(R_sessTimes_1_cell); % table
R_sessTimes_1 = table2array(R_sessTimes_1_table); % array

% Initialize array to store trimmed session times (only min, sec, millisec)
trimmed_R_sessTimes_1 = strings(size(R_sessTimes_1));

% Loop through each time string
for R_time_i = 1:length(R_sessTimes_1)
    % Convert string to datetime
    dt = datetime(R_sessTimes_1(R_time_i), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z');

    % Format to keep only minute, second, and millisecond
    trimmed_R_sessTimes_1(R_time_i) = datestr(dt, 'MM:SS.FFF');
end

% Convert the trimmed times to duration
trimmed_R_sessDurs_1 = duration(trimmed_R_sessTimes_1, 'InputFormat', 'mm:ss.SSS');

% Find the offset from the first time
R_sessTimeOffsets_1 = trimmed_R_sessDurs_1 - trimmed_R_sessDurs_1(1);

% Initialize a base time
baseTime = duration(0, 0, 0);

% Apply the offset to the base time
R_uniformTimes_1 = baseTime + R_sessTimeOffsets_1;

% Calculate differences between consecutive times
R_sessDurations_1 = diff(R_uniformTimes_1);

% Convert the duration array to seconds
R_sessDurations_1_seconds = seconds(R_sessDurations_1);



%%

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
outputDir = [mainDir2 filesep 'tempTests'];
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
    hemID = nameParts{8};
    matName_title = [dateID , '-' , sessID, '-', hemID];

    % Find and load corresponding dlcDAT MAT file
    matTempfind = [dateID , '_' , sessID];
    matInd = contains(mainMAT2 , matTempfind);
    matName = mainMAT2{matInd};
    load(matName)

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
    % euclidall = euclidall * pixels_to_mm; % converting euclidean distances to mm

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

end

% cd(outputDir)


%% 

%%

tab_vidObj = VideoReader('20230912_idea08_session001_rightCam-0000DLC_resnet50_Clin_2023-09-12_LSTN_v2Oct3shuffle1_100000_labeled.mp4');
tab_vid = struct('cdata',zeros(tab_vidObj.Height,tab_vidObj.Width,3,'uint8'),'colormap',[]);

%% Convert Tablet Video to dataframe

frami = 1;
while hasFrame(tab_vidObj)
   tab_vid(frami).cdata = readFrame(tab_vidObj);
%    imshow(frame)
   frami = frami+1;
end
whos tab_vid
v1who = whos('tab_vid');
round(v1who.bytes/1000000/1000,2)
disp('Video1 done!')

%% Determine which frames to trim on based presence of tablet
% Loop through / Plot / Title with Frame #

tabTab_vid = struct2table(tab_vid);

% for fi = 2900:height(tabTab_vid)
%     imshow(tabTab_vid.cdata{fi})
%     title(num2str(fi))
%     pause
% end

% Start 196
% Stop 3059

%% Load in experimental camera

% dlc_vidObj = VideoReader('20230912_idea08_session001_leftCam-0000.mp4');
% dlc_vid = struct('cdata',zeros(dlc_vidObj.Height,dlc_vidObj.Width,3,'uint8'),'colormap',[]);
% 
% frami = 1;
% while hasFrame(dlc_vidObj)
%    dlc_vid(frami).cdata = readFrame(dlc_vidObj);
% %    imshow(frame)
%    frami = frami+1;
% end
% 
% dlcTab_vid = struct2table(dlc_vid);
% disp('Video1 done!')

%% Trim DLC video by start and stop from Tablet video

dlcTab_vid2use = tabTab_vid(206:3010,:);

%% 

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
outputDir = [mainDir2 filesep 'tempTests'];
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
    hemID = nameParts{8};
    matName_title = [dateID , '-' , sessID, '-', hemID];

    % Find and load corresponding dlcDAT MAT file
    matTempfind = [dateID , '_' , sessID];
    matInd = contains(mainMAT2 , matTempfind);
    matName = mainMAT2{matInd};
    load(matName)

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
    % euclidall = euclidall * pixels_to_mm; % converting euclidean distances to mm

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

end


%% Load and Trim MAT file

% load('dlcDAT_20230912_session001_idea08_resnet50_rightCam-0000DLC.mat','outDATA')
% d1_labDAT = outDATA.labelTab.leftCam;

load('dlcDAT_20230912_session001_idea08_resnet50_rightCam-0000DLC.mat')
d1_labDAT = outDATA.Properties.VariableNames;



%%

dlc_lab2use = d1_labDAT(206:3010,:);

% offsetSamples = 195;
% offsetSecs = 195/60;
totalNumSamplesV = height(dlc_lab2use);
totalNumSecs = totalNumSamplesV/60; % 60 fps

totalNumSampsLFP = floor(totalNumSecs*250);

%%
% Original signal at 60 Hz
ts_DLC = 0:1/60:(height(dlc_lab2use)-1)/60;
% Target sampling rate at 250 Hz
ts_LFP = 0:1/250:(height(streamOfInt)-1)/250;

allColNs = dlc_lab2use.Properties.VariableNames;
dlc_lab2use2int = table;
for coli = 1:width(dlc_lab2use)

    % Tmp col
    tmpCol = allColNs{coli};
    
    % Interpolation
    x_250 = interp1(ts_DLC, transpose(dlc_lab2use.(tmpCol)),...
        ts_LFP, 'spline');

    dlc_lab2use2int.(tmpCol) = transpose(x_250);
end


%% Kinematic and LFP plot
close all

tiledlayout(2,1,"TileSpacing","tight")
xTime = ts_LFP;
% plot(dlc_lab2use2int.Elbow_x)
% hold on
% yyaxis left
elbowXsm = smoothdata(dlc_lab2use2int.Elbow_x,'gaussian',70);
nexttile
plot(xTime,elbowXsm);
xlim([0 round(max(ts_LFP))])
ylabel('X deflection')
xlabel('Time in seconds')
ecg = perceive_ecg(transpose(streamOfInt),250,0);
% plot(ecg.cleandata);
% yyaxis right
nexttile
plot(xTime,ecg.cleandata);
xlim([0 round(max(ts_LFP))])
ylabel('uV')
xlabel('Time in seconds')



%% 
close all
for di = 1:height(dlcTab_vid)
     imshow(dlcTab_vid.cdata{di})
     pause

end


%%

for ffi = 1:height(dlcTab_vid2use.cdata)
    imshow(dlcTab_vid2use.cdata{ffi})
    hold on
    plot(dlc_lab2use.PalmBase_x(ffi),dlc_lab2use.PalmBase_y(ffi),'ro')
    pause(0.01)
    cla
end

%% To do

% 1. clean up interpolated point of interest

% 2. clean up LFP with artifact removal

% 3. instanteous theta , beta and gamma


%% Check DLC labeled video
dlcLab_vidLoc = mainDIR;
cd(dlcLab_vidLoc)

dlc_lab_vidObj = VideoReader('20230912_idea08_session001_leftCam-0000DLC_resnet50_Clin_2023-09-12_session1_labeled.mp4');
dlc_lab_vid = struct('cdata',zeros(dlc_lab_vidObj.Height,dlc_lab_vidObj.Width,3,'uint8'),'colormap',[]);

frami = 1;
while hasFrame(dlc_lab_vidObj)
   dlc_lab_vid(frami).cdata = readFrame(dlc_lab_vidObj);
%    imshow(frame)
   frami = frami+1;
end
dlcLabTab_vid = struct2table(dlc_lab_vid);
disp('Video1 done!')


dlc_lablab2use = dlcLabTab_vid(206:3010,:);


%% Animated plot

plotXaxis = xTime;
kinematicY = elbowXsm;
lfpY = ecg.cleandata;

% dlc_lablab2use


% 4.1 samples per frame
stePS = round(linspace(1,11688,2805));

for fi = 1:height(dlc_lablab2use)
    subplot(4,1,1:2)

    imshow(dlc_lablab2use.cdata{fi})

    subplot(4,1,3)
    plot(plotXaxis(1:stePS(fi)),kinematicY(1:stePS(fi)))
    ylim([min(kinematicY) max(kinematicY)])
    xlim([0 max(plotXaxis)])

    subplot(4,1,4)
    plot(plotXaxis(1:stePS(fi)),lfpY(1:stePS(fi)))
    ylim([min(lfpY) max(lfpY)])
    xlim([0 max(plotXaxis)])

    if fi == 1
        pause
    end

    pause(0.01)

end


%%
for ffi = 1:height(dlc_lablab2use.cdata)
    imshow(dlc_lablab2use.cdata{ffi})
    hold on
    % plot(dlc_lab2use.PalmBase_x(ffi),dlc_lab2use.PalmBase_y(ffi),'ro')
    pause(0.01)
    cla
end
