% USE only wiht 2024a

% Load video
% Loop through frames
% Color code masks 
% Place number at centroid
% Explore masks

%%

pcNAME = getenv('COMPUTERNAME');

% Combine Percept LFP with DLC Video
% Determine frames for video based on recording start
switch pcNAME
    case 'DESKTOP-I5CPDO7' % work pc
        mainDIR = 'D:\Dropbox\erin_clinic_video\Clin_2023-09-12_session1_soffoffLi_Lcam-emr-2023-09-16\videos\SegTESTtwo';
        % Tablet video
    case'DESKTOP-FAGRV5G'
        mainDIR = 'E:\Dropbox\erin_clinic_video\Clin_2023-09-12_session1_soffoffLi_Lcam-emr-2023-09-16\videos';
end
cd(mainDIR)

%% Save video files at Matlab file - UNLABEL

unLABELvid_obj = VideoReader('20230912_idea08_session001_leftCam-0000c.mp4');
UNLabel_vid = struct('cdata',zeros(unLABELvid_obj.Height,unLABELvid_obj.Width,3,'uint8'),'colormap',[]);

frami = 1;
while hasFrame(unLABELvid_obj)
   UNLabel_vid(frami).cdata = readFrame(unLABELvid_obj);
   frami = frami+1;
end
disp('Video1 done!')

save('UNLABELvideo.mat','UNLabel_vid','-v7.3');

%% Save video file to Matlab file - LABEL

LABELvid_obj = VideoReader('20230912_idea08_session001_leftCam-0000DLC_resnet50_Clin_2023-09-12_session1_labeledc.mp4');
Label_vid = struct('cdata',zeros(LABELvid_obj.Height,LABELvid_obj.Width,3,'uint8'),'colormap',[]);

frami = 1;
while hasFrame(LABELvid_obj)
   Label_vid(frami).cdata = readFrame(LABELvid_obj);
   frami = frami+1;
end
disp('Video2 done!')

save('LABELvideo.mat','Label_vid','-v7.3');

%% Plot points from CSV over UNLABEL images to CHECK
% Part 1 - convert DLC csv to mat

csv_vidLoc = mainDIR;
save_matLoc = mainDIR;
user_CSV = 'MoveKinematicLabels_20230912_leftCam_session1.csv';


dlcIO_processCSV('dirLOC',1,'userLOC',csv_vidLoc,'saveLOC',save_matLoc,...
    'selCSV',user_CSV,'USERid','None')

%% Plot points from CSV over UNLABEL images to CHECK
% Part 2 - load and correct frames

load("dlcDAT_MoveKinematicLabels_20230912_leftCam_session1.mat")
tmpFRAMES = outDATA.labelTab.session1(2:end-1,:);

%% Plot points from CSV over UNLABEL images to CHECK
% Part 3 - Overlay and check
close all
tiledlayout(1,2,"TileSpacing","tight","Padding","tight")
nexttile()
imshow(Label_vid(100).cdata)
title('Ground Truth')

nexttile()
imshow(UNLabel_vid(100).cdata)
hold on
[x_VALUES , y_VALUES] = getCoordDLCframe(tmpFRAMES, 100);
scaTtest = scatter(x_VALUES,y_VALUES,30,"magenta","filled");
scaTtest.MarkerFaceAlpha = 0.5;
title('Label Test')


%% Create bounding around borders of points - CHECK

close all
tiledlayout(1,3,"TileSpacing","tight","Padding","tight")
nexttile()
imshow(Label_vid(100).cdata)
title('Ground Truth')

nexttile()
imshow(UNLabel_vid(100).cdata)
hold on
[x_VALUES , y_VALUES] = getCoordDLCframe(tmpFRAMES, 100);
scaTtest = scatter(x_VALUES,y_VALUES,30,"magenta","filled");
scaTtest.MarkerFaceAlpha = 0.5;
title('Label Test')

nexttile()
imshow(UNLabel_vid(100).cdata)
hold on
[x_VALUES , y_VALUES] = getCoordDLCframe(tmpFRAMES, 100);
scaTtest = scatter(x_VALUES,y_VALUES,30,"red","filled");
scaTtest.MarkerFaceAlpha = 0.5;

[xB_min,yB_min,B_width,B_height] = getBOUNDboxVs(x_VALUES,y_VALUES,0.2);

rectangle('Position', [xB_min, yB_min, B_width, B_height],...
    'EdgeColor', 'y', 'LineWidth', 1, 'Curvature',0.2);

% rectangle('Position', [x_MIN, y_MIN, width_box, height_box],...
%     'EdgeColor', 'y', 'LineWidth', 1, 'Curvature',0.2);

title('Bounding Box')


%% Check points in all frames

for uli = 1: size(UNLabel_vid,2)

    imshow(UNLabel_vid(uli).cdata)
    title('Ground Truth')
    hold on
    [x_VALUES , y_VALUES] = getCoordDLCframe(tmpFRAMES, uli);
    scaTtest = scatter(x_VALUES,y_VALUES,30,"red","filled");
    scaTtest.MarkerFaceAlpha = 0.5;
        title([num2str(uli)])

    pause(0.01)


end

%% Create SETPOINTS - check -- centroid and four points around

close all

centroidPOINT = mean([x_VALUES y_VALUES]);

% Define the percentage for expansion
expCEN_factor = 0.05;

% Calculate the offsets
cenX_offset = centroidPOINT(1) * expCEN_factor;
cenY_offset = centroidPOINT(2) * expCEN_factor;

% Calculate new coordinates
cenX_left =  centroidPOINT(1) - cenX_offset;
cenX_right = centroidPOINT(1) + cenX_offset;
cenY_below = centroidPOINT(2) - cenY_offset;
cenY_above = centroidPOINT(2) + cenY_offset;

allSETpoints = [centroidPOINT;...
                cenX_left , centroidPOINT(2);...
                cenX_right , centroidPOINT(2);...
                centroidPOINT(1), cenY_below;...
                centroidPOINT(1), cenY_above];

close all
tiledlayout(1,4,"TileSpacing","tight","Padding","tight")
nexttile()
imshow(Label_vid(100).cdata)
title('Ground Truth')

nexttile()
imshow(UNLabel_vid(100).cdata)
hold on
[x_VALUES , y_VALUES] = getCoordDLCframe(tmpFRAMES, 100);
scaTtest = scatter(x_VALUES,y_VALUES,30,"magenta","filled");
scaTtest.MarkerFaceAlpha = 0.5;
title('Label Test')

nexttile()
imshow(UNLabel_vid(100).cdata)
hold on
[x_VALUES , y_VALUES] = getCoordDLCframe(tmpFRAMES, 100);
scaTtest = scatter(x_VALUES,y_VALUES,30,"red","filled");
scaTtest.MarkerFaceAlpha = 0.5;

[xB_min,yB_min,B_width,B_height] = getBOUNDboxVs(x_VALUES,y_VALUES,0.2);

rectangle('Position', [xB_min, yB_min, B_width, B_height],...
    'EdgeColor', 'y', 'LineWidth', 1, 'Curvature',0.2);
title('Bounding Box')

nexttile()
imshow(UNLabel_vid(100).cdata)
hold on
scaTtest = scatter(allSETpoints(:,1),allSETpoints(:,2),30,"cyan","filled");
scaTtest.MarkerFaceAlpha = 0.5;

%% CHECK SUBSET
close all
for fii = 1:50:3000

    tmp2uFRAME = fii;

    imshow(UNLabel_vid(tmp2uFRAME).cdata)
    hold on
    [x_VALUES , y_VALUES] = getCoordDLCframe(tmpFRAMES, tmp2uFRAME);

    centroidPOINT = mean([x_VALUES y_VALUES]);

    % Define the percentage for expansion
    expCEN_factor = 0.05;

    % Calculate the offsets
    cenX_offset = centroidPOINT(1) * expCEN_factor;
    cenY_offset = centroidPOINT(2) * expCEN_factor;

    % Calculate new coordinates
    cenX_left =  centroidPOINT(1) - cenX_offset;
    cenX_right = centroidPOINT(1) + cenX_offset;
    cenY_below = centroidPOINT(2) - cenY_offset;
    cenY_above = centroidPOINT(2) + cenY_offset;

    allSETpoints = [centroidPOINT;...
        cenX_left , centroidPOINT(2);...
        cenX_right , centroidPOINT(2);...
        centroidPOINT(1), cenY_below;...
        centroidPOINT(1), cenY_above];

    scaTtest = scatter(allSETpoints(:,1),allSETpoints(:,2),30,"cyan","filled");
    scaTtest.MarkerFaceAlpha = 0.5;


    hold off
    pause

    cla


end

%% Add elbow and palm points to cell array
allPOINTS3 = cell(size(UNLabel_vid,2),1);

for ali = 1:size(UNLabel_vid,2)

    tmpROW = tmpFRAMES(ali,:);
    % Palm, MidForeArm , Elbow
    pemCols = colNAMES(contains(colNAMES,{'Palm','Elbow','MidForeArm'}));
    pemData = table2array(tmpROW(:,contains(colNAMES,{'Palm','Elbow','MidForeArm'})));
    %
    pemXvals = transpose(pemData(contains(pemCols,'x')));
    pemYvals = transpose(pemData(contains(pemCols,'y')));

    xyVALS = [pemXvals , pemYvals];

    allPOINTS3{ali} = xyVALS;
end


%% RUN SetPointDefine

%% Load in allpoints from GUI
load('PointFrameStart.mat')

nonNanInd = find(~isnan(allPoints(:,1)));
allPoints2 = allPoints;
for nii = 1:length(nonNanInd)

    tmpNi = nonNanInd(nii);
    if nii == length(nonNanInd)
        tmpNi2 = height(allPoints);
    else
        tmpNi2 = nonNanInd(nii+1)-1;
    end
    allPoints2(tmpNi:tmpNi2,:) = repmat(allPoints(tmpNi,:), length(tmpNi:tmpNi2),1);

end

allPoints2(:,3) = allPoints2(:,1) > 30;

%% Initialize Segment model

obj = segmentAnythingModel;


%%

UNlabel_vidMasks = UNLabel_vid;
for frameI = 1:size(UNLabel_vid,2)

    embeddings = extractEmbeddings(obj,UNLabel_vid(frameI).cdata);
    % promptPointROI = allPoints2(frameI,1:2);
    promptPointROI = allPOINTS3{frameI};

    [masks, ~] = segmentObjectsFromEmbeddings(obj, embeddings,...
        size(UNLabel_vid(frameI).cdata),...
        ForegroundPoints = promptPointROI,...
        ReturnMultiMask = true);
    UNlabel_vidMasks(frameI).cdata = masks;

    disp(frameI)

end



%% DEAL WITH MASKS ---- SegMasks

load("SegMasks.mat")
load("OrigVid.mat")
load("PointFrameRefine.mat")
%%
tab_vidARM = UNlabel_vidMasks;
for sii = 1:size(UNlabel_vidMasks,2)

    if isempty(UNlabel_vidMasks(sii).cdata)
        tab_vidARM(sii).cdata = [];
        continue
    else
        tmpIM = UNlabel_vidMasks(sii).cdata(:,:,2);
        tmpFRAME = UNLabel_vid(sii).cdata;
        tmpPoint = allPoints2(sii,:);

        % BW = imbinarize(tmpIM);

        [B,L] = bwboundaries(tmpIM,'noholes');
        LL = logical(L);

        % [B2 , L2 , L2M] = getRidOfsmall(B,L);

        stats = regionprops("table",LL,"Area","Centroid","PixelIdxList", ...
            "MajorAxisLength","MinorAxisLength","BoundingBox");

        LLBW = struct;
        LLBW.ImageSize = size(L);
        LLBW.PixelIdxList = stats.PixelIdxList;
        LLBW.Connectivity = 8;
        LLBW.NumObjects = length(LLBW.PixelIdxList);

        areaS = stats.Area;
        cenS = stats.Centroid;
        selection = (areaS > 5000) & (cenS(:,1) < 450);
        BW2 = cc2bw(LLBW, ObjectsToKeep = selection);
        [B3,L3] = bwboundaries(BW2,'noholes');
        L3mask = logical(L3);

        % CC = bwpropfilt(LL,"Area",[5000 Inf]);

        stats2 = regionprops("table",BW2,"Area","Centroid", ...
            "MajorAxisLength","MinorAxisLength","BoundingBox");

        % % Image (mask) dimensions
        % imageHeight = size(BW2,1);
        % imageWidth = size(BW2,2);
        % 
        % % Rectangle coordinates [x, y, width, height]
        % rect = stats2.BoundingBox;
        % 
        % % Initialize the mask
        % mask = false(imageHeight, imageWidth);
        % 
        % % Fill in the rectangle
        % xStart = round(rect(1));
        % yStart = round(rect(2));
        % xEnd = xStart + round(rect(3)) - 1;
        % yEnd = yStart + round(rect(4)) - 1;
        % mask(yStart:yEnd, xStart:xEnd) = true;
        % 
        % BW8 = activecontour(tmpFRAME,mask,100,'edge');
        % imshow(BW8)

        % grayIM = rgb2gray(tmpFRAME);
        % points = detectFASTFeatures(grayIM,'ROI', stats2.BoundingBox);

        if isempty(stats2)
            continue
        else

            % USE mask
            tmpARM = tmpFRAME;

            invertMASK = ~L3mask;
            invertedMask3D = repmat(invertMASK, [1, 1, 3]);
            tmpARM(invertedMask3D) = 255;

            imshow(tmpARM)
            title(num2str(sii))

            tab_vidARM(sii).cdata = tmpARM;
            tab_vidARM(sii).B3 = B3;
            tab_vidARM(sii).L3 = L3;
            tab_vidARM(sii).Rprops = stats2;

            % imshow(label2rgb(L3, @jet, [.5 .5 .5]))
            hold on

            grayIM = rgb2gray(tmpARM);
            grayIM2 = rgb2gray(tmpFRAME);

            % ----------PLOT BOUNDARY AROUND ARM FOR EMPHASIS
            % for k = 1:length(B3)
            %     boundary = B3{k};
            %     plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
            % end

            %--------------- PLOT RECTANGLE TO ADD INFO 
            % rectangle('Position',stats2.BoundingBox,'Curvature',0.2)

            % plot(pointsEF.selectStrongest(5))
            % plot(pointsFAST.selectStrongest(5))

        end

        % TO DO FUTURE
        % USE BOUNDING BOX - USE TO REFINE OUTLINE






    end
    pause(0.01)
    cla

end

test = 1;


%%

I = imread("rice.png");
BW = imbinarize(I);
imshow(BW)
CC = bwconncomp(BW); 
stats = regionprops("table",CC,"Area","BoundingBox");
area = stats.Area;
bbox = stats.BoundingBox;
selection = (area > 50) & (bbox(:,3) < 15) & (bbox(:,4) >= 20);
%%













% function [B2 , L2] = getRidOfsmall(B,L)
% 
% threshS = 60;
% keepIDlog = cellfun(@(x) round(numel(x)/2), B, "UniformOutput",true);
% B2 = B(keepIDlog > threshS);
% ridIDnum = find(keepIDlog < threshS);
% 
% L2 = L;
% L2(ismember(L,ridIDnum)) = 0;
% 
% end