% USE only wiht 2024a

% Load video
% Loop through frames
% Color code masks 
% Place number at centroid
% Explore masks

%%

<<<<<<< Updated upstream
pcNAME = getenv('COMPUTERNAME');


% Combine Percept LFP with DLC Video

% Determine frames for video based on recording start

switch pcNAME
    case 'jst'
        mainDIR = 'D:\Dropbox\erin_clinic_video\Clin_2023-09-12_session1_soffoffLi_Lcam-emr-2023-09-16\videos';
        % Tablet video
    case'DESKTOP-FAGRV5G'

        mainDIR = 'E:\Dropbox\erin_clinic_video\Clin_2023-09-12_session1_soffoffLi_Lcam-emr-2023-09-16\videos';

end
=======
%% Combine Percept LFP with DLC Video

% Determine frames for video based on recording start

mainDIR = 'D:\Dropbox\erin_clinic_video\Clin_2023-09-12_session1_soffoffLi_Lcam-emr-2023-09-16\videos';
% Tablet video

>>>>>>> Stashed changes
tab_vidLoc = [mainDIR, filesep , 'TabletVideo'];
cd(tab_vidLoc)

tab_vidObj = VideoReader('20230912_idea08_session001_rightCam-0000.mp4');
tab_vid = struct('cdata',zeros(tab_vidObj.Height,tab_vidObj.Width,3,'uint8'),'colormap',[]);


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

<<<<<<< Updated upstream
save('OrigVid.mat','tab_vid','-v7.3');

=======
<<<<<<< Updated upstream
>>>>>>> Stashed changes
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

tab_vidMasks = tab_vid;
for frameI = 1:size(tab_vid,2)

    if allPoints2(frameI,3) == 1

        embeddings = extractEmbeddings(obj,tab_vid(frameI).cdata);

        promptPointROI = allPoints2(frameI,1:2);

        [masks, ~] = segmentObjectsFromEmbeddings(obj, embeddings,...
            size(tab_vid(frameI).cdata),...
            ForegroundPoints = promptPointROI,...
            ReturnMultiMask = true);
        tab_vidMasks(frameI).cdata = masks;
    else
        tab_vidMasks(frameI).cdata = [];
        continue
    end
    disp(frameI)

end


=======
>>>>>>> Stashed changes
%% 

imshow(tab_vid(300).cdata)
roi = drawpoint;
promptPoint = roi.Position;

embeddings = extractEmbeddings(obj,tab_vid(300).cdata);
[masks, scores] = segmentObjectsFromEmbeddings(obj, embeddings, size(tab_vid(300).cdata),ForegroundPoints=promptPoint, ReturnMultiMask=true)

%%
ii = tab_vid(300).cdata;
showImage = labeloverlay(ii,masks(:,:,3));
imshow(showImage)
% imagesc(masks) % 3 matrix mas

%%

obj = segmentAnythingModel;
I = imread('TEST_FASTsamImage.png');
imshow(I);
roi = drawpoint;

promptPoint = roi.Position;
embeddings = extractEmbeddings(obj,I);

[masks, scores] = segmentObjectsFromEmbeddings(obj, embeddings, size(I), ForegroundPoints=promptPoint);

showImage = labeloverlay(I,masks);

imshow(showImage)

title("SAM(predict score:"+string(scores)+")");

%%
I2 = rgb2gray(I);
pout_imadjust = imadjust(I2);
pout_histeq = histeq(I2);

%%
embeddings = extractEmbeddings(obj,pout_histeq);
[masks2, scores] = segmentObjectsFromEmbeddings(obj, embeddings, size(pout_histeq), ForegroundPoints=promptPoint);
%%
showImage = labeloverlay(pout_histeq,masks2);

imshow(showImage)

%%

imshow(masks)


%% DEAL WITH MASKS ---- SegMasks

load("SegMasks.mat")
load("OrigVid.mat")
load("PointFrameRefine2.mat")
%%
tab_vidARM = tab_vidMasks;
for sii = 1:size(tab_vidMasks,2)

    if isempty(tab_vidMasks(sii).cdata)
        tab_vidARM(sii).cdata = [];
        continue
    else
        tmpIM = tab_vidMasks(sii).cdata(:,:,2);
        tmpFRAME = tab_vid(sii).cdata;
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

            tab_vidARM(sii).cdata = tmpARM;
            tab_vidARM(sii).B3 = B3;
            tab_vidARM(sii).L3 = L3;
            tab_vidARM(sii).Rprops = stats2;

            % imshow(label2rgb(L3, @jet, [.5 .5 .5]))
            hold on

            grayIM = rgb2gray(tmpARM);
            grayIM2 = rgb2gray(tmpFRAME);
            pointsEF = detectMinEigenFeatures(grayIM,'ROI', stats2.BoundingBox);
            pointsFAST = detectFASTFeatures(grayIM,'ROI', stats2.BoundingBox);
            [regions,mserCC] = detectMSERFeatures(grayIM2,'ROI', stats2.BoundingBox,...
                'ThresholdDelta',8);
            % ----------PLOT BOUNDARY AROUND ARM FOR EMPHASIS
            % for k = 1:length(B3)
            %     boundary = B3{k};
            %     plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
            % end

            %--------------- PLOT RECTANGLE TO ADD INFO 
            % rectangle('Position',stats2.BoundingBox,'Curvature',0.2)

            plot(regions,'showPixelList',true,'showEllipses',false)
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