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
    case 'jst'
        mainDIR = 'D:\Dropbox\erin_clinic_video\Clin_2023-09-12_session1_soffoffLi_Lcam-emr-2023-09-16\videos';
        % Tablet video
    case'DESKTOP-FAGRV5G'

        mainDIR = 'E:\Dropbox\erin_clinic_video\Clin_2023-09-12_session1_soffoffLi_Lcam-emr-2023-09-16\videos';

end
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