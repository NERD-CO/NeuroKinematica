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
    case 'DESKTOP-I5CPDO7'
        mainDIR = 'D:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data\SegmentTEST\TabletVideo';
        % Tablet video
    case'DESKTOP-FAGRV5G'
        mainDIR = 'E:\Dropbox\erin_clinic_video\Clin_2023-09-12_session1_soffoffLi_Lcam-emr-2023-09-16\videos';

end
cd(mainDIR)

%% DEAL WITH MASKS ---- SegMasks

load("SegMasks.mat")
load("OrigVid.mat")
load("PointFrameRefine.mat")
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

        % [B2 , L2 , L2M] = getRidOfsmall_NN(B,L);

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

            % ----------PLOT BOUNDARY AROUND ARM FOR EMPHASIS
            % for k = 1:length(B3)
            %     boundary = B3{k};
            %     plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
            % end

            %--------------- PLOT RECTANGLE TO ADD INFO 
            rectangle('Position',stats2.BoundingBox,'Curvature',0.2)

            % plot(pointsEF.selectStrongest(5))
            % plot(pointsFAST.selectStrongest(5))

        end

        % TO DO FUTURE
        % USE BOUNDING BOX - USE TO REFINE OUTLINE






    end
    pause(0.01)
    % if mod(sii,20) == 0
    %     cd('D:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data\SegmentTEST\armImages')
    %     frameNAME = ['frameN_',num2str(sii),'.png'];
    %     exportgraphics(gca,frameNAME,'BackgroundColor','none')
    % end
    cla

end

test = 1;


%%
