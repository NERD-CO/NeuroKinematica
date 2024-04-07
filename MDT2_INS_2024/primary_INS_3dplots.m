%%
% generateINS_DatasetS_v2('MDT8' , 'left')

cd('E:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data\MDT7\RightBody\FOR_INS_3D_plot')
load("MDT8_GroupB_Test2_HAND_PS.mat")
% finger tap
% videoID = 't1_20230816_idea08_session002_rightCam-0000DLC_resnet50_INS_2024_MPR8_LSTNMar20shuffle1_100000_labeled.mp4';
% Hand OC
% videoID = 't1_20230816_idea08_session005_rightCam-0000DLC_resnet50_INS_2024_MPR7_LSTNMar20shuffle1_100000_labeled.mp4'
% Hand PS
videoID = 't2_20230816_idea08_session003_rightCam-0000DLC_resnet50_INS_2024_MPR8_LSTNMar20shuffle1_100000_labeled.mp4';
%% HAND OPEN CLOSE
% PCA plot of movement kin
% REMOVE LIKELIHOOD and PUSH THROUGH
colNames = dlc_lab2use2int.Properties.VariableNames; %
colNames2 = cellfun(@(x) split(x,'_'), colNames,...
    'UniformOutput',false);
colNames3 = unique(cellfun(@(x) x{1}, colNames2,...
    'UniformOutput',false));
colNames4 = colNames3(~matches(colNames3,'frames'));

tmpMOve.MoveType{1}; %%%% GROUP A ---- HAND OC

% startINDi2 = startINDi + 50;
% endINDi2 = endINDi - 50;
trimFRAMES = dlc_lab2use2int(startINDi2:endINDi2,:);

% REMOVE likelihood
colKEEP = contains(trimFRAMES.Properties.VariableNames,{'_x','_y'});
trimFRAMES2 = trimFRAMES(:,colKEEP);
[~,pcSCORES] = pca(table2array(trimFRAMES2));

top3 = pcSCORES(:,1:3);
% figure;
% plot3(top3(:,1),top3(:,2),top3(:,3))
figure;
plot(top3(:,1))

% figure('NumberTitle', 'off', 'Name', 'test')
%%%%%% THIS IS THE STATIC VERSION OF 3D TIMEPLOT
% % remove none fingers;
% getVnames0 = trimFRAMES.Properties.VariableNames;
% getVnames = getVnames0(~contains(getVnames0,'frames'));
% getVnames2 = getVnames(~contains(getVnames,{'MidForeArm','Elbow'}));

% frameTable = dataTABLE.frames;
% pointTable = trimFRAMES(:,getVnames2);

% cVnames0 = cellfun(@(x) strsplit(x,'_'), getVnames2, 'UniformOutput', false);
% cVnames = cellfun(@(x) x{1}, cVnames0, 'UniformOutput', false);
% uniNames = unique(cVnames);

% plasCMP = colormap(plasma);
% plasCMP2u = plasCMP(round(linspace(1,256,length(uniNames))),:);
% for fi = 1:length(uniNames)
%
%     tUiName = uniNames{fi};
%     uiNindex = ismember(cVnames,tUiName);
%     uiTab = pointTable(:,uiNindex);
%
%     % Get X and Y
%     tPnames0 = cellfun(@(x) strsplit(x,'_'),...
%         uiTab.Properties.VariableNames, 'UniformOutput', false);
%     tPnames = cellfun(@(x) x{2}, tPnames0, 'UniformOutput', false);
%     xyIND = contains(tPnames,{'x','y'});
%     xyiTab = table2array(uiTab(:,xyIND));
%
%     allFrames = xyiTab;
%     allColors = repmat(plasCMP2u(fi,:),height(uiTab) ,1);
%
%     % scatter(allFrames(:,1),allFrames(:,2),30,allColors,'filled')
%     scatter3(transpose(1:length(allFrames)),allFrames(:,1),allFrames(:,2),30,allColors,'filled')
%     hold on
%
% end
% % set view
% view([41 34])

%% LOAD IN CORRECT VIDEO

% plot the video
% t1_20230816_idea08_session005_rightCam-0000DLC_resnet50_INS_2024_MPR7_LSTNMar20shuffle1_100000_labeled.mp4
% dlcLab_vidLoc = vidLoc;
% cd(dlcLab_vidLoc)

dlc_lab_vidObj = VideoReader(videoID);
dlc_lab_vid = struct('cdata',zeros(dlc_lab_vidObj.Height,dlc_lab_vidObj.Width,3,'uint8'),'colormap',[]);

frami = 1;
while hasFrame(dlc_lab_vidObj)
    dlc_lab_vid(frami).cdata = readFrame(dlc_lab_vidObj);
    %    imshow(frame)
    frami = frami+1;
end
dlcLabTab_vid = struct2table(dlc_lab_vid);
disp('Video1 done!')


dlc_lablab2useVIDE = dlcLabTab_vid(startIND:endIND,:);

%%

figure;
getVnames0 = trimFRAMES.Properties.VariableNames;
getVnames = getVnames0(~contains(getVnames0,'frames'));
getVnames2 = getVnames(~contains(getVnames,{'MidForeArm','Elbow'}));

% frameTable = dataTABLE.frames;
pointTable = trimFRAMES(:,getVnames2);

cVnames0 = cellfun(@(x) strsplit(x,'_'), getVnames2, 'UniformOutput', false);
cVnames = cellfun(@(x) x{1}, cVnames0, 'UniformOutput', false);
uniNames = unique(cVnames);

plasCMP = colormap(plasma);
plasCMP2u = plasCMP(round(linspace(1,256,length(uniNames))),:);

% numFRAMES = floor(1695/433);
startINDS = transpose(round(linspace(1,height(pointTable),height(dlc_lablab2useVIDE))));
stopINDS = startINDS(2:end)-1;
startINDS = startINDS(1:end-1);
allXvaLS = table2array(pointTable(:,contains(pointTable.Properties.VariableNames,{'_x'})));
allYvaLS = table2array(pointTable(:,contains(pointTable.Properties.VariableNames,{'_y'})));
figure
set(gcf,'Position',[1000 124 706 1114])
for vii = 1:height(dlc_lablab2useVIDE)-1
    tmpFrame = dlc_lablab2useVIDE.cdata{vii,1};
    subplot(2,1,1)
    imshow(tmpFrame)

    if vii == 1
        pause
    end


    % 4 points
    subplot(2,1,2)
    frameIIx = allXvaLS(startINDS(1):stopINDS(vii),:);
    frameIIx2 = reshape(frameIIx,numel(frameIIx),1);
    frameIIy = allYvaLS(startINDS(1):stopINDS(vii),:);
    frameIIy2 = reshape(frameIIy,numel(frameIIy),1);
    frameIIz = repmat(transpose(1:stopINDS(vii)),1,11);
    frameIIz2 = reshape(frameIIz,numel(frameIIz),1);


    rowIndices = repmat(1:size(plasCMP2u,1), stopINDS(vii), 1);
    rowIndices = rowIndices(:);
    plasCMP2u2u = plasCMP2u(rowIndices, :);


    scatter3(frameIIz2,frameIIx2,frameIIy2,60,plasCMP2u2u,'filled')
    view([41 34])


    pause(0.01)

end

%%







% 3D plot of movement kin

% cool animation plot

% wavelet cluster analysis ---- see pubicaiton
















function [xDims , yDims] = getDataDims(pointTable)

tPnames0 = cellfun(@(x) strsplit(x,'_'),...
    pointTable.Properties.VariableNames, 'UniformOutput', false);
tPnames = cellfun(@(x) x{2}, tPnames0, 'UniformOutput', false);
xyIND = contains(tPnames,{'x','y'});

xyiTab = table2array(pointTable(:,xyIND));

frameTx = xyiTab(:,1:2:size(xyiTab,2));
frameTy = xyiTab(:,2:2:size(xyiTab,2));

xDims.min = min(frameTx,[],'all');
xDims.max = max(frameTx,[],'all');

yDims.min = min(frameTy,[],'all');
yDims.max = max(frameTy,[],'all');


end