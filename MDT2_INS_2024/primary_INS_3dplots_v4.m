%%
% generateINS_DatasetS_v2('MDT8' , 'left')
close all
clear
cd('E:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data\MDT7\RightBody\FOR_INS_3D_plot')
movePICK = 1;
colorSWf = 4;
switch movePICK
    case 1 % HAND OC ---- MOVIE
        load("MDT7_GroupA_Test1_HAND_OC.mat")
        videoID = 't1_20230816_idea08_session005_rightCam-0000DLC_resnet50_INS_2024_MPR7_LSTNMar20shuffle1_100000_labeled.mp4';
        videoDLC = 'dlcDAT_t1_idea08_session005_LSTN_RBODY.mat';
        cropRECT = [84.5100   38.5100  162.9800  340.9800];
        view2use = [100.3367 83.3221];
    case 2 % HAND PS
        load("MDT8_GroupB_Test2_HAND_PS.mat")
        videoID = 't2_20230816_idea08_session003_rightCam-0000DLC_resnet50_INS_2024_MPR8_LSTNMar20shuffle1_100000_labeled.mp4';
        useALPHA = 0.8;
    case 3 % Finger
        load("MDT8_GroupA_Test1_Finger.mat");
        videoID = 't1_20230816_idea08_session002_rightCam-0000DLC_resnet50_INS_2024_MPR8_LSTNMar20shuffle1_100000_labeled.mp4';
        useALPHA = 0.2;
end



%%
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

load(videoDLC,"outDATA")
dlcVIDm = outDATA.labelTab.session005;
orgtrimFr = dlcVIDm(startIND:endIND,:);


% REMOVE likelihood
% colKEEP = contains(trimFRAMES.Properties.VariableNames,{'_x','_y'});
% trimFRAMES2 = trimFRAMES(:,colKEEP);
% [~,pcSCORES] = pca(table2array(trimFRAMES2));

trimFRAMES = orgtrimFr
colKEEP = contains(orgtrimFr.Properties.VariableNames,{'_x','_y'});
trimFRAMES2 = orgtrimFr(:,colKEEP);

% top3 = pcSCORES(:,1:3);
% % figure;
% % plot3(top3(:,1),top3(:,2),top3(:,3))
% figure;
% plot(top3(:,1))




%%

% figure;
getVnames0 = trimFRAMES.Properties.VariableNames;
getVnames = getVnames0(~contains(getVnames0,'frames'));
getVnames2 = getVnames(~contains(getVnames,{'MidForeArm','Elbow'}));

% frameTable = dataTABLE.frames;
pointTable = trimFRAMES(:,getVnames2);

cVnames0 = cellfun(@(x) strsplit(x,'_'), getVnames2, 'UniformOutput', false);
cVnames = cellfun(@(x) x{1}, cVnames0, 'UniformOutput', false);
uniNames = unique(cVnames);

plasCMP = colormap(plasma);
plasCMPfl = flipud(plasCMP);
close all

switch colorSWf
    case 1
        plasCMP2u = plasCMP; % normal
    case 2
        plasCMP2u = plasCMPfl; % flipped;
    case 3 % ATTEMPT TO TRIM flipped
        plasCMP2u = plasCMPfl(56:256,:);
    case 4 % ATTEMPT TO TRIM normal
        plasCMP2u = plasCMP(1:200,:);
        % plasCMP2u = flipud(plasCMP2u);
end

plasCMP2u = plasCMP2u(round(linspace(1,height(plasCMP2u),length(uniNames))),:);

% numFRAMES = floor(1695/433);
% startINDS = transpose(round(linspace(1,height(pointTable),height(dlc_lablab2useVIDE))));
% startINDS = round(linspace(1,height(pointTable),1695/5));
% stopINDS = startINDS(2:end)-1;
% startINDS = startINDS(1:end-1);
allXvaLS = table2array(pointTable(:,contains(pointTable.Properties.VariableNames,{'_x'})));
allYvaLS = table2array(pointTable(:,contains(pointTable.Properties.VariableNames,{'_y'})));

% flip
allYvaLS = -allYvaLS;
allXvaLS = -allXvaLS;

if movePICK == 1

    figure
    set(gcf,'Position',[1000 124 706 1114])
    for vii = 1:height(dlc_lablab2useVIDE)-1
        tmpFrame = dlc_lablab2useVIDE.cdata{vii,1};
        subplot(2,1,1)
        % imshow(tmpFrame)
        tmpFrameTrim = imcrop(tmpFrame,cropRECT);
        imshow(tmpFrameTrim)
        axis fill

        if vii == 2
            pause
        end

        % 4 points
        subplot(2,1,2)
        % [xDims , yDims] = getDataDims(pointTable);
        [xDims , yDims] = getDataDims2(allXvaLS , allYvaLS);

        frameIIx = allXvaLS(1:vii,:);
        frameIIx2 = reshape(frameIIx,numel(frameIIx),1);
        frameIIy = allYvaLS(1:vii,:);
        frameIIy2 = reshape(frameIIy,numel(frameIIy),1);
        frameIIz = repmat(transpose(1:vii),1,11);
        frameIIz2 = reshape(frameIIz,numel(frameIIz),1);

        rowIndices = repmat(1:size(plasCMP2u,1), vii, 1);
        rowIndices = rowIndices(:);
        plasCMP2u2u = plasCMP2u(rowIndices, :);
        % X-Axis is frame number or frameIIz2
        % Y-Axis is X-value or frameIIx2
        
        ssALL = scatter3(frameIIz2,frameIIx2,frameIIy2,60,plasCMP2u2u,'filled');
        ssALL.MarkerFaceAlpha = 0.6;
        view(view2use)

        xticklabels([])
        yticklabels([])
        zticklabels([])
        ylabel('X pixel location')
        xlabel('Frame number')
        zlabel('Y pixel location')
        xlim([1 height(pointTable)])
        zlim([floor(yDims.min) ceil(yDims.max)])
        ylim([floor(xDims.min) ceil(xDims.max)])

        box off

        pause(0.001)

    end

end


if movePICK ~= 1

    close all
    startINDS = transpose(round(linspace(1,height(pointTable),height(dlc_lablab2useVIDE))));
    stopINDS = startINDS(2:end)-1;
    startINDS = startINDS(1:end-1);
    allXvaLS = table2array(pointTable(:,contains(pointTable.Properties.VariableNames,{'_x'})));
    allYvaLS = table2array(pointTable(:,contains(pointTable.Properties.VariableNames,{'_y'})));
    figure
    set(gcf,'Position',[178   311   382   991])
    tiledlayout(3,1,"TileSpacing","tight","Padding","tight")

    frameIIx2f = reshape(allXvaLS,numel(allXvaLS),1);
    frameIIy2f = reshape(allYvaLS,numel(allYvaLS),1);
    frameIIz = repmat(transpose(1:height(allYvaLS)),1,11);
    frameIIz2f = reshape(frameIIz,numel(frameIIz),1);

    rowIndices = repmat(1:size(plasCMP2u,1), height(allYvaLS), 1);
    rowIndices = rowIndices(:);
    plasCMP2u2u = plasCMP2u(rowIndices, :);

    nexttile
    ssALL1 = scatter3(frameIIz2f,frameIIx2f,frameIIy2f,60,plasCMP2u2u,'filled');
    ssALL1.MarkerFaceAlpha = useALPHA;
    view([41 34])
    xticklabels([])
    yticklabels([])
    zticklabels([])
    ylabel('X')
    xlabel('Frame number')
    zlabel('Y')
    box off
    axis tight


    nexttile
    ssALL2 = scatter3(frameIIz2f,frameIIx2f,frameIIy2f,60,plasCMP2u2u,'filled');
    ssALL2.MarkerFaceAlpha = useALPHA;
    view([9 10])
    xticklabels([])
    yticklabels([])
    zticklabels([])
    ylabel('X')
    xlabel('Frame number')
    zlabel('Y')
    box off
    % axis vis3d


    nexttile
    ssALL3 = scatter3(frameIIz2f,frameIIx2f,frameIIy2f,60,plasCMP2u2u,'filled');
    ssALL3.MarkerFaceAlpha = useALPHA;
    view([91 45])
    xticklabels([])
    yticklabels([])
    zticklabels([])
    ylabel('X')
    xlabel('Frame number')
    zlabel('Y')
    box off
    % axis vis3d





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



function [xDims , yDims] = getDataDims2(xVals , yVals)


xDims.min = min(xVals,[],'all');
xDims.max = max(xVals,[],'all');

yDims.min = min(yVals,[],'all');
yDims.max = max(yVals,[],'all');


end