fs = 250;

[pxxMM,fMM] = pwelch(pointX2sm,250,125,250,250);
p2bmm = pow2db(pxxMM);
p2bMMs = smoothdata(p2bmm,'gaussian',2);
fMMi = fMM > 1 & fMM < 10;
fMMf = fMM(fMMi);
pxxMf = p2bMMs(fMMi);

% figure;
% plot(fMMf,pxxMf)




% plot(f,10*log10(pxx))
%
% xlabel('Frequency (Hz)')
% ylabel('PSD (dB/Hz)')
%
% xlim([1 10])


%%

[Pupvxx,upVfxx] = pwelch(transpose(movLFP), hanning(250), 125, 256, 250, 'onesided');
uVp_t = sqrt(Pupvxx).*rms(hanning(250)).*sqrt(2).*2.*250/256;
% PxxP = 10*log10(Pxx);
% uPv_A = uVp_t;

uPv_As = smoothdata(uVp_t,'gaussian',5);
% pxxps = smoothdata(PxxP,'gaussian',5);

upv50 = upVfxx < 50;
upVfxxU = upVfxx(upv50);
uVp_tU = uPv_As(upv50);


figure;
plot(upVfxxU,uVp_tU)
% xlim([0 50])

% figure;
% plot(pxxps)
% xlim([0 50])
hold on
plot(uPvRESTfreq , uPvRESTpow)
legend('epoch','rest')


%%
generateINS_DatasetS_v2('MDT7' , 'left')
% PCA plot of movement kin
% REMOVE LIKELIHOOD and PUSH THROUGH
colNames = dlc_lab2use2int.Properties.VariableNames; %
colNames2 = cellfun(@(x) split(x,'_'), colNames,...
    'UniformOutput',false);
colNames3 = unique(cellfun(@(x) x{1}, colNames2,...
    'UniformOutput',false));
colNames4 = colNames3(~matches(colNames3,'frames'));

optItemInd = colNames4{tmpMOve.OptItem};
tmpMOve.MoveType{1}
xNAMEot = [optItemInd,'_x'];
yNAMEot = [optItemInd,'_y'];
startINDi = startINDi + 50;
endINDi = endINDi - 50;
trimFRAMES = dlc_lab2use2int(startINDi:endINDi,:);

% REMOVE likelihood
colKEEP = contains(trimFRAMES.Properties.VariableNames,{'_x','_y'});
trimFRAMES2 = trimFRAMES(:,colKEEP);

[a,b,c,d] = pca(table2array(trimFRAMES2));

top3 = b(:,1:3);

figure;
plot3(top3(:,1),top3(:,2),top3(:,3))
figure;
plot(top3(:,1))







figure('NumberTitle', 'off', 'Name', 'test')


% remove none fingers;

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
for fi = 1:length(uniNames)

    tUiName = uniNames{fi};
    uiNindex = ismember(cVnames,tUiName);
    uiTab = pointTable(:,uiNindex);

    % Get X and Y
    tPnames0 = cellfun(@(x) strsplit(x,'_'),...
        uiTab.Properties.VariableNames, 'UniformOutput', false);
    tPnames = cellfun(@(x) x{2}, tPnames0, 'UniformOutput', false);
    xyIND = contains(tPnames,{'x','y'});
    xyiTab = table2array(uiTab(:,xyIND));

    allFrames = xyiTab;
    allColors = repmat(plasCMP2u(fi,:),height(uiTab) ,1);

    % scatter(allFrames(:,1),allFrames(:,2),30,allColors,'filled')
    scatter3(transpose(1:length(allFrames)),allFrames(:,1),allFrames(:,2),30,allColors,'filled')
    hold on

end
% set view
view([41 34])

% plot the thing





%%

% plot the video
% t1_20230816_idea08_session005_rightCam-0000DLC_resnet50_INS_2024_MPR7_LSTNMar20shuffle1_100000_labeled.mp4
% dlcLab_vidLoc = vidLoc;
% cd(dlcLab_vidLoc)

dlc_lab_vidObj = VideoReader('t1_20230816_idea08_session005_rightCam-0000DLC_resnet50_INS_2024_MPR7_LSTNMar20shuffle1_100000_labeled.mp4');
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

numFRAMES = floor(1695/433);
startINDS = transpose(round(linspace(1,1695,433)));
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


for fi = 1:length(uniNames)

    tUiName = uniNames{fi};
    uiNindex = ismember(cVnames,tUiName);
    uiTab = pointTable(:,uiNindex);

    % Get X and Y
    tPnames0 = cellfun(@(x) strsplit(x,'_'),...
        uiTab.Properties.VariableNames, 'UniformOutput', false);
    tPnames = cellfun(@(x) x{2}, tPnames0, 'UniformOutput', false);
    xyIND = contains(tPnames,{'x','y'});
    xyiTab = table2array(uiTab(:,xyIND));

    allFrames = xyiTab;
    allColors = repmat(plasCMP2u(fi,:),height(uiTab) ,1);

    % scatter(allFrames(:,1),allFrames(:,2),30,allColors,'filled')
    scatter3(transpose(1:length(allFrames)),allFrames(:,1),allFrames(:,2),30,allColors,'filled')
    hold on

end
% set view
view([41 34])



%%
tailLen = 10;
figure;
[xDims , yDims] = getDataDims(pointTable);

frameTa = nan(round(height(pointTable)*length(uniNames)),2);
colorsTa = nan(round(height(pointTable)*length(uniNames)),3);

start = 1;
stop = length(uniNames);

for fi = 1:height(pointTable)

    % Extract each frame
    tFrame = pointTable(fi,:);

    tPnames0 = cellfun(@(x) strsplit(x,'_'),...
        tFrame.Properties.VariableNames, 'UniformOutput', false);
    tPnames = cellfun(@(x) x{2}, tPnames0, 'UniformOutput', false);
    xyIND = contains(tPnames,{'x','y'});

    xyiTab = table2array(tFrame(:,xyIND));

    frameTx = xyiTab(1:2:length(xyiTab));
    frameTy = xyiTab(2:2:length(xyiTab));

    % Add to ALL frames

    frameTa(start:stop,:) = transpose([frameTx ; frameTy]);
    colorsTa(start:stop,:) = plasCMP2u;

    start = stop + 1;
    stop = stop + length(uniNames);

    % Remove NANs
    frameTaP = frameTa(~isnan(frameTa(:,1)),:);
    colorsTaP = colorsTa(~isnan(colorsTa(:,1)),:);

    if fi == 1

        %         scT = scatter(frameTaP(:,1),frameTaP(:,2),30,colorsTaP,'filled');
        scatter(frameTaP(:,1),frameTaP(:,2),30,colorsTaP,'filled');
        xlim([xDims.min xDims.max]);
        ylim([yDims.min yDims.max]);

        pause(0.02)


    else

        maxLENG = round(tailLen*length(uniNames)) + 1; % This is where tail length is determined
        curLEN = size(frameTaP,1);

        if size(frameTaP,1) > maxLENG

            newStart = size(frameTaP,1) - round(tailLen*length(uniNames));
            frameTaP = frameTaP(newStart:curLEN,:);
            colorsTaP = colorsTaP(newStart:curLEN,:);

        end

        %         scT.XData = frameTaP(:,1);
        %         scT.YData = frameTaP(:,2);
        %         scT.CData = colorsTaP;

        scatter(frameTaP(:,1),frameTaP(:,2),30,colorsTaP,'filled');
        xlim([xDims.min xDims.max]);
        ylim([yDims.min yDims.max]);

        pause(0.01)

    end % End of number 1 check
end % End of for loop

%%
% Get X and Y
figure;
tPnames0 = cellfun(@(x) strsplit(x,'_'),...
    pointTable.Properties.VariableNames, 'UniformOutput', false);
tPnames = cellfun(@(x) x{2}, tPnames0, 'UniformOutput', false);
xIND = contains(tPnames,'x');
yIND = contains(tPnames,'y');

xiTab = table2array(pointTable(:,xIND));
yiTab = table2array(pointTable(:,yIND));

xyAiTab = [reshape(xiTab,numel(xiTab),1) , reshape(yiTab,numel(yiTab),1)];

% START OF ANIMATION

[xDims , yDims] = getDataDims(pointTable);

frameTa = nan(round(height(pointTable)*length(uniNames)),2);
colorsTa = nan(round(height(pointTable)*length(uniNames)),3);

start = 1;
stop = length(uniNames);

for fi = 1:height(pointTable)

    % Extract each frame
    tFrame = pointTable(fi,:);

    tPnames0 = cellfun(@(x) strsplit(x,'_'),...
        tFrame.Properties.VariableNames, 'UniformOutput', false);
    tPnames = cellfun(@(x) x{2}, tPnames0, 'UniformOutput', false);
    xyIND = contains(tPnames,{'x','y'});

    xyiTab = table2array(tFrame(:,xyIND));

    frameTx = xyiTab(1:2:length(xyiTab));
    frameTy = xyiTab(2:2:length(xyiTab));

    % Add to ALL frames

    frameTa(start:stop,:) = transpose([frameTx ; frameTy]);
    colorsTa(start:stop,:) = plasCMP2u;

    start = stop + 1;
    stop = stop + length(uniNames);

    % Remove NANs
    frameTaP = frameTa(~isnan(frameTa(:,1)),:);
    colorsTaP = colorsTa(~isnan(colorsTa(:,1)),:);

    if fi == 1

        grSc = scatter(xyAiTab(:,1),xyAiTab(:,2),25,[0.5 0.5 0.5], 'filled');
        grSc.MarkerFaceAlpha = 0.5;
        hold on
        scT = scatter(frameTaP(:,1),frameTaP(:,2),30,colorsTaP,'filled');

        xlim([xDims.min xDims.max]);
        ylim([yDims.min yDims.max]);

        pause(0.02)

    else

        scT.XData = frameTaP(:,1);
        scT.YData = frameTaP(:,2);
        scT.CData = colorsTaP;

        xlim([xDims.min xDims.max]);
        ylim([yDims.min yDims.max]);

        pause(0.01)

    end % End of number 1 check
end % End of for loop

%%
figure;
[xDims , yDims] = getDataDims(pointTable);

frameTa = nan(round(height(pointTable)*length(uniNames)),2);
colorsTa = nan(round(height(pointTable)*length(uniNames)),3);

start = 1;
stop = length(uniNames);

for fi = 1:height(pointTable)

    % Extract each frame
    tFrame = pointTable(fi,:);

    tPnames0 = cellfun(@(x) strsplit(x,'_'),...
        tFrame.Properties.VariableNames, 'UniformOutput', false);
    tPnames = cellfun(@(x) x{2}, tPnames0, 'UniformOutput', false);
    xyIND = contains(tPnames,{'x','y'});

    xyiTab = table2array(tFrame(:,xyIND));

    frameTx = xyiTab(1:2:length(xyiTab));
    frameTy = xyiTab(2:2:length(xyiTab));

    % Add to ALL frames

    frameTa(start:stop,:) = transpose([frameTx ; frameTy]);
    colorsTa(start:stop,:) = plasCMP2u;

    start = stop + 1;
    stop = stop + length(uniNames);

    % Remove NANs
    frameTaP = frameTa(~isnan(frameTa(:,1)),:);
    colorsTaP = colorsTa(~isnan(colorsTa(:,1)),:);

    if fi == 1

        scT = scatter(frameTaP(:,1),frameTaP(:,2),30,colorsTaP,'filled');

        xlim([xDims.min xDims.max]);
        ylim([yDims.min yDims.max]);

        pause(0.02)

    else

        scT.XData = frameTaP(:,1);
        scT.YData = frameTaP(:,2);
        scT.CData = colorsTaP;

        xlim([xDims.min xDims.max]);
        ylim([yDims.min yDims.max]);

        pause(0.01)

    end % End of number 1 check
end % End of for loop


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