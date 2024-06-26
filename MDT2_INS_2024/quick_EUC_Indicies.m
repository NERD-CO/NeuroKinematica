

csv_vidLoc = 'E:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data\MDT8\RightBody';
save_matLoc = 'E:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data\MDT8\RightBody';
GROUP_Video_CSV = 't2_20230816_idea08_session003_rightCam-0000DLC_resnet50_INS_2024_MPR8_LSTNMar20shuffle1_100000.csv';
dlcIO_processCSV('dirLOC',1,'userLOC',csv_vidLoc,'saveLOC',save_matLoc,...
    'selCSV',GROUP_Video_CSV,'USERid','JAT','hemiS','L','bodyS','R')


%% mdt7 baseline
matFileN = 'dlcDAT_20230908_session003_leftCam_LSTN_RBODY.mat';
moveCSV = '20230908_idea08_session003_leftCam-0000DLC_resnet50_MDT2_MPR9Sep24shuffle1_100000_labeled_MoveIndex.csv';
%

load(matFileN , 'outDATA')
tFname = fieldnames(outDATA.labelTab);
tmpTab = outDATA.labelTab.(tFname{1});

colNames = tmpTab.Properties.VariableNames; %
colNames2 = cellfun(@(x) split(x,'_'), colNames,...
    'UniformOutput',false);
colNames3 = unique(cellfun(@(x) x{1}, colNames2,...
    'UniformOutput',false));
colNames4 = colNames3(~matches(colNames3,'frames'));

% Initialize 'euclidall' to store Euclidean distances between successive points
euclidall = zeros(height(tmpTab)-1,length(colNames4));

% Iterate over each label and compute Euclidean distance for each frame
for label_i = 1:length(colNames4)

    tmpLabel_x = [colNames4{label_i} , '_x'];
    tmpLabel_y = [colNames4{label_i} , '_y'];

    tmpXdata = tmpTab.(tmpLabel_x);
    tmpYdata = tmpTab.(tmpLabel_y);

    labelData = [tmpXdata , tmpYdata];

    for frame_i = 1:height(labelData)
        if frame_i ~= height(labelData)
            point1 = labelData(frame_i,:);
            point2 = labelData(frame_i + 1,:);
            euclidall(frame_i , label_i) = pdist2(point1 , point2);
        end
    end
end

euclidalltmp = euclidall;

mall = mean(mean(euclidalltmp));
sall = std(std(euclidalltmp));
artThresh = mall + (sall*2);

euclidalltmp2 = euclidalltmp;
euclidalltmp2(euclidalltmp > artThresh) = 0;

rangeAll = zeros(1,width(euclidalltmp2));
for rrii = 1:width(euclidalltmp2)
    rangeAll(rrii) = range(euclidalltmp2(:,rrii));
end

% Create y data matrix
yPLOTDATA = euclidalltmp2;
pCURRENT = 0;
yLABELSi = zeros(width(euclidalltmp2),1);
yLABELSn = cell(width(euclidalltmp2),1);
for piip = 1:width(euclidalltmp2)

    if piip ~= 1
        pBump = rangeAll(piip-1);
        pCURRENT = pCURRENT + pBump;
    end

    tmpDATA = euclidalltmp2(:,piip) + pCURRENT;
    yPLOTDATA(:,piip) = tmpDATA;

    yLABELSi(piip) = pCURRENT;
    yLABELSn{piip} = ['Marker ', num2str(piip)];
end

MarkerPlotData = yPLOTDATA;
MarkerYlabels = yLABELSn;
MarkerYindicies = yLABELSi;

plot(MarkerPlotData,'Color',[0 0 0 0.5])

xlim([1 height(MarkerPlotData)])

% fix y lim
ylim([0 max(MarkerPlotData,[],'all')+5])

% fix labels
yticks(MarkerYindicies)
yticklabels(MarkerYlabels)

% load CSV

moveTAB = readtable(moveCSV);
emptyROWS = moveTAB.MoveN == 0;
moveTAB2 = moveTAB(~emptyROWS,:);
% Loop through everything but LFP TABLET and REST
moveTypeR = moveTAB2(~matches(moveTAB2.MoveType,{'LFP TABLET','REST'}),:);


colorS = rand(15,3);

hold on
for mii = 1:height(moveTypeR)

    xline(moveTypeR.BeginF(mii),'-','Color',colorS(mii,:),'LineWidth',3)
    xline(moveTypeR.EndF(mii),'-','Color',colorS(mii,:),'LineWidth',3)

end
