


dlc_lab2use2int



colNames = dlc_lab2use2int.Properties.VariableNames; %
colNames2 = cellfun(@(x) split(x,'_'), colNames,...
    'UniformOutput',false);
colNames3 = unique(cellfun(@(x) x{1}, colNames2,...
    'UniformOutput',false));
colNames4 = colNames3(~matches(colNames3,'frames'));

% Initialize 'euclidall' to store Euclidean distances between successive points
euclidall = zeros(height(dlc_lab2use2int)-1,length(colNames4));

% Iterate over each label and compute Euclidean distance for each frame
for label_i = 1:length(colNames4)

    tmpLabel_x = [colNames4{label_i} , '_x'];
    tmpLabel_y = [colNames4{label_i} , '_y'];

    tmpXdata = dlc_lab2use2int.(tmpLabel_x);
    tmpYdata = dlc_lab2use2int.(tmpLabel_y);

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

colorSSS = rand(height(moveIItable),3);

for mii = 1:height(moveIItable)

    tmpMOve = moveIItable(mii,:);
    startIND = tmpMOve.BeginF;
    endIND = tmpMOve.EndF;

    [~ , startINDi] = min(abs(startIND - trimFrame_int));
    [~ , endINDi] = min(abs(endIND - trimFrame_int));

    startINDi2 = round(startIND*4.73) - round(180*4.73);
    endINDi2 = round(endIND*4.73) - round(180*4.73);

    xline(startINDi,'-','Color',colorSSS(mii,:))
    xline(endINDi,'-','Color',colorSSS(mii,:))

    xline(startINDi2,'-','Color',colorSSS(mii,:) ,'LineWidth',3)
    xline(endINDi2,'-','Color',colorSSS(mii,:),'LineWidth',3)



end