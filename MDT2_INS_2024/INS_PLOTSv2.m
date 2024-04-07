generateINS_DatasetS_v2('MDT9' , 'right')
clear
close all
% cd('E:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data\finalResults')
% load("MDT7_left.mat")

%%


%%%%% PLOT 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1)
allBetaA = zeros(height(outDATAFin.lfpKINdata.GroupA.betapsdall2),...
    size(outDATAFin.lfpKINdata.GroupA.betapsdall2{1},1));

allBetaB = zeros(height(outDATAFin.lfpKINdata.GroupB.betapsdall2),...
    size(outDATAFin.lfpKINdata.GroupB.betapsdall2{1},1));

for biA = 1:height(outDATAFin.lfpKINdata.GroupA.betapsdall2)
    tbeA = outDATAFin.lfpKINdata.GroupA.betapsdall2{biA}(:,1);
    allBetaA(biA,:) = tbeA;
end

for biB = 1:height(outDATAFin.lfpKINdata.GroupB.betapsdall2)
    tbeB = outDATAFin.lfpKINdata.GroupA.betapsdall2{biB}(:,1);
    allBetaB(biB,:) = tbeB;
end

betaNormalized1 = [allBetaA ; allBetaB];

betaNormalized2 = (betaNormalized1 - mean(betaNormalized1,'all'))./std(betaNormalized1,[],'all');

allBetaAn = betaNormalized2(1:height(allBetaA),:);
allBetaBn = betaNormalized2(height(allBetaA)+1:height(betaNormalized2),:);

xDATAbet =  outDATAFin.lfpKINdata.GroupA.betapsdall2{biA}(:,2);

plot(xDATAbet,transpose(allBetaAn),'r')
hold on
plot(xDATAbet,mean(allBetaAn),'r','LineWidth',3);
plot(xDATAbet,transpose(allBetaBn),'b')
plot(xDATAbet,mean(allBetaBn),'b','LineWidth',3);

legend({'','','','',outDATAFin.conditionID.GroupA ,...
             '','','',outDATAFin.conditionID.GroupB})

xlim([0.5 50])

%%%%% PLOT 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2)
allMoveA = zeros(height(outDATAFin.lfpKINdata.GroupA.movepsdall2),...
    size(outDATAFin.lfpKINdata.GroupA.movepsdall2{1},1));

allMoveB = zeros(height(outDATAFin.lfpKINdata.GroupB.movepsdall2),...
    size(outDATAFin.lfpKINdata.GroupB.movepsdall2{1},1));

for miA = 1:height(outDATAFin.lfpKINdata.GroupA.movepsdall2)
    tmvA = outDATAFin.lfpKINdata.GroupA.movepsdall2{miA}(:,1);
    allMoveA(miA,:) = tmvA;
end

for miB = 1:height(outDATAFin.lfpKINdata.GroupB.movepsdall2)
    tmvB = outDATAFin.lfpKINdata.GroupB.movepsdall2{miB}(:,1);
    allMoveB(miB,:) = tmvB;
end

moveNormalized1 = [allMoveA ; allMoveB];

moveNormalized2 = (moveNormalized1 - mean(moveNormalized1,'all'))./std(moveNormalized1,[],'all');

allmoveAn = moveNormalized2(1:height(allMoveA),:);
allmoveBn = moveNormalized2(height(allMoveA)+1:height(moveNormalized2),:);

xDATAmov =  outDATAFin.lfpKINdata.GroupA.movepsdall2{miA}(:,2);

plot(xDATAmov,transpose(allmoveAn),'r')
hold on
plot(xDATAmov,mean(allmoveAn),'r','LineWidth',3);
plot(xDATAmov,transpose(allmoveBn),'b')
plot(xDATAmov,mean(allmoveBn),'b','LineWidth',3);

legend({'','','','',outDATAFin.conditionID.GroupA ,...
             '','','',outDATAFin.conditionID.GroupB})