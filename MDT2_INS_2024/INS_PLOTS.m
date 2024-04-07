generateINS_DatasetS_v2('MDT9' , 'right')
load("MDT9_right.mat")

close all

groupAkin = outDATAFin.lfpKINdata.GroupA.groupKinAmp;
groupBkin = outDATAFin.lfpKINdata.GroupB.groupKinAmp;

groupAlfp = outDATAFin.lfpKINdata.GroupA.groupLFPAmp;
groupBlfp = outDATAFin.lfpKINdata.GroupB.groupLFPAmp;

tiledlayout(1,2)

nexttile
s1 = swarmchart(ones(length(groupAkin),1),groupAkin,'red');
s1.XJitter = "density";
s1.XJitterWidth = 0.2;
hold on
s2 = swarmchart(ones(length(groupBkin),1)+1,groupBkin,'blue');
s2.XJitter = "density";
s2.XJitterWidth = 0.2;

xticks(1:2)
xlim([0.5 2.5])
xticklabels({outDATAFin.conditionID.GroupA ,...
             outDATAFin.conditionID.GroupB})

ylabel('mean euclidean distance')

[pval1,~,~] = ranksum(groupAkin,groupBkin);

text(1.5 , round(mean([groupAkin ; groupBkin])),...
    ['p = ',num2str(round(pval1,3))])

nexttile
s1 = swarmchart(ones(length(groupAlfp),1),groupAlfp,'red');
s1.XJitter = "density";
s1.XJitterWidth = 0.2;
hold on
s2 = swarmchart(ones(length(groupBlfp),1)+1,groupBlfp,'blue');
s2.XJitter = "density";
s2.XJitterWidth = 0.2;

xticks(1:2)
xlim([0.5 2.5])
xticklabels({outDATAFin.conditionID.GroupA ,...
             outDATAFin.conditionID.GroupB})

ylabel('mean beta power')

[pval2,b,c] = ranksum(groupAlfp,groupBlfp);

text(1.5 , round(mean([groupAlfp ; groupBlfp])),...
    ['p = ',num2str(round(pval2,3))])



%%
close all
subplot(1,2,1)
allBetaA = zeros(height(outDATAFin.lfpKINdata.GroupA.betapsdall),...
    size(outDATAFin.lfpKINdata.GroupA.betapsdall{1},1));

allBetaB = zeros(height(outDATAFin.lfpKINdata.GroupB.betapsdall),...
    size(outDATAFin.lfpKINdata.GroupB.betapsdall{1},1));

for biA = 1:height(outDATAFin.lfpKINdata.GroupA.betapsdall)
    tbeA = outDATAFin.lfpKINdata.GroupA.betapsdall{biA}(:,1);
    allBetaA(biA,:) = tbeA;
end

for biB = 1:height(outDATAFin.lfpKINdata.GroupB.betapsdall)
    tbeB = outDATAFin.lfpKINdata.GroupA.betapsdall{biB}(:,1);
    allBetaB(biB,:) = tbeB;
end

betaNormalized1 = [allBetaA ; allBetaB];

betaNormalized2 = (betaNormalized1 - mean(betaNormalized1,'all'))./std(betaNormalized1,[],'all');

allBetaAn = betaNormalized2(1:height(allBetaA),:);
allBetaBn = betaNormalized2(height(allBetaA)+1:height(betaNormalized2),:);

xDATAbet =  outDATAFin.lfpKINdata.GroupA.betapsdall{biA}(:,2);

plot(xDATAbet,transpose(allBetaAn),'r')
hold on
plot(xDATAbet,mean(allBetaAn),'r','LineWidth',3);
plot(xDATAbet,transpose(allBetaBn),'b')
plot(xDATAbet,mean(allBetaBn),'b','LineWidth',3);

legend({'','','','',outDATAFin.conditionID.GroupA ,...
             '','','',outDATAFin.conditionID.GroupB})

xlim([0.5 50])

subplot(1,2,2)
allMoveA = zeros(height(outDATAFin.lfpKINdata.GroupA.movepsdall),...
    size(outDATAFin.lfpKINdata.GroupA.movepsdall{1},1));

allMoveB = zeros(height(outDATAFin.lfpKINdata.GroupB.movepsdall),...
    size(outDATAFin.lfpKINdata.GroupB.movepsdall{1},1));

for miA = 1:height(outDATAFin.lfpKINdata.GroupA.movepsdall)
    tmvA = outDATAFin.lfpKINdata.GroupA.movepsdall{miA}(:,1);
    allMoveA(miA,:) = tmvA;
end

for miB = 1:height(outDATAFin.lfpKINdata.GroupB.movepsdall)
    tmvB = outDATAFin.lfpKINdata.GroupB.movepsdall{miB}(:,1);
    allMoveB(miB,:) = tmvB;
end

moveNormalized1 = [allMoveA ; allMoveB];

moveNormalized2 = (moveNormalized1 - mean(moveNormalized1,'all'))./std(moveNormalized1,[],'all');

allmoveAn = moveNormalized2(1:height(allMoveA),:);
allmoveBn = moveNormalized2(height(allMoveA)+1:height(moveNormalized2),:);

xDATAmov =  outDATAFin.lfpKINdata.GroupA.movepsdall{miA}(:,2);

plot(xDATAmov,transpose(allmoveAn),'r')
hold on
plot(xDATAmov,mean(allmoveAn),'r','LineWidth',3);
plot(xDATAmov,transpose(allmoveBn),'b')
plot(xDATAmov,mean(allmoveBn),'b','LineWidth',3);

legend({'','','','',outDATAFin.conditionID.GroupA ,...
             '','','',outDATAFin.conditionID.GroupB})