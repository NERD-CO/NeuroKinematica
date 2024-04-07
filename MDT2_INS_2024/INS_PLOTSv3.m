% generateINS_DatasetS_v2('MDT9' , 'right')
% clear
% close all
cd('E:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data\finalResults')
% load("MDT7_left.mat")

%%
matdir1 = dir('*.mat');
matdir2 = {matdir1.name};

meanAallbeta = zeros(4,52);
meanBallbeta = zeros(4,52);
meanAallkin = zeros(4,8);
meanBallkin = zeros(4,8);
for midi = 1:length(matdir2)

    tmpF = matdir2{midi};
    load(tmpF,'outDATAFin');
     

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

    meanAallkin(midi,:) = mean(allmoveAn);
    meanBallkin(midi,:) =  mean(allmoveBn);



    betaNormalized1 = [allBetaA ; allBetaB];

    betaNormalized2 = (betaNormalized1 - mean(betaNormalized1,'all'))./std(betaNormalized1,[],'all');

    allBetaAn = betaNormalized2(1:height(allBetaA),:);
    allBetaBn = betaNormalized2(height(allBetaA)+1:height(betaNormalized2),:);

    meanAallbeta(midi,:) = mean(allBetaAn);
    meanBallbeta(midi,:) = mean(allBetaBn);

    xDATAbet =  outDATAFin.lfpKINdata.GroupA.betapsdall2{biA}(:,2);

end

%%

close all
%%%%% PLOT 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tiledlayout(1,2,"TileSpacing","tight")

set(gcf,'Position',[161 342 934 724])

nexttile
plot(xDATAbet,meanAallbeta,'Color',[1 0 0 0.4],'LineWidth',1);
hold on
plot(xDATAbet,mean(meanAallbeta),'r','LineWidth',3)
plot(xDATAbet,meanBallbeta,'Color',[0 0 1 0.4],'LineWidth',1);
plot(xDATAbet,mean(meanBallbeta),'b','LineWidth',3)

legend({'','','',outDATAFin.conditionID.GroupA ,...
            outDATAFin.conditionID.GroupB});

groupAmean = mean(meanAallbeta);
groupBmean = mean(meanBallbeta);
xlabel('Frequency Hz')
ylabel('Normalized power')
title('LFP')
[~,beta_pval,~] = kstest2(groupAmean,groupBmean);

xlim([0.5 50])
axis square

text(25,1.5,['p = ',num2str(round(beta_pval,3))],'FontSize',14)
box off

%%%%% PLOT 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile

plot(xDATAmov,meanAallkin,'Color',[1 0 0 0.4],'LineWidth',1);
hold on
plot(xDATAmov,mean(meanAallkin),'r','LineWidth',3)
plot(xDATAmov,meanBallkin,'Color',[0 0 1 0.4],'LineWidth',1);
plot(xDATAmov,mean(meanBallkin),'b','LineWidth',3)

% legend({outDATAFin.conditionID.GroupA ,...
%             outDATAFin.conditionID.GroupB});

groupAmean = mean(meanAallkin);
groupBmean = mean(meanBallkin);

[~,move_pval,~] = kstest2(groupAmean,groupBmean);
xlabel('Frequency Hz')
ylabel('Normalized power')
xlim([2 9])
title('Kinematic')
axis square
box off

text(5,-1.5,['p = ',num2str(round(move_pval,3))],'FontSize',14)