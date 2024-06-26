% generateINS_DatasetS_v2('MDT9' , 'right')
% clear
% close all

getPCid = getenv('COMPUTERNAME');
switch getPCid
    case 'DESKTOP-FAGRV5G'
        cd('E:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data\finalResults')
    otherwise
        cd('D:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data\finalResults')
end
% load("MDT7_left.mat")

%%
matdir1 = dir('*.mat');
matdir2 = {matdir1.name};

%%
stackAallbeta = nan(100,52);
stackBallbeta = nan(100,52);
stackAallkin = nan(100,8);
stackBallkin = nan(100,8);
meanAallbeta = zeros(4,52);
meanBallbeta = zeros(4,52);
meanAallkin = zeros(4,8);
meanBallkin = zeros(4,8);
stackAallKINids = cell(100,1);
stackBallKINids = cell(100,1);
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

    curNANbA = find(isnan(stackAallbeta(:,1)),1,'first');
    curNANbB = find(isnan(stackBallbeta(:,1)),1,'first');
    stackAallbeta(curNANbA:curNANbA+height(allBetaAn)-1,:) = allBetaAn;
    stackBallbeta(curNANbB:curNANbB+height(allBetaBn)-1,:) = allBetaBn;

    meanAallbeta(midi,:) = mean(allBetaAn);
    meanBallbeta(midi,:) = mean(allBetaBn);

    curNANkA = find(isnan(stackAallkin(:,1)),1,'first');
    curNANkB = find(isnan(stackBallkin(:,1)),1,'first');
    stackAallkin(curNANkA:curNANkA+height(allmoveAn)-1,:) = allmoveAn;
    stackBallkin(curNANkB:curNANkB+height(allmoveBn)-1,:) = allmoveBn;

    curMOVEidsA = outDATAFin.lfpKINdata.GroupA.moveID;
    curMOVEidsB = outDATAFin.lfpKINdata.GroupB.moveID;

    emCELLSA = find(cellfun(@(x) isempty(x), stackAallKINids,'UniformOutput',true),...
        1,'first');
    stackAallKINids(emCELLSA:emCELLSA+length(curMOVEidsA)-1,1) = curMOVEidsA;

    emCELLSB = find(cellfun(@(x) isempty(x), stackBallKINids,'UniformOutput',true),...
        1,'first');
    stackBallKINids(emCELLSB:emCELLSB+length(curMOVEidsB)-1,1) = curMOVEidsB;

    xDATAbet =  outDATAFin.lfpKINdata.GroupA.betapsdall2{biA}(:,2);

end

%% clean up
stackAallbeta = stackAallbeta(~isnan(stackAallbeta(:,1)),:);
stackBallbeta = stackBallbeta(~isnan(stackBallbeta(:,1)),:);
stackAallkin = stackAallkin(~isnan(stackAallkin(:,1)),:);
stackBallkin = stackBallkin(~isnan(stackBallkin(:,1)),:);
stackAallKINids = stackAallKINids(cellfun(@(x) ~isempty(x), stackAallKINids,...
    'UniformOutput',true));
stackBallKINids = stackBallKINids(cellfun(@(x) ~isempty(x), stackBallKINids,...
    'UniformOutput',true));

%% Move into Move buckets

% Complete 5/5/2024
GroupAKINids = unique(stackAallKINids);
GroupBKINids = unique(stackBallKINids);

logINDcellA = cell(1,3);
for uniI = 1:3
    logINDcellA{uniI} = matches(stackAallKINids,GroupAKINids{uniI});
end
groupAidTab = cell2table(logINDcellA,'VariableNames',GroupAKINids);

logINDcellB = cell(1,3);
for uniI = 1:3
    logINDcellB{uniI} = matches(stackBallKINids,GroupBKINids{uniI});
end
groupBidTab = cell2table(logINDcellB,'VariableNames',GroupBKINids);

%% NEW SUBmove PLOT - 5/6/2024

% TO DO - fix LFP y-axis
% To Do - fix KIN y-axis
% to do - add subtitle

MBSColor = [219, 175, 33]/255;
TMRColor = [17, 17, 15]/255;

close all
tiledlayout(2,3,"TileSpacing","tight")

set(gcf,'Position',[161 342 934 724])

nexttile % LFP 1 - finger tap
fingerTapAbeta = stackAallbeta(groupAidTab.("FINGER TAP"){1},:);
fingerTapBbeta = stackBallbeta(groupBidTab.("FINGER TAP"){1},:);
plot(xDATAbet,fingerTapAbeta,'Color',TMRColor,'LineWidth',1);
hold on
plot(xDATAbet,mean(fingerTapAbeta),'Color',TMRColor,'LineWidth',3)
plot(xDATAbet,fingerTapBbeta,'Color',MBSColor,'LineWidth',1);
plot(xDATAbet,mean(fingerTapBbeta),'Color',MBSColor,'LineWidth',3)
ax = gca;
title('LFP: Finger tap');
ax.TitleHorizontalAlignment = 'left'; 
ylim([-1 5])

legend({'','','','',outDATAFin.conditionID.GroupA ,...
    '','','',outDATAFin.conditionID.GroupB});

fingertAmean = mean(fingerTapAbeta);
fingertBmean = mean(fingerTapBbeta);
xlabel('Frequency Hz')
ylabel('Z-scored power')
[~,beta1_pval,~] = kstest2(fingertAmean,fingertBmean);

xlim([0.5 50])
xticks([1 13 30 50])
axis square
text(25,1.5,['p = ',num2str(round(beta1_pval,3))],'FontSize',14)
box off

nexttile % LFP 2 - hand oc
handOCAbeta = stackAallbeta(groupAidTab.("HAND OC"){1},:);
handOCBbeta = stackBallbeta(groupBidTab.("HAND OC"){1},:);
plot(xDATAbet,handOCAbeta,'Color',TMRColor,'LineWidth',1);
hold on
plot(xDATAbet,mean(handOCAbeta),'Color',TMRColor,'LineWidth',3)
plot(xDATAbet,handOCBbeta,'Color',MBSColor,'LineWidth',1);
plot(xDATAbet,mean(handOCBbeta),'Color',MBSColor,'LineWidth',3)
ax = gca;
title('LFP: Hand Open/Close');
ax.TitleHorizontalAlignment = 'left'; 
ylim([-1 5])

handOCAmean = mean(handOCAbeta);
handOCBmean = mean(handOCBbeta);
[~,beta2_pval,~] = kstest2(handOCAmean,handOCBmean);
xlim([0.5 50])
xticks([1 13 30 50])
axis square
text(25,1.5,['p = ',num2str(round(beta2_pval,3))],'FontSize',14)
box off

nexttile % LFP 3 - hand ps
handPSAbeta = stackAallbeta(groupAidTab.("HAND PS"){1},:);
handPSBbeta = stackBallbeta(groupBidTab.("HAND PS"){1},:);
plot(xDATAbet,handPSAbeta,'Color',TMRColor,'LineWidth',1);
hold on
plot(xDATAbet,mean(handPSAbeta),'Color',TMRColor,'LineWidth',3)
plot(xDATAbet,handPSBbeta,'Color',MBSColor,'LineWidth',1);
plot(xDATAbet,mean(handPSBbeta),'Color',MBSColor,'LineWidth',3)
ax = gca;
title('LFP: Hand Pronation/Supination');
ax.TitleHorizontalAlignment = 'left'; 
ylim([-1 5])

handPSAmean = mean(handPSAbeta);
handPSBmean = mean(handPSBbeta);
[~,beta3_pval,~] = kstest2(handPSAmean,handPSBmean);
xlim([0.5 50])
xticks([1 13 30 50])
axis square
text(25,1.5,['p = ',num2str(round(beta2_pval,3))],'FontSize',14)
box off


moveMENTs = {"FINGER TAP","HAND OC","HAND PS"};

for kii = 1:3

    nexttile % KIN 1 - finger tap
    moveAKin = stackAallkin(groupAidTab.(moveMENTs{kii}){1},:);
    moveBKin = stackBallkin(groupBidTab.(moveMENTs{kii}){1},:);
    plot(xDATAmov,moveAKin,'Color',TMRColor,'LineWidth',1);
    hold on
    plot(xDATAmov,mean(moveAKin),'Color',TMRColor,'LineWidth',3)
    plot(xDATAmov,moveBKin,'Color',MBSColor,'LineWidth',1);
    plot(xDATAmov,mean(moveBKin),'Color',MBSColor,'LineWidth',3)
    ax = gca;
    switch kii
        case 1
            title('Kinematic: Finger tap');
        case 2
            title('Kinematic: Hand Open/Close');
        case 3
            title('Kinematic: Hand Pronation/Supination');
    end
    ax.TitleHorizontalAlignment = 'left';
    % ylim([-1 5])

    if kii == 1

        legend({'','','','',outDATAFin.conditionID.GroupA ,...
            '','','',outDATAFin.conditionID.GroupB});

        xlabel('Frequency Hz')
        ylabel('Z-scored power')

    end

    moveAmean = mean(moveAKin);
    moveBmean = mean(moveBKin);

    [~,kin_pval,~] = kstest2(moveAmean,moveBmean);

    xlim([2 9])
    % xticks([1 13 30 50])
    axis square
    text(5,-1.5,['p = ',num2str(round(kin_pval,3))],'FontSize',14)
    box off

end


%%

MBSColor = [219, 175, 33]/255;
TMRColor = [17, 17, 15]/255;

close all
%%%%% PLOT 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tiledlayout(1,2,"TileSpacing","tight")

set(gcf,'Position',[161 342 934 724])

nexttile
plot(xDATAbet,meanAallbeta,'Color',TMRColor,'LineWidth',1);
hold on
plot(xDATAbet,mean(meanAallbeta),'Color',TMRColor,'LineWidth',3)
plot(xDATAbet,meanBallbeta,'Color',MBSColor,'LineWidth',1);
plot(xDATAbet,mean(meanBallbeta),'Color',MBSColor,'LineWidth',3)

legend({'','','','',outDATAFin.conditionID.GroupA ,...
    '','','','',outDATAFin.conditionID.GroupB});

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

plot(xDATAmov,meanAallkin,'Color',TMRColor,'LineWidth',1);
hold on
plot(xDATAmov,mean(meanAallkin),'Color',TMRColor,'LineWidth',3)
plot(xDATAmov,meanBallkin,'Color',MBSColor,'LineWidth',1);
plot(xDATAmov,mean(meanBallkin),'Color',MBSColor,'LineWidth',3)

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