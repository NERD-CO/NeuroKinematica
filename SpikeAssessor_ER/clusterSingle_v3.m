function [ nextInds , exceedflag , updatedID] = clusterSingle_v3( waveForms, clustInd, clustID, Fs, newID)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

tempWaves1 = transpose(waveForms);

numCclusts = numel(clustID);

memIND = clustInd;

for ci = 1:numCclusts
    
    clu2use = clustID(ci);
    
    tempWaves2 = tempWaves1(clustInd == clu2use,:);
    
    shMemIND = clustInd == clu2use;
    
    [ xFeats ] = createSpkFeatures2( transpose(tempWaves2),  Fs , 0 );
    
    [~,b] = pca(xFeats);

    [idxF,~] = kmeans(b(:,1:2),2,'Distance','cityblock',...
        'Replicates',6);
    
    [nextInds , exceedflag , updatedID] = getNextInd(clu2use , memIND , shMemIND, idxF, newID);
    
    memIND = nextInds;
    
end



end

function [nextInds , exceedflag ,updatedID] = getNextInd(orgClust , memIND , locOrg, locNew , newID)

% Original Cluster Index
% clustInd

% Location of original cluster
% locOrg
% Location of new clusters
% locNew
% IDs for new clusters = original ID + min of available from 1:6

locOrgRep(locOrg,:) = locNew;

allValues = unique(memIND);

if numel(allValues) == 6
    
    nextInds = memIND;
    exceedflag = 1;
else
    
    newIDS = [orgClust , newID];
    newID = newID + 1;
    
    newLocInds = unique(locNew);
    
    for ni = 1:numel(newLocInds)
        
        curLocs = locOrgRep == newLocInds(ni);
        
        memIND(curLocs) = newIDS(ni);
        
    end
    
    nextInds = memIND;
    exceedflag = 0;
    
end

updatedID = newID;

end














