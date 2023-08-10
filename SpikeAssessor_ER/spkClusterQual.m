function [ nidx , lmetrics , oCen] = spkClusterQual( featData, clustID, clustCens,  stage )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



if stage == 1
    
    [lmetrics] = initialQual(featData, clustID);
    
    nidx = nan;
    oCen = nan;
    
elseif stage == 2
    
    if size(featData,1) < 30
        lmetrics.IsDs = NaN;
        lmetrics.LRat = NaN;
        nidx = ones(size(featData,1),1);
        oCen = computeCentroids(nidx, featData);
        return
    else
        
        [fmetrics] = initialQual(featData, clustID);
        
        [lmetrics , nidx, oCen] = finalQual(featData, fmetrics, clustID, clustCens);
        
    end
    
end




end

%%%% INITAL Quality Parameters function
function [metrics] = initialQual(normalX, clusID)

nSpikes = size(normalX,1);

IsolDist = zeros(max(clusID),1); % > 20
Lratio = zeros(max(clusID),1); % < 0.4

for i = 1:max(clusID)
    
    ClusterSpikes = find(clusID == i);
    
    if length(ClusterSpikes) < size(normalX,2)
        Lratio(i) = 0.8;
        IsolDist(i) = 5;
        continue
    end
    
    nClusterSpikes = length(ClusterSpikes);
    NoiseSpikes = setdiff(1:nSpikes, ClusterSpikes);
    m = mahal(normalX, normalX(ClusterSpikes,:));
    
%     mCluster = m(ClusterSpikes);
    mNoise = m(NoiseSpikes);
    
    df = size(normalX,2);
    
    L = sum(1-chi2cdf(m(NoiseSpikes),df));
    Lratio(i) = L/nClusterSpikes;
    
    
%     InClu = ismember(1:nSpikes, ClusterSpikes);
    if (nClusterSpikes < nSpikes/2)
        [sorted, ~] = sort(mNoise);
        IsolDist(i) = sorted(nClusterSpikes);
    else
        IsolDist(i) = NaN;
    end
    
    metrics.IsDs = IsolDist;
    metrics.LRat = Lratio;
    
end

end


%%%% FINAL Quality Parameters function
function [ metrics , nclID , oCens ] = finalQual(normalX, inMetrics, idxF, inCens)

nSpikes = size(normalX,1);

IsolDist = inMetrics.IsDs;
Lratio = inMetrics.LRat;

clustREFINE = true;

while clustREFINE
    
    if all(IsolDist < 20 | Lratio > 0.4)
        
        [nclID, oCens] = kmeans(normalX,2,'Distance','cityblock',...
                              'Replicates',5);
        
        [metrics] = initialQual(normalX, nclID);
        
        disp('All clusters are bad')
        
        return
        
    
    elseif any(IsolDist < 20 | Lratio > 0.4)
        
        % Take bad cluster points and reassign two kept centroids
        % Step 1
        
        badClustInd = find(IsolDist < 20 | Lratio > 0.4);
        badClustNum = numel(badClustInd);
        
        goodCens = inCens(~(IsolDist < 20 | Lratio > 0.4),:);
        goodCenID = find(~(IsolDist < 20 | Lratio > 0.4));
        inCens = inCens(goodCenID,:);
        %
        badClustInfo = struct;
        for bci = 1:badClustNum
            
            bciId = badClustInd(bci);
            
            badCltPts = find(idxF == bciId);
            
            reaSSign = zeros(length(badCltPts),3);
            for reAs = 1:numel(badCltPts)
                
                tempPoint = normalX(badCltPts(reAs),:);
                reaSSign(reAs,1) = idxF(badCltPts(reAs));
                
                cenDists = zeros(length(goodCenID),1);
                for cens = 1:size(goodCens,1)
                    
                    compDist = [goodCens(cens,:); tempPoint];
                    
                    cenDists(cens,1) = pdist(compDist,'euclidean');
                    
                end
                
                % Find minDist - assign dist and new cluster
                [minDist,minDid] = min(cenDists);
                reaSSign(reAs,2) = minDist;
                reaSSign(reAs,3) = goodCenID(minDid);
                
            end
            
            badClustInfo.(strcat('bdC_', num2str(bciId))).bdClID = bciId;
            badClustInfo.(strcat('bdC_', num2str(bciId))).bdPts = badCltPts;
            badClustInfo.(strcat('bdC_', num2str(bciId))).reCls = reaSSign;
            
        end
        
        % Reassign cluster values
        for nci = 1:length(fieldnames(badClustInfo))
            
            bdID = badClustInd(nci);
            ptInd = badClustInfo.(strcat('bdC_', num2str(bdID))).bdPts;
            idxF(ptInd) = badClustInfo.(strcat('bdC_', num2str(bdID))).reCls(:,3);
            
        end
        
        % Recompute cluster metrics
        IsolDist = zeros(numel(goodCenID),1); % > 20
        Lratio = zeros(numel(goodCenID),1); % < 0.4
        
        for i = 1:numel(goodCenID)
            
            cID = goodCenID(i);
            
            ClusterSpikes = find(idxF == cID);
            nClusterSpikes = length(ClusterSpikes);
            NoiseSpikes = setdiff(1:nSpikes, ClusterSpikes);
            m = mahal(normalX, normalX(ClusterSpikes,:));
            
%             mCluster = m(ClusterSpikes);
            mNoise = m(NoiseSpikes);
            
            df = size(normalX,2);
            
            L = sum(1-chi2cdf(m(NoiseSpikes),df));
            Lratio(i) = L/nClusterSpikes;
            
            %InClu = ismember(1:nSpikes, ClusterSpikes);
            if (nClusterSpikes < nSpikes/2)
                [sorted, ~] = sort(mNoise);
                IsolDist(i) = sorted(nClusterSpikes);
            else
                IsolDist(i) = 21;
            end
            
        end

        % Check whether new clusters pass muster
        if any(IsolDist < 20 | Lratio > 0.4)
            
            if sum((IsolDist < 20 | Lratio > 0.4)) < 2
                
                nclID = idxF;
                metrics.IsDs = IsolDist;
                metrics.LRat = Lratio;
                
                oCens = computeCentroids(idxF, normalX);
                
                clustREFINE = false;
                
            elseif sum((IsolDist < 20 | Lratio > 0.4)) >= 2
                
                % Fix idxf
                presNums = unique(idxF);
                newNumsC = 1:1:length(presNums);
                for ppi = 1:length(presNums)
                    
                    changeIND = idxF == presNums(ppi);
                    idxF(changeIND) = newNumsC(ppi);
                    
                end
                clustREFINE = true;
                
            end
        else
            
            nclID = idxF;
            metrics.IsDs = IsolDist;
            metrics.LRat = Lratio;
            
            oCens = computeCentroids(idxF, normalX);
            
            clustREFINE = false;
            
        end
        
    else
        
        nclID = idxF;
        metrics.IsDs = IsolDist;
        metrics.LRat = Lratio;
        
        oCens = computeCentroids(idxF, normalX);
        
        clustREFINE = false;

    end
    
    
end




end



function newCens = computeCentroids(idx, featData)

idxIDS = unique(idx);
newCens = zeros(numel(unique(idx)),size(featData,2));
for fi = 1:numel(unique(idx))

    tempID = idxIDS(fi);
    tempFeats = featData(idx == tempID,:);
    newCens(fi,:) = mean(tempFeats);
    
end

end
