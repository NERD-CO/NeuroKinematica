function [corInds] = computeASCI_Clz_v05(waveforms)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% if numel(unique(idx)) == 1;
%      corInds = nan(length(idx),1);
%      return
% end


tempWaves = transpose(waveforms);

corVals = zeros(size(tempWaves,1),1);

% Create Rp, Rn and Rz from clust template

meanCw = mean(tempWaves);

template = meanCw/max(abs(meanCw));

% get start number of waveforms
startNumWaves = size(tempWaves,1);
acceptReduc = round(startNumWaves*0.85);

numCheck = true;

crit = 0.625;

while numCheck
    
    % Start While loop
    
    Rp = nan(size(template));
    Rn = nan(size(template));
    
    for mi = 1:length(template)
        
        tRz = template(mi);
        
        if abs(tRz) <= 0.25
            Rp(mi) = 0.25;
            Rn(mi) = -0.25;
        elseif tRz > 0.25
            Rp(mi) = max([tRz/2, 0.25]);
            Rn(mi) = Rp(mi) - 0.25;
        else
            Rn(mi) = min([tRz/2 , -0.25]);
            Rp(mi) = Rn(mi) + 0.25;
        end
        
    end
    
    
    for ai = 1:startNumWaves
        
        tWave = tempWaves(ai,:);
        tWaveN = tWave/max(abs(meanCw));
        
        SY = nan(1,size(tWaveN,2));
        SX = template;
        
        for ti = 1:size(tWaveN,2)
            
            yi = tWaveN(ti);
            if yi <= Rp(ti) && yi >= Rn(ti)
                SY(ti) = 0;
            elseif yi >= Rp(ti)
                SY(ti) = 1;
            elseif yi <= Rn(ti)
                SY(ti) = -1;
            else
                test = 1;
            end
            
            xi = template(ti);
            if xi <= Rp(ti) && xi >= Rn(ti)
                SX(ti) = 0;
            elseif xi >= Rp(ti)
                SX(ti) = 1;
            elseif xi <= Rn(ti)
                SX(ti) = -1;
            end
            
        end
        
        %------------------------------------------------------------
        % To satisfy dot prod makes indicies where both zeros equal 1
        bzInd = SY == 0 & SX == 0;
        SX(bzInd) = 1;
        SY(bzInd) = 1;
        %------------------------------------------------------------
        
        SCC = sum(SY.*SX)/length(SY);
        
        corVals(ai) = SCC;
        
    end

    corInds = false(size(corVals));
    
    corInds(corVals >= crit) = true;
    
    endNumWaves = sum(corInds);
    
    if endNumWaves < acceptReduc
        if crit >= 0.575
            numCheck = true;
            crit = crit - 0.025;
        else
            numCheck = false;
        end
    else
        numCheck = false;
    end
    
end
% end while loop







end

