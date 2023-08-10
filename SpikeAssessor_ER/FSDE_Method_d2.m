function [FSDE_ds] = FSDE_Method_d2(tempWaves)

% FD(n) = s(n) - s(n - 1)
% SD(n) = FD(n) - FD(n - 1)

% s = spike waveform points
% n = n number of points

% Get Index where FDmin,SDmin,SDmax occur

FSDEvals = zeros(size(tempWaves,1),4);
for wavi = 1:length(FSDEvals)
    tempWave = tempWaves(wavi,:);
    detrendW = detrend(tempWave);
    tempFD = firstDer(detrendW);
    tempSD = secondDer(tempFD);
    [~,tempFDMinLoc] = min(tempFD);
    [~,tempSDminLoc] = min(tempSD);
    [~,tempSDmaxLoc] = max(tempSD);
    [~,tempFDmaxLoc] = max(tempFD);
    
    FSDEvals(wavi,1) = detrendW(tempFDMinLoc);
    FSDEvals(wavi,2) = detrendW(tempSDminLoc);
    FSDEvals(wavi,3) = detrendW(tempSDmaxLoc);
    FSDEvals(wavi,4) = detrendW(tempFDmaxLoc);
end

FSDE_ds = mat2dataset(FSDEvals,'VarNames',{'FDmin','SDmin','SDmax','FDmax'});

% hold on
% plot(tempWaves(:,wavi),'.--')
% plot(tempFD,'r.--')
% plot(tempSD,'k.--')
% line([0 32],[0 0],'LineStyle','--','Color','k');


return

% s = tempWaves(:,1);


function FD = firstDer(s)
FD = zeros(length(s),1);
for n = 2:length(s)
    FD(n) = s(n) - s(n - 1);
end

return


function SD = secondDer(s)
SD = zeros(length(s),1);
for n2 = 2:length(s)
    SD(n2) = s(n2) - s(n2 - 1);
end

return

