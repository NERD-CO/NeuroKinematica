function [ expData ] = extractWaveforms_Clz_v06(spikeData, spkThresh, noisThresh, minDist, Fs , pnFlag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Initial Vars
sampLen = round(Fs/1000);
% boundary = sampLen*2;

% Positive waves
[pWavesLocs ,  ~] = peakseek(spikeData , minDist , spkThresh);

% Negative waves
ridP = spikeData;
ridP(ridP > 0) = 0;
[nWavesLocs ,  ~] = peakseek(abs(ridP) , minDist , spkThresh);

% Positive Noise
if sum(spikeData > noisThresh) == 0
    pNoiseLocs = [];
else
    [pNoiseLocs ,  ~] = peakseek(spikeData , 1 , noisThresh);
end

pNoiseInds = ~ismember(pWavesLocs,pNoiseLocs);
pWavesLocs = pWavesLocs(pNoiseInds);
pWavesLocs = pWavesLocs(pWavesLocs > sampLen);

% Negative Noise
if sum(abs(ridP) > noisThresh) == 0
    nNoiseLocs = [];
else
    [nNoiseLocs ,  ~] = peakseek(spikeData , 1 , noisThresh);
end

nNoiseInds = ~ismember(nWavesLocs,nNoiseLocs);
nWavesLocs = nWavesLocs(nNoiseInds);
nWavesLocs = nWavesLocs(nWavesLocs > sampLen);



[posWaves , posWpost] = computeWaves(1,spikeData,pWavesLocs,sampLen);
[negWaves , negWpost] = computeWaves(1,spikeData,nWavesLocs,sampLen);

% allWaves = [negWaves , posWaves];
if pnFlag == 1
    allWaves = posWaves;
    allLocs = pWavesLocs;
    allPositions = posWpost;
    expData.posWaveInfo.waves = posWaves;
    expData.posWaveInfo.positions = posWpost;
    expData.posWaveInfo.numWaves = size(posWaves,2);
    expData.posWaveInfo.posLocs = pWavesLocs;
else
    allWaves = negWaves;
    allLocs = nWavesLocs;
    allPositions = negWpost;
    expData.negWaveInfo.waves = negWaves;
    expData.negWaveInfo.positions = negWpost;
    expData.negWaveInfo.numWaves = size(negWaves,2);
    expData.negWaveInfo.neglocs = nWavesLocs;
end

expData.waveAnylsInfo.spkThresh = spkThresh;
expData.waveAnylsInfo.noiseThresh = noisThresh;

expData.allWaves.waves = allWaves;
expData.allWaves.positions = allPositions;
expData.allWavesInfo.numWaves = size(allWaves,2);
expData.allWavesInfo.alllocs = allLocs;







end




% Subfunction
function [WavesOut , PositionOut] = computeWaves(nOp,spkData,waves,sampLen)

totLen = length(1 - (sampLen/2):1 + (round(sampLen*0.8)-1));


WavesOut = zeros(totLen, length(waves));
PositionOut = zeros(totLen, length(waves));
for ai = 1:length(waves)
    
    locP = waves(ai);
    if nOp == -1
        sampPeriod = locP - (round(sampLen*0.8)-1):locP + (sampLen/2);
    elseif nOp == 1
        sampPeriod = locP - (sampLen/2):locP + (round(sampLen*0.8)-1);
    end

    
    if any(sampPeriod > length(spkData))
        
        WavesOut(:,ai) = [];
        PositionOut(:,ai) = [];
        
        break
    else

        WavesOut(:,ai) = spkData(sampPeriod);
        PositionOut(:,ai) = sampPeriod;

    end
    
end



end