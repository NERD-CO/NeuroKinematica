function [ outfeats ] = createSpkFeatures( waveforms, spkFS )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

tempWaves = transpose(waveforms);

features = struct;

% Peak
features.Peak = max(tempWaves,[],2);
% Valley
features.Valley = min(tempWaves,[],2);
% Energy
features.Energy = trapz(abs(tempWaves),2);
% Combine for WavePCA analysis

% Peak to Peak width
[~, pInd] = max(tempWaves,[],2);
[~, vInd] = min(tempWaves,[],2);

features.widthMS = (abs(pInd-vInd)/round(spkFS))*1000;

% features.WavePC1 = pcScores(:,1);

features.FSDE_Values =...
    FSDE_Method(tempWaves);

intialX = [features.Peak , abs(features.Valley) ,...
    features.Energy , features.widthMS,...
    features.FSDE_Values.FDmin , features.FSDE_Values.SDmax ,...
    features.FSDE_Values.SDmin];

% Eliminate negative values
ridNegX = abs(intialX);

% Normalize to 0-1
outfeats = bsxfun(@rdivide, ridNegX, max(ridNegX));



end

