function [ outfeats ] = createSpkFeatures2( waveforms, spkFS, normFlag )
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
[pVal, pInd] = max(tempWaves,[],2);
[vVal, vInd] = min(tempWaves,[],2);

features.widthMS = (abs(pInd-vInd)/round(spkFS))*1000;
features.Amp = pVal - vVal;

% features.WavePC1 = pcScores(:,1);

features.FSDE_Values =...
    FSDE_Method_d2(tempWaves);

intialX = [features.Peak , abs(features.Valley) ,...
    features.Energy , features.widthMS,...
    features.FSDE_Values.FDmax , features.FSDE_Values.SDmax ,...
    features.Amp];

% Eliminate negative values

if normFlag
    
    ridNegX = abs(intialX);
    
    % Normalize to 0-1
    outfeats = bsxfun(@rdivide, ridNegX, max(ridNegX));
    
else
    
    outfeats = intialX;
    
end


end

