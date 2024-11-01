function [x_VALUES , y_VALUES] = getCoordDLCframe(inputFRAME, frameNUM)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


colNAMES = inputFRAME.Properties.VariableNames;
colNAMESxy = contains(colNAMES,{'x','y'});
xyCOLSall = inputFRAME(:,colNAMESxy);
xyCOLnames = xyCOLSall.Properties.VariableNames;

frameROW = table2array(xyCOLSall(frameNUM,:));

x_VALUES = transpose(frameROW(contains(xyCOLnames,'x')));

y_VALUES = transpose(frameROW(contains(xyCOLnames,'y')));



end