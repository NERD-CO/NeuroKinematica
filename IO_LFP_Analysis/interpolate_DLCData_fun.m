function dlcDataInterpolated = interpolate_DLCData_fun(dlcData, tsOriginal, tsTarget)
% Function to interpolate DLC data from original timestamps to target LFP timestamps
variableNames = dlcData.Properties.VariableNames;
dlcDataInterpolated = table();
for varIdx = 1:length(variableNames)
    dlcDataInterpolated.(variableNames{varIdx}) = interp1(tsOriginal, dlcData.(variableNames{varIdx}), tsTarget, 'linear', 'extrap');
end
end