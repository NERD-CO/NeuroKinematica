function [ttlCount,ttlLog,ttlMlist] = matfileTTLcheck(matdir)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

cd(matdir)

matdirlist1 = dir('*.mat');
matdirlist2 = {matdirlist1.name};

ttlLog = zeros(length(matdirlist2),1,'logical');
ttlMlist = cell(length(matdirlist2),1);
ttlCount = zeros(length(matdirlist2),1);
for mi = 1:length(matdirlist2)

    tmpName = matdirlist2{mi};
    tmpMatfile = matfile(tmpName);
    tmpVarName = whos(tmpMatfile);
    VarNameCk = {tmpVarName.name};
    ttlCHECK = matches('CDIG_IN_1_KHz',VarNameCk);
    if any(ttlCHECK)
        ttlLog(mi) = true;
        ttlMlist{mi} = tmpName;
        load(tmpName,"CDIG_IN_1_Down")
        ttlCount(mi) = length(CDIG_IN_1_Down);
    end
end

ttlMlist = ttlMlist(cellfun(@(x) ~isempty(x), ttlMlist, 'UniformOutput',true));
ttlCount = ttlCount(ttlLog);


end