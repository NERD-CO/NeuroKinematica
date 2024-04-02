
function [B2 , L2] = getRidOfsmall_NN(B,L)

threshS = 60;
keepIDlog = cellfun(@(x) round(numel(x)/2), B, "UniformOutput",true);
B2 = B(keepIDlog > threshS);
ridIDnum = find(keepIDlog < threshS);

L2 = L;
L2(ismember(L,ridIDnum)) = 0;

end