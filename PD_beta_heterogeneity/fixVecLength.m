function [newSTRUCT] = fixVecLength(filename)

% load('MPR2_L_11112022.mat')

load(filename)

fiNames = {'MS','ML','BS','BL'};
subNames = {'S1', 'S2', 'S3'};
newSTRUCT = outStruct;
for fi = 1:length(fiNames)

    for si = 1:length(subNames)
        tmpF = outStruct.(fiNames{fi}).(subNames{si});

        newFs = fieldnames(tmpF);

        for ni = 1:length(newFs)

            tmpFs = tmpF.(newFs{ni})(1:13400);

            newSTRUCT.(fiNames{fi}).(subNames{si}).(newFs{ni}) = tmpFs;
        end
    end

end


%% do infos