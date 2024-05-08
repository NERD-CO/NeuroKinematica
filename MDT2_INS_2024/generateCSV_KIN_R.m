
matList = dir("*.mat");
matLIST2 = {matList.name};


for mi = 1:length(matLIST2)

    load(matLIST2{mi},'outDATAFin')
    filePARTS = split(matLIST2{mi},{'_','.'});
    subID = filePARTS{1};
    hemiID = filePARTS{2};

    groupIDs = {'GroupA','GroupB'};
    for grI = 1:2

        tmpCSVs = outDATAFin.lfpKINdata.(groupIDs{grI}).moveXY;
        tmpMOVEa = outDATAFin.lfpKINdata.(groupIDs{grI}).moveID;

        for ttI = 1:length(tmpCSVs)

            tmpCSVi = tmpCSVs{ttI};
            tmpMOVEi = tmpMOVEa{ttI};

            % generate timeVector = 60fps

            signalLength = height(tmpCSVi);
            ts = transpose(0:1/60:((signalLength-1)/60));
            
            % create table
            outTABLE = array2table([ts tmpCSVi],'VariableNames',{'ts','x','y'});


            % Save table with TS, X, and Y
            % Name with Subject, hemisphere, Group, Movement

            tmpNAME = [subID , '_' , hemiID,'_',groupIDs{grI},'_',...
                tmpMOVEi,'.csv'];
            writetable(outTABLE,tmpNAME)

        end
    end




end
