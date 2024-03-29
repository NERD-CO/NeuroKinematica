function [] = accelAOcheck(excELloC , allDataLoc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

cd(excELloC)
inTable = readtable("Subject_AO.xlsx");

subjectNUM = unique(inTable.StudyNum);
accFileCheck = zeros(height(inTable),1,'logical');
fCount = 1;
for si = 1:length(subjectNUM)

    tmpSubject = subjectNUM(si);
    subTable = inTable(ismember(inTable.StudyNum,tmpSubject),:);

    fileNAMES = subTable.ao_MAT_file;
    caseDATE = unique(subTable.CaseFolder);
    caseDATE1 = caseDATE{1};
    for fi = 1:length(fileNAMES)

        tmpFile = fileNAMES{fi};
        fileLOC = [allDataLoc , filesep , caseDATE1 , filesep ,...
            'Raw Electrophysiology MATLAB', filesep, tmpFile ];

        accMatfile = matfile(fileLOC);
        accMatList = whos(accMatfile);
        accVarCheck = {accMatList.name};

        accfname = 'CACC_3___01___Sensor_1___X';
        % Is present in Matfile
        firstCheck = ismember(accfname,accVarCheck);

        if ~firstCheck
            fCount = fCount + 1;
            continue
        end

        % Is not empty
        load(fileLOC, accfname);
        secondCheck = isempty(CACC_3___01___Sensor_1___X);
        if secondCheck % is empty
            fCount = fCount + 1;
            continue
        else
            accFileCheck(fCount) = true;
            fCount = fCount + 1;
        end
    end

    disp(['Subject ', num2str(si), ' out of ' num2str(length(subjectNUM)), ' Complete'])
end


cd(excELloC)
saveName = '';
blah = inTable;
blah.IMUCheck = accFileCheck;
save(saveName, "inTable");


end