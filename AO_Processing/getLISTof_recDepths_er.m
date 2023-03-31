function [] = getLISTof_recDepths_er(studyID,matLOC, saveLOC)

cd(matLOC)
% cd('C:\Users\Admin\Downloads\Raw Electrophysiology AO\Raw Electrophysiology Mat')
matDir = dir('*.mat');
matDIR2 = {matDir.name};


fulldepthName = cell(30,1);
trimdepthName = cell(30,1);
numdepthName = cell(30,1);
framesTTLc = nan(30,1);
trajectorY = nan(30,1);
hemisphere = cell(30,1);
countC = 1;
% Loop through files
for mi = 1:length(matDIR2)
    
    matTemp = matDIR2{mi};

    tempHemi = matTemp(1);
    tempTraj = str2double(matTemp(3));

    disp(['Assessing file: ', matTemp])

    depthName1 = extractAfter(matTemp,4);
    depthName2 = extractBefore(depthName1,length(depthName1)-3);
    % Look for TTLs at depth
    deciLoc = strfind(depthName2,'.');
    mmDepthN = round(str2double(depthName2(1:deciLoc+2)),2);

    if mmDepthN < 7
        disp(['File ', matTemp ,' was below 7mm from target'])
        % Look for TTLs 
        matftemp1 = whos(matfile(matTemp));
        matVarList = {matftemp1.name};

        ttlCHECK = matches('CDIG_IN_1_KHz',matVarList);

        if ttlCHECK
            disp(['File ', matTemp ,' had TTLs'])
            % Look for TTLs for a certain duration
            load(matTemp,'CDIG_IN_1_Up')

            % check for 2 seconds of recording
            if numel(CDIG_IN_1_Up) > 120
                disp(['File ', matTemp ,' had 120 frames - SAVING'])

                fulldepthName{countC} = matTemp;
                trimdepthName{countC} = depthName2;
                numdepthName{countC} = mmDepthN;
                framesTTLc(countC) = numel(CDIG_IN_1_Up);
                trajectorY(countC) = tempTraj;
                hemisphere{countC} = tempHemi;

                countC = countC + 1;
            else
                disp(['File ', matTemp ,' had ', num2str(numel(CDIG_IN_1_Up)) ,...
                    ' frames - NOT SAVING'])
            end
        else
            disp(['File ', matTemp ,' did not have TTLs'])
        end
    end
end

% Create table
outTable1 = table(fulldepthName,trimdepthName,numdepthName,framesTTLc,trajectorY,...
    hemisphere,'VariableNames',{'FullFile','DepthName','DepthMM','FrameCount',...
    'TrajNum','Hemi'});
outTableF = outTable1(1:find(isnan(outTable1.FrameCount),1,'first')-1,:);

cd(saveLOC);

saveNAME = [studyID , '_DLClist.mat'];

save(saveNAME,"outTableF");




end



