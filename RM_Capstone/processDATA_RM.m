%% 03_09_2023


%%

dLOCmove = 'Z:\Martini_R_2024\DLC_Analysis\IO_2023-03-09_P1_RSTN\processedMovement';
cd(dLOCmove)

dlcDepths = {'0p4','2p0'};

csvList1 = dir('*.csv');
csvList2 = {csvList1.name};
for dlcId = 1:length(dlcDepths)

    dlcDepL = transpose(csvList2(contains(csvList2,dlcDepths{dlcId})));

end

%%

dLOCmer = 'Z:\Martini_R_2024\DLC_MER\DLC_Mat2\03_09_2023\ClusteredSpikeTimes';
cd(dLOCmer)

spkDEpths = {'0.453','2.059'};

cd(dLOCmer)
matList1 = dir('*.mat');
matList2 = {matList1.name};
for spkId = 1:length(spkDEpths)

    spkDepL = transpose(matList2(contains(matList2,spkDEpths{spkId})));

end

%% 05 11 2023

%%

dLOCmove = 'Z:\Martini_R_2024\DLC_Analysis\IO_2023-05-11_P4_LSTN\processedMovement';
cd(dLOCmove)

dlcDepths = {'3p7','5p6'};

csvList1 = dir('*.csv');
csvList2 = {csvList1.name};
dlcDepL = {};
for dlcId = 1:length(dlcDepths)

    dlcDepLtmp = transpose(csvList2(contains(csvList2,dlcDepths{dlcId})));
    dlcDepL = [dlcDepL ; dlcDepLtmp]

end

%%

dLOCmer = 'Z:\Martini_R_2024\DLC_MER\DLC_Mat2\05_11_2023\ClusteredSpikeTimes';
cd(dLOCmer)

spkDEpths = {'3.701','5.640'};

cd(dLOCmer)
matList1 = dir('*.mat');
matList2 = {matList1.name};
spkDepL = {};
for spkId = 1:length(spkDEpths)

    spkDepLtmp = transpose(matList2(contains(matList2,spkDEpths{spkId})));
    spkDepL = [spkDepL ; spkDepLtmp]

end

%% 06_08_2023_LSTN

%%

dLOCmove = 'Z:\Martini_R_2024\DLC_Analysis\IO_06_08_2023_LSTN\processedMovement';
cd(dLOCmove)

dlcDepths = {'0p75','1p4','3p5'};

csvList1 = dir('*.csv');
csvList2 = {csvList1.name};
dlcDepL = {};
for dlcId = 1:length(dlcDepths)

    dlcDepLtmp = transpose(csvList2(contains(csvList2,dlcDepths{dlcId})));
    dlcDepL = [dlcDepL ; dlcDepLtmp]

end

%%

dLOCmer = 'Z:\Martini_R_2024\DLC_MER\DLC_Mat2\06_08_2023_L\ClusteredSpikeTimes';
cd(dLOCmer)

spkDEpths = {'0.757','1.402','3.504'};

cd(dLOCmer)
matList1 = dir('*.mat');
matList2 = {matList1.name};
spkDepL = {};
for spkId = 1:length(spkDEpths)

    spkDepLtmp = transpose(matList2(contains(matList2,spkDEpths{spkId})));
    spkDepL = [spkDepL ; spkDepLtmp]

end


%% 06_08_2023_RSTN

%%

dLOCmove = 'Z:\Martini_R_2024\DLC_Analysis\IO_06_08_2023_RSTN\processedMovement';
cd(dLOCmove)

dlcDepths = {'0p3','1p9'};

csvList1 = dir('*.csv');
csvList2 = {csvList1.name};
dlcDepL = {};
for dlcId = 1:length(dlcDepths)

    dlcDepLtmp = transpose(csvList2(contains(csvList2,dlcDepths{dlcId})));
    dlcDepL = [dlcDepL ; dlcDepLtmp]

end

%%

dLOCmer = 'Z:\Martini_R_2024\DLC_MER\DLC_Mat2\06_08_2023_R\ClusteredSpikeTimes';
cd(dLOCmer)

spkDEpths = {'0.324','1.927'};

cd(dLOCmer)
matList1 = dir('*.mat');
matList2 = {matList1.name};
spkDepL = {};
for spkId = 1:length(spkDEpths)

    spkDepLtmp = transpose(matList2(contains(matList2,spkDEpths{spkId})));
    spkDepL = [spkDepL ; spkDepLtmp]

end























%%

% load table
cd('C:\Users\Admin\Documents\Github\NeuroKinematica\RM_Capstone');
dataTable = readtable("depthLOCS.xlsx");

for ddi = 1:height(dataTable)

    spkD = dataTable.spkDepth{ddi};
    dlcD = dataTable.dlcDepth{ddi};

    cleanSpk = spkD(2:end-1);
    cleanDlc = dlcD(2:end-1);

    cd(dLOCmove)
    tmpDlctab = readtable(cleanDlc);

    cd(dLOCmer)
    load(cleanSpk,'spikeClInfo')

    % Do something


end

%%

cd('C:\Users\Admin\Documents\Github\NeuroKinematica\RM_Capstone')
locTable = readtable('depthLOCS.xlsx');


merLOC = 'I:\Ruby_M_Capstone\MER_Data';
dlcLOC = 'I:\Ruby_M_Capstone\DLC_Data';



















