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

dirLetter = 'F';

merLOC = [dirLetter,':\Ruby_M_Capstone\MER_Data'];
dlcLOC = [dirLetter,':\Ruby_M_Capstone\DLC_Data'];

frAll = zeros(1000, 1);
betPall = zeros(1000, 1);
cluCount = 1;
for fi = 1:height(locTable)

    % Load in DLC
    tmpFOLDn = locTable.foldN{fi};
    tmpdlcLOC = locTable.dlcDepth{fi};
    tmpdlcLOCc = replace(tmpdlcLOC,'''','');
    dlcLoadt = [dlcLOC , filesep , tmpFOLDn ,filesep, tmpdlcLOCc];
    dlctmptab = readtable(dlcLoadt);

    % Load in MER
    tmpmerLOC = locTable.spkDepth1{fi};
    tmpmerLOCc = replace(tmpmerLOC,'''','');
    merLoadt = [merLOC , filesep , tmpFOLDn ,filesep, tmpmerLOCc];
    load(merLoadt,'spikeClInfo');

    % firing rate during movements
    indiceSfirst = dlctmptab.Locations(1);
    indiceFcon = round((indiceSfirst/120)*44000);
    [~ , spikeIndF] = min(abs(indiceFcon - spikeClInfo.SpikeTSindex));

    indiceSlast = dlctmptab.Locations(height(dlctmptab));
    indiceLcon = round((indiceSlast/120)*44000);
    [~ , spikeIndL] = min(abs(indiceLcon - spikeClInfo.SpikeTSindex));

    % Determine number of clusters
    cluAll = unique(spikeClInfo.clusterIDS);
    cluUni = sum(cluAll < 100000);
    cluUids = cluAll(cluAll < 100000);

    % loop through clusters
    for cluI = 1:cluUni

        cluSind = spikeClInfo.clusterIDS == cluUids(cluI);
        cluFrame = cluSind(spikeIndF:spikeIndL);

        % 1. Extract 1s between FRAME indicies
        frameD = spikeClInfo.SpikeTSindex(spikeIndF:spikeIndL);
        spikesInFrame = frameD(cluFrame);

        % 2. Compute total duration


        totDur = (spikesInFrame(end) - spikesInFrame(1))/44000;

        if totDur < 1
            frAll(cluCount) = NaN;
            betPall(cluCount) = NaN;
            continue
        else

            % 3. Compute Firing rate
            spkFS = numel(spikesInFrame) / totDur;

            % 4. Save Firing rate
            frAll(cluCount) = spkFS;

            % 5. Compute beta freq
            % Convert spike times to spike rate signal
            bin_size = 0.01; % 1 ms bins
            edges = 0:bin_size:totDur; % Create edges for histogram bins
            spikes2secs = spikesInFrame/44000;
            spike_rate = histcounts(spikes2secs, edges);

            % Compute FFT and frequencies
            fft_result = fft(spikes2secs);
            n = length(spikes2secs); % number of points
            frequencies = (0:n-1)*(1/totDur)/n;

            half_n = ceil(n/2);
            positive_freqs = frequencies(1:half_n);
            power_spectrum = abs(fft_result(1:half_n)).^2;

            beta_indices = find(positive_freqs >= 13 & positive_freqs <= 30);
            beta_frequencies = positive_freqs(beta_indices);
            beta_power = power_spectrum(beta_indices);

            if isempty(beta_power)
                betPall(cluCount) = NaN;
            else
                keyboard;
            end

            cluCount = cluCount + 1;

        end




    end









end



















