
cd('C:\Users\Admin\Documents\Github\NeuroKinematica\RM_Capstone')
locTable = readtable('depthLOCS.xlsx');

dirLetter = 'F';

merLOC = [dirLetter,':\Ruby_M_Capstone\MER_Data'];
dlcLOC = [dirLetter,':\Ruby_M_Capstone\DLC_Data'];

frAll = zeros(1000, 1);
betPall = zeros(1000, 1);
cluCount = 1;
subID = zeros(1000,1);
moveID = cell(1000,1);
distID = cell(1000,1);
allmovePS = zeros(1000,1);
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

    movePS = (dlctmptab.Timepoints(end) - dlctmptab.Timepoints(1))/height(dlctmptab);

    % loop through clusters
    for cluI = 1:cluUni

        cluSind = spikeClInfo.clusterIDS == cluUids(cluI);
        cluFrame = cluSind(spikeIndF:spikeIndL);

        % 1. Extract 1s between FRAME indicies
        frameD = spikeClInfo.SpikeTSindex(spikeIndF:spikeIndL);
        spikesInFrame = frameD(cluFrame);

        % 2. Compute total duration
        if isempty(spikesInFrame)
            frAll(cluCount) = NaN;
            betPall(cluCount) = NaN;
            subID(cluCount) = locTable.subID(fi);
            moveID{cluCount} = locTable.MoveT{fi};
            distID{cluCount} = locTable.dist{fi};
            allmovePS(cluCount) = movePS;
            cluCount = cluCount + 1;
            continue

        else

            totDur = (spikesInFrame(end) - spikesInFrame(1))/44000;

            if totDur < 1
                frAll(cluCount) = NaN;
                betPall(cluCount) = NaN;
                subID(cluCount) = locTable.subID(fi);
                moveID{cluCount} = locTable.MoveT{fi};
                distID{cluCount} = locTable.dist{fi};
                allmovePS(cluCount) = movePS;
                cluCount = cluCount + 1;
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
                fft_result = fft(spike_rate);
                n = length(spike_rate); % number of points
                % frequencies = (0:n-1)*(1/totDur)/n;
                frequencies = (0:n-1) / n * (1/bin_size);

                half_n = ceil(n/2);
                positive_freqs = frequencies(1:half_n);
                power_spectrum = abs(fft_result(1:half_n)).^2;

                beta_indices = find(positive_freqs >= 13 & positive_freqs <= 30);
                beta_frequencies = positive_freqs(beta_indices);
                beta_power = power_spectrum(beta_indices);

                if isempty(beta_power)
                    betPall(cluCount) = NaN;
                else
                    betPall(cluCount) = mean(beta_power);
                end

                subID(cluCount) = locTable.subID(fi);
                moveID{cluCount} = locTable.MoveT{fi};
                distID{cluCount} = locTable.dist{fi};
                allmovePS(cluCount) = movePS;
                cluCount = cluCount + 1;

            end
        end
    end
end

% create table

outTable = table(subID,moveID,distID,frAll,betPall,allmovePS,...
    'VariableNames',{'Subject','Movement','STN_Loc','FiringRate','BetaPower',...
    'MoveSpeedSecs'});


cutOFFtab = find(frAll == 0,1,'first');

finalTable = outTable(1:cutOFFtab - 1,:);

writetable(finalTable,'RM_finalData.csv')









