function final_MaxDur = calcMaxSpkDuration(All_SpikesPerMove_Tbl)

maxDurAll = zeros(height(All_SpikesPerMove_Tbl), 1);

for row_i = 1:height(All_SpikesPerMove_Tbl)
    tempSpk = All_SpikesPerMove_Tbl.C1{row_i};
    if isempty(tempSpk) || numel(tempSpk) < 3
        maxDurAll(row_i) = NaN;
    else
        maxDurAll(row_i) = max(tempSpk - All_SpikesPerMove_Tbl.TTL_spk_idx_Start(row_i));
    end


end

final_MaxDur = max(maxDurAll); % sample num

% MaxDur_sec = final_MaxDur/AO_spike_fs;

spikesMatrix = zeros(height(All_SpikesPerMove_Tbl), final_MaxDur);

for row_i = 1:height(All_SpikesPerMove_Tbl)
    tempSpk = All_SpikesPerMove_Tbl.C1{row_i};
    if isempty(tempSpk) || numel(tempSpk) < 3
        continue;
    else
        spike_indices = tempSpk - All_SpikesPerMove_Tbl.TTL_spk_idx_Start(row_i);
        spikesMatrix(row_i, spike_indices) = 1;
    end

end

spikesMatrix = spikesMatrix(143:156,:);
[nTrials, T] = size(spikesMatrix); % T (time in samp)
bin = 4400; % 10 ms
M = floor(T/bin);
c = reshape(sum(spikesMatrix(:,1:M*bin),1), bin, M);
counts_bin = sum(c,1);
psth_bin_Hz = (counts_bin / nTrials) * (1000/bin);
time_bin_ms = (0:M-1)*bin + bin/2;

figure;
subplot(2,1,1)
[row, col] = find(spikesMatrix);   % row = trial index, col = time (ms)
scatter(col, row, 8, 'r', 'filled');

figure;
plot(time_bin_ms,psth_bin_Hz,'K','LineWidth',3)


end