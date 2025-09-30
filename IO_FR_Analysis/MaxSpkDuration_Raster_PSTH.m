function [Max_SpikeDuration_samples, spikesMatrix] = MaxSpkDuration_Raster_PSTH(CaseDate, All_SpikesPerMove_Tbl, AO_spike_fs, Case_FRKin_Dir)

% Calculate the max spike segment duration (# of samples in longest segment)
max_SpkDuration_All = zeros(height(All_SpikesPerMove_Tbl), 1);

for row_i = 1:height(All_SpikesPerMove_Tbl)
    tempSpk_vec = All_SpikesPerMove_Tbl.C1{row_i};
    if isempty(tempSpk_vec) || numel(tempSpk_vec) < 3
        max_SpkDuration_All(row_i) = NaN;
    else
        max_SpkDuration_All(row_i) = max(tempSpk_vec - All_SpikesPerMove_Tbl.TTL_spk_idx_Start(row_i));
    end
end

Max_SpikeDuration_samples = max(max_SpkDuration_All); % samples in longest segment
% Max_SpkDur_seconds = Max_SpikeDuration_samples/AO_spike_fs; % seconds
% Max_SpkDus_ms = Max_SpkDur_seconds * 1000; % milliseconds


spikesMatrix = zeros(height(All_SpikesPerMove_Tbl), Max_SpikeDuration_samples);
for row_i = 1:height(All_SpikesPerMove_Tbl)
    tempSpk_vec = All_SpikesPerMove_Tbl.C1{row_i};
    if isempty(tempSpk_vec) || numel(tempSpk_vec) < 3
        continue;
    else
        spike_indices = tempSpk_vec - All_SpikesPerMove_Tbl.TTL_spk_idx_Start(row_i);
        spikesMatrix(row_i, spike_indices) = 1;
    end

end


% Calculate the peri-stimulus time histogram (PSTH) and prepare for plotting
binSize_ms = 10; % 10 ms 
bin_sec = binSize_ms ./ 1000; % 0.01
bin_samp = AO_spike_fs * bin_sec; % 440

[nTrials, Time_samples] = size(spikesMatrix); % T (time in samp)
M = floor(Time_samples/bin_samp); % samples

spk_counts_per_sample = sum(spikesMatrix(:,1:M*bin_samp), 1); % total across trials
spk_counts = reshape(spk_counts_per_sample, bin_samp, M);
counts_bin = sum(spk_counts,1);  % spikes per bin across all trials

% psth_bin_Hz = (counts_bin / nTrials) * (1000/bin_samp);
psth_bin_Hz = (counts_bin / nTrials) / bin_sec;

time_bin_samp = (0:M-1) * bin_samp + (bin_samp/2); % bin centers (samples)
time_bin_ms = (time_bin_samp / AO_spike_fs) * 1000;   % bin centers (ms)

% fig = figure('Position',[100 100 800 600]);

% Plot the scatter of spikes and the PSTH in the same figure
[row, col] = find(spikesMatrix);      % row = trial index, col = spike sample index
col_ms = (col / AO_spike_fs) * 1000;  % convert to ms (col = spike time)
figure; 
subplot(2,1,1);
scatter(col_ms, row, 8, 'b', 'filled');
xlabel('Time (ms)');
ylabel('Trial Index');
title('Spike Raster Plot');
xlim([0, max(col_ms)]);
hold on;

% Peri-Stimulus Time Histogram
subplot(2,1,2);
plot(time_bin_ms,psth_bin_Hz,'k','LineWidth',3)
xlabel('Time (ms)');
ylabel('PSTH (Hz)');
title('Peri-Stimulus Time Histogram');
xlim([0, max(time_bin_ms)]);
grid on;


%% Save figure

% save_filename = fullfile(Case_FRKin_Dir, [CaseDate '_Raster_PSTH.png']);
% saveas(fig, save_filename);
% fprintf('Saved figure to: %s\n', save_filename);


end