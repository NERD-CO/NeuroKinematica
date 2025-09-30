function [] = plot_Raster_PSTH(CaseDate, All_SpikesPerMove_Tbl, AO_spike_fs, Case_FRKin_Dir)
% plot_Raster_PSTH per Move Rep

%% depth and movement ids

depth_ids = {'t','c','b'};
depth_labels = {'dorsal STN','central STN','ventral STN'};
depth_colors = [0 0.6 0; 0.6 0 0.8; 0 0.4 0.8];
rest_color = [0.5 0.5 0.5];

sorted_moveTypes = {'HAND OC','HAND PS','ARM EF','REST'};
active_movements = {'HAND OC','HAND PS','ARM EF'};
move_types = intersect(sorted_moveTypes, unique(All_SpikesPerMove_Tbl.MoveType),'stable');
hasREST = any(strcmp(move_types,'REST'));


%% Main Function

for moveT_i = 1:numel(move_types)
    move_type = move_types{moveT_i};
    for depth_i = 1:numel(depth_ids)
        depth_code = depth_ids{depth_i};
        depth_name = depth_labels{depth_i};

        % Extract the relevant trials for this MoveType and depth
        move_tbl = All_SpikesPerMove_Tbl(strcmp(All_SpikesPerMove_Tbl.MoveType, move_type) & ...
            contains(All_SpikesPerMove_Tbl.move_trial_ID, depth_code), :);
        if isempty(move_tbl), continue; end

        % Calculate the max spike segment duration (# of samples in longest segment)
        max_SpkDuration_All = zeros(height(move_tbl), 1);

        % Calculate the max spike segment duration (# of samples in longest segment)
        % max_SpkDuration_All = zeros(height(All_SpikesPerMove_Tbl), 1);

        for row_i = 1:height(move_tbl)
            tempSpk_vec = move_tbl.C1{row_i};
            if isempty(tempSpk_vec) || numel(tempSpk_vec) < 3
                max_SpkDuration_All(row_i) = NaN;
            else
                max_SpkDuration_All(row_i) = max(tempSpk_vec - move_tbl.TTL_spk_idx_Start(row_i));
            end
        end

        Max_SpikeDuration_samples = max(max_SpkDuration_All); % samples in longest segment
        Max_SpkDur_seconds = Max_SpikeDuration_samples/AO_spike_fs; % seconds
        Max_SpkDus_ms = Max_SpkDur_seconds * 1000; % milliseconds


        spikesMatrix = zeros(height(move_tbl), Max_SpikeDuration_samples);
        for row_i = 1:height(move_tbl)
            tempSpk_vec = move_tbl.C1{row_i};
            if isempty(tempSpk_vec) || numel(tempSpk_vec) < 3
                continue;
            else
                spike_indices = tempSpk_vec - move_tbl.TTL_spk_idx_Start(row_i);
                spikesMatrix(row_i, spike_indices) = 1;
            end

        end

        
        % Calculate the peri-stimulus time histogram (PSTH) and prepare for plotting
        binSize = 10; % 10 ms
        binSize_sec = binSize ./ 1000; % 0.01
        binSize_samp = AO_spike_fs * binSize_sec; % 440

        [nTrials, Time_samples] = size(spikesMatrix); % T (time in samp)
        M = floor(Time_samples/binSize_samp); % samples

        spk_counts_per_sample = sum(spikesMatrix(:,1:M*binSize_samp), 1); % total across trials
        spk_counts = reshape(spk_counts_per_sample, binSize_samp, M);
        counts_bin = sum(spk_counts,1);  % spikes per bin across all trials

        %psth_bin_Hz = (counts_bin / nTrials) * (1000/bin_samp);
        % psth_bin_Hz = (counts_bin / nTrials) / binSize_sec;
        psth_bin_Hz = (counts_bin / nTrials) * (binSize_samp/binSize);

        time_bin_samp = (0:M-1) * binSize_samp + (binSize_samp/2); % bin centers (samples)
        time_bin_ms = (time_bin_samp / AO_spike_fs) * 1000;   % bin centers (ms)

        %fig = figure('Position',[100 100 800 600]);

        % Plot the scatter of spikes and the PSTH in the same figure
        [row, col] = find(spikesMatrix);      % row = trial index, col = spike sample index
        col_ms = (col / AO_spike_fs) * 1000;  % convert to ms (col = spike time)
        figure;
        subplot(2,1,1);
        scatter(col_ms, row, 8, 'b', 'filled');
        xlabel('Time (ms)');
        ylabel('Trial Index');
        title('Spike Raster Plot', sprintf('%s | %s ', depth_name, move_type));
        % xlim([0, max(col_ms)]);
        xlim([0, 400]);
        hold on;

        % Peri-Stimulus Time Histogram
        subplot(2,1,2);
        plot(time_bin_ms,psth_bin_Hz,'k','LineWidth',3)
        xlabel('Time (ms)');
        ylabel('PSTH (Hz)');
        title('Peri-Stimulus Time Histogram', sprintf('%s | %s ', depth_name, move_type));
        % xlim([0, max(time_bin_ms)]);
        xlim([0, 400]);
        grid on;

    end
end


%% Save figure
% 
% save_filename = fullfile(Case_FRKin_Dir, [CaseDate '_Raster_PSTH.png']);
% saveas(fig, save_filename);
% fprintf('Saved figure to: %s\n', save_filename);
