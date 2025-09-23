%% MaxSpkDuration_Raster_PSTH

% Sampling rates
TTL_fs = 44000;         % Hz, Alpha Omega TTL clock
AO_spike_fs = 44000;    % Hz, Alpha Omega spike sampling rate
AO_LFP_fs = 1375;       % Hz, Alpha Omega LFP sampling rate
DLC_fs = 100;           % fps, Video/DLC frame rate

% Patient specific file dirs
CaseDate = '03_23_2023'; % Adjust as needed
% '03_23_2023'; % NER 2025
% '04_05_2023'; % NER 2025
% '05_18_2023_b_bilateral'; % NER 2025
% '05_31_2023';
% '06_08_2023_bilateral'; % NER 2025
% '08_23_2023'; % NANS 2026

MoveDir_CaseID = 'IO_03_23_2023_LSTN'; % Adjust as needed
% 'IO_03_23_2023_LSTN'; % NER 2025
% 'IO_04_05_2023_RSTN'; % NER 2025
% 'IO_05_18_2023_b_LSTN'; % NER 2025
% 'IO_05_31_2023_LSTN';
% 'IO_06_08_2023_LSTN'; % NER 2025
% 'IO_08_23_2023_RSTN'; % NANS 2026


% Data folder paths
Case_DataDir = fullfile(IO_DataDir, CaseDate);
% ephysTbl_Dir = fullfile(Case_DataDir,'DLC_Ephys');
MoveDataDir = fullfile(IO_DataDir, 'Processed DLC');
KinematicsDir = fullfile(IO_DataDir, 'Kinematic Analyses');
Ephys_Kin_Dir = fullfile(IO_DataDir, 'Ephys_Kinematics');
FR_Kin_Dir = fullfile(Ephys_Kin_Dir, 'FR_Kinematic_Analyses');
Case_FRKin_Dir = fullfile(FR_Kin_Dir, CaseDate); % input dir (new ephysTbleDir) and output (results) dir


%% Handle Bilateral Cases

% Specify hemisphere in command window if bilateral
isBilateral = contains(CaseDate,'bilateral','IgnoreCase',true);
if isBilateral
    fprintf('[INFO] Bilateral case detected: %s\n',CaseDate);
    CaseDate_hem = input('Enter hemisphere (LSTN or RSTN): ','s');
    if ~ismember(CaseDate_hem,{'LSTN','RSTN'})
        error('Invalid hemisphere');
    end
    % Append hemisphere-specific folder
    Case_FRKin_Dir = fullfile(Case_FRKin_Dir, CaseDate_hem);
else
    CaseDate_hem = '';
end
fprintf('[INFO] Loading spike data from: %s\n', Case_FRKin_Dir);


%% Load All_SpikesPerMove_Tbl

cd(Case_FRKin_Dir)
Tbl_list = dir('*Spikes*.mat');
Tbl_names = {Tbl_list.name};
spk_case = Tbl_names{contains(Tbl_names, 'offset')}; % offset version preferred
load(spk_case, 'All_SpikesPerMove_Tbl');


%% OPTIONAL: Case-specific cleaning (remove duplicates if needed)

% for 3_23_2023 case - remove duplicates / only plot for primary electrode
if CaseDate == '03_23_2023'
    All_SpikesPerMove_Tbl = All_SpikesPerMove_Tbl(168:end,1:13); % Comment or adjust as needed
else
    All_SpikesPerMove_Tbl = All_SpikesPerMove_Tbl;
end


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
        xlim([0, max(col_ms)]);
        hold on;

        % Peri-Stimulus Time Histogram
        subplot(2,1,2);
        plot(time_bin_ms,psth_bin_Hz,'k','LineWidth',3)
        xlabel('Time (ms)');
        ylabel('PSTH (Hz)');
        title('Peri-Stimulus Time Histogram', sprintf('%s | %s ', depth_name, move_type));
        xlim([0, max(time_bin_ms)]);
        grid on;

    end
end


%% Save figure
% 
% save_filename = fullfile(Case_FRKin_Dir, [CaseDate '_Raster_PSTH.png']);
% saveas(fig, save_filename);
% fprintf('Saved figure to: %s\n', save_filename);
