function TEST__IO_FR_Compute_and_Plot(case_ID, plot_ID, ephys_offset, IO_DataDir)

% comps_IO_plotting('raster_fig', '03_23_2023', 1, 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative')

% case_ID = 
% multiHem = ; % 0 (unilat), 1 (bilat)
plot_ID = 'raster_fig';
ephys_offset = 1; % 0 (no offset), 1 (offset)


%% specify directory where case-specific data files are located
curPCname = getenv('COMPUTERNAME');

switch curPCname
    case 'DESKTOP-I5CPDO7'  % PC_1
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'DSKTP-JTLAB-EMR'  % Lab Desktop
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'NSG-M-H8J3X34'    % PC_2
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
end

cd(IO_DataDir)
Subject_AO = readtable('Subject_AO.xlsx');

%% Specify the case-specific data directory
Case_DataDir = fullfile(IO_DataDir, case_ID);
ephysTbl_Dir = fullfile(Case_DataDir, 'DLC_Ephys');

% % For Bialteral cases, specify hemisphere:
CaseDate_hem = 'LSTN'; % comment out when N/A
% CaseDate_hem = 'RSTN'; % comment out when N/A

% % For Bialteral cases:
ephysTbl_Dir = [ephysTbl_Dir, filesep, CaseDate_hem]; % comment out when N/A

%% Define sampling rates

TTL_fs = 44000; % Hz
AO_spike_fs = 44000; % Hz
AO_LFP_fs = 1375; % Hz
DLC_fs = 100; % fps

%% Main function

% Add GitHub library path
addpath 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\AO_Processing';

% PSTH and Firing Rate parameters
binSize_PSTH = 0.05; % 50 ms bin size
window_PSTH = [-0.05 0.45]; % 50ms before to 450ms after movement onset
edges_PSTH = window_PSTH(1):binSize_PSTH:window_PSTH(2); % bin edges for 50ms bins

binSize_FR = 0.01; % 10 ms bin size
window_FR = [-0.05 0.45];
edges_FR = window_FR(1):binSize_FR:window_FR(2);

% Process data based on the specified plot_ID
switch plot_ID
    case 'raster_fig'
        % Navigate to ephys directory
        cd(ephysTbl_Dir);
        Tbl_list = dir('*.mat');
        Tbl_names = {Tbl_list.name};

        % Filter for spike tables
        if ephys_offset
            search_index = contains(Tbl_names, 'Spikes') & contains(Tbl_names, 'offset');
        else
            search_index = contains(Tbl_names, 'Spikes') & ~contains(Tbl_names, 'offset');
        end

        if ~any(search_index)
            error('No matching files found in ephys directory.');
        end

        spk_case = Tbl_names{search_index};
        load(spk_case, 'All_SpikesPerMove_Tbl');

        % Process movement types and depths
        moveType_ids = unique(All_SpikesPerMove_Tbl.MoveType);
        results = struct();

        % Process each movement type
        for movT_i = 1:numel(moveType_ids)
            moveType = moveType_ids{movT_i};
            move_subtbl = All_SpikesPerMove_Tbl(matches(All_SpikesPerMove_Tbl.MoveType, moveType), :);

            % Skip empty tables
            if isempty(move_subtbl)
                warning('No data available for movement type %s. Skipping...', moveType);
                continue;
            end

            % Process each STN depth
            depth_num = cellfun(@(x) x(1), move_subtbl.move_trial_ID, 'UniformOutput', false);
            depth_ids = unique(depth_num);

            % Initialize storage for PSTH and FR results
            PSTH_all_trials = [];
            FR_all_trials = [];

            for depth_i = 1:length(depth_ids)
                % Filter table for current depth
                temp_depth_idx = matches(depth_num, depth_ids{depth_i});
                move_depth_tbl = move_subtbl(temp_depth_idx, :);

                % Skip if depth table is empty
                if isempty(move_depth_tbl)
                    warning('No data available for depth %s. Skipping...', depth_ids{depth_i});
                    continue;
                end

                % Process each trial for the current depth
                for row_i = 1:height(move_depth_tbl)
                    temp_depth_row = move_depth_tbl(row_i, :);
                    temp_depth_spks = temp_depth_row.C1{1};

                    % Skip empty spike data
                    if isempty(temp_depth_spks)
                        warning('No spike data for trial %d. Skipping...', row_i);
                        continue;
                    end

                    % Convert spike times to seconds
                    temp_depth_seconds = transpose(((temp_depth_spks - temp_depth_row.TTL_spk_idx_Start) / AO_spike_fs) - 0.05);

                    % Compute PSTH and FR
                    spike_counts_PSTH = histcounts(temp_depth_seconds, edges_PSTH);
                    spike_counts_FR = histcounts(temp_depth_seconds, edges_FR);

                    PSTH_all_trials = [PSTH_all_trials; spike_counts_PSTH];
                    FR_all_trials = [FR_all_trials; spike_counts_FR / binSize_FR];
                end
            end

            % Compute mean and standard deviation of PSTH
            if ~isempty(PSTH_all_trials)
                mean_PSTH = mean(PSTH_all_trials, 1);
                std_PSTH = std(PSTH_all_trials, 0, 1);
            else
                mean_PSTH = [];
                std_PSTH = [];
                warning('No PSTH data for movement type %s.', moveType);
            end

            % Compute mean and standard deviation of firing rates
            if ~isempty(FR_all_trials)
                mean_FR = mean(FR_all_trials, 1);
                std_FR = std(FR_all_trials, 0, 1);
            else
                mean_FR = [];
                std_FR = [];
                warning('No firing rate data for movement type %s.', moveType);
            end


            % Plot results
            if ~isempty(mean_PSTH)
                figure;
                bar(edges_PSTH(1:end-1), mean_PSTH, 'FaceColor', 'b');
                title(sprintf('PSTH for Movement Type %s', moveType));
                xlabel('Time (s)');
                ylabel('Spike Count');
            end

            if ~isempty(mean_FR)
                figure;
                bar(edges_FR(1:end-1), mean_FR, 'FaceColor', 'r');
                title(sprintf('Firing Rate for Movement Type %s', moveType));
                xlabel('Time (s)');
                ylabel('Firing Rate (spikes/s)');
            end
        end
end
end
