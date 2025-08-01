function [] = IO_plotting_v3(case_ID, plot_ID, ephys_offset)

clear; clc;
addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis'

% comps_IO_plotting(CaseDate, 'raster_fig', 1)

%% Directory Setup

% specify directory where case-specific data files are located
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


%% Inputs:

plot_ID = 'raster_fig';
ephys_offset = 1;

%% Isolate specific CaseDate / studyID (StudyNum in Subject_AO csv)

% CaseDate = '03_09_2023'; % studyID = 1, ptID 1

% CaseDate = '03_23_2023'; % studyID = 2, ptID 2    *
% CaseDate = '04_05_2023'; % studyID = 3, ptID 2    *

% CaseDate = '04_13_2023_bilateral'; % studyID = 4(L*), 5(R), ptID 3

% CaseDate = '05_11_2023'; % studyID = 6, ptID 4
% CaseDate = '05_18_2023_a'; % studyID = 7, ptID 4

 CaseDate = '05_18_2023_b_bilateral'; % studyID = 8(L*), 9(R), ptID = 5

% CaseDate = '05_31_2023';  % studyID = 10, ptID 6

% CaseDate = '06_08_2023_bilateral'; % studyID = 11(L*), 12(R), ptID = 7

% CaseDate = '07_13_2023_bilateral'; % studyID = 15(L), 16(R), ptID = 9

case_ID = CaseDate;

%% Define Case-specific Data Directories

Case_DataDir = [IO_DataDir, filesep, CaseDate];
ephysTbl_Dir = [Case_DataDir, filesep, 'DLC_Ephys'];

% Handle bilateral cases and hemisphere selection
isBilateral = contains(CaseDate, 'bilateral', 'IgnoreCase', true);

if isBilateral
    fprintf('[INFO] Bilateral case detected: %s\n', CaseDate);
    
    % Prompt user for hemisphere choice (LSTN or RSTN)
    CaseDate_hem = input('Enter hemisphere (LSTN or RSTN): ', 's');
    
    % Validate input
    validHems = {'LSTN','RSTN'};
    if ~ismember(CaseDate_hem, validHems)
        error('Invalid input. Please enter LSTN or RSTN.');
    end
    
    % Append hemisphere-specific folder
    ephysTbl_Dir = fullfile(ephysTbl_Dir, CaseDate_hem);
    fprintf('[INFO] Hemisphere-specific directory set: %s\n', ephysTbl_Dir);
else
    CaseDate_hem = ''; % No hemisphere for unilateral cases
    fprintf('[INFO] Using base processed directory: %s\n', ephysTbl_Dir);
end


%% Define sampling rates (Alpha Omega and Video sampling fs)

TTL_fs = 44000; % Hz
AO_spike_fs = 44000; % Hz
AO_LFP_fs = 1375; % Hz
DLC_fs = 100; % fps


%% Main function

addpath 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\AO_Processing'

% define FR parameters (spikes/sec)
binSize_FR = 0.01; % 10 ms bin size
window_FR = [-0.05 0.45]; % 50ms before to 450ms after movement onset
edges_FR = window_FR(1):binSize_FR:window_FR(2); % bin edges for 10ms bins

% define PSTH parameters
binSize_PSTH = 0.05; % 50 ms bin size
window_PSTH = [-0.05 0.45]; % 50ms before to 450ms after movement onset
edges_PSTH = window_PSTH(1):binSize_PSTH:window_PSTH(2); % bin edges for 50ms bins

% switch case setup
switch plot_ID
    case 'raster_fig'   % spike response during movement reps
        % Navigate to ephys directory
        cd(ephysTbl_Dir)
        Tbl_list = dir('*.mat');
        Tbl_names = {Tbl_list.name};

        % Filter for spike tables
        if ephys_offset         % offset = 50ms before trial start
            search_index = contains(Tbl_names, 'Spikes') & contains(Tbl_names, 'offset');
        else
            search_index = contains(Tbl_names, 'Spikes') & ~contains(Tbl_names, 'offset');
        end

        if ~any(search_index)
            error('No matching files found in ephys directory.');
        end
        spk_case = Tbl_names{search_index};
        load(spk_case, 'All_SpikesPerMove_Tbl');


        %%% Case-specific modifications %%%
        % for 3_23_2023 case only: Remove duplicates / only plot for primary electrode  %%%% comment out / adjust for other cases
        % if CaseDate == '03_23_2023'
        %     All_SpikesPerMove_Tbl = All_SpikesPerMove_Tbl(165:end,1:13);  % 3_23_2023: Remove duplicates / only plot for primary electrode
        % else
        %     All_SpikesPerMove_Tbl = All_SpikesPerMove_Tbl;
        % end
         

        % Automatically extract depth-specific tables sub-indexed by movement type
        depth_specific_tables = extract_STN_depth_tables(All_SpikesPerMove_Tbl);

        % Determine count and unique ID(s) of Motor task reps
        moveType_num = numel(unique(All_SpikesPerMove_Tbl.MoveType));
        moveType_ids = unique(All_SpikesPerMove_Tbl.MoveType);

        % Plot FR per MovRep for each MovType at all STN depths
        figure;
        tiledlayout(1, moveType_num)
        for fig_i = 1:moveType_num
            move_subtbl = All_SpikesPerMove_Tbl(matches(All_SpikesPerMove_Tbl.MoveType, moveType_ids{fig_i}),:);

            % Determine count and unique ID(s) of STN depths in All_SpikesPerMove_Tbl
            depth_num = cellfun(@(x) x(1), All_SpikesPerMove_Tbl.move_trial_ID, 'UniformOutput', false);
            depth_ids = unique(depth_num);

            % uisetcolor, coloors.com
            color_vec = [0    0.4471    0.7412; ...
                0.7176    0.2745    1.0000; ...
                0.4667    0.6745    0.1882];

            % % initialize legend info and color mapping
            % legend_info = {};
            % depth_color_idx = [];
            % legend_handles = []; % gobjects(0) % empty array of graphic objects

            nexttile
            % loop through MovReps for each MovType per STN depth
            for movRep_row = 1:height(move_subtbl)
                temp_row = move_subtbl(movRep_row,:);
                temp_spks = temp_row.C1{1};
                depth_idx = find(matches(depth_ids,temp_row.move_trial_ID{1}(1)));
                temp_seconds = transpose(((temp_spks - temp_row.TTL_spk_idx_Start)/AO_spike_fs) - 0.05); % incorporate this in All_SpikesPerMove_Tbl
                m_line = line([temp_seconds; temp_seconds], ...
                    [zeros(1,numel(temp_seconds)) + movRep_row; ...
                    zeros(1,numel(temp_seconds)) + movRep_row+1], ...
                    'color', color_vec(depth_idx,:));

                % store legend info, only if depth hasn't been plotted yet
                % if ~ismember(depth_idx, depth_color_idx)
                %     legend_info{end+1} = sprintf('depth %s', depth_ids{depth_idx});
                %     depth_color_idx(end+1) = depth_idx;
                %     depth_color = color_vec(depth_idx,:);
                %     legend_handles = [legend_handles; m_line];
                % end
            end

            title(moveType_ids{fig_i})
            xlabel('Time (s)')
            ylabel('Movement Repetition')
            xlim([-0.1 0.4])
            % ylim([0 50])

            % resize and move figure(s) to stay within bounds of display screen
            set(gcf, 'Units', 'normalized', 'OuterPosition', [0.2, 0.2, 0.6, 0.6]);
        end


        PSTH_all_trials = [];
        FR_all_trials = [];

        % loop through MovReps for each MovType at all STN depths
        for movT_i = 1:moveType_num
            moveType = moveType_ids{movT_i};
            move_subtbl = All_SpikesPerMove_Tbl(matches(All_SpikesPerMove_Tbl.MoveType, moveType_ids{movT_i}),:);

            % Determine count and unique ID(s) of STN depths in All_SpikesPerMove_Tbl
            depth_num = cellfun(@(x) x(1), All_SpikesPerMove_Tbl.move_trial_ID, 'UniformOutput', false);
            depth_ids = unique(depth_num);

            % Determine count of STN depth trial reps in move_subtbl
            depth_num2 = cellfun(@(x) x(1), move_subtbl.move_trial_ID, 'UniformOutput', false);

            % loop through each STN depth, compute FR per Mov rep at each STN depth trial in move_subtbl
            for depth_i = 1:length(depth_ids)
                % Filter table for current depth
                temp_trialID = [move_subtbl.move_trial_ID{depth_i},'_',move_subtbl.MoveType{depth_i},'_',num2str(move_subtbl.MoveN(depth_i))];

                temp_depth_idx = matches(depth_num2, depth_ids{depth_i}); % logical idx
                move_depth_tbl = move_subtbl(temp_depth_idx,:);

                temp_depth_row = move_subtbl(depth_i,:);
                temp_depth_spks = temp_depth_row.C1{1};
                temp_depth_seconds = transpose(((temp_depth_spks - temp_depth_row.TTL_spk_idx_Start)/AO_spike_fs) - 0.05); % incorporate this in All_SpikesPerMove_Tbl

                % figure;
                % histogram(cell2mat(move_depth_tbl.C1_ts))
                figure;
                nbins = 10;
                BinWidth = 0.05;
                BinLimits = [-0.05 0.45];
                histogram(temp_depth_seconds,nbins,"BinWidth",0.05,"BinLimits",[-0.05 0.45])
                % need y-axis to show tick values as the number of spikes
                % per 50ms bin (converted to number of spikes per second)
                % save rasters and PSTHs as PDFs - make editale in ADOBE

                % resize and move figure(s) to stay within bounds of display screen
                set(gcf, 'Units', 'normalized', 'OuterPosition', [0.2, 0.2, 0.6, 0.6]);

                % Extract spike times for current trial
                spikeTimes = temp_depth_seconds;

                % Create histogram of spike times
                spike_counts_PSTH = histcounts(spikeTimes, edges_PSTH);
                spike_counts_FR = histcounts(spikeTimes, edges_FR);

                % Store PSTH for this trials
                PSTH_all_trials = [PSTH_all_trials; spike_counts_PSTH];

                % Convert spike counts to FR (spikes per second)
                FR_moverep = spike_counts_FR/binSize_FR;

                % Store FR for move trial
                FR_all_trials = [FR_all_trials; FR_moverep];
            end

            % Average PSTH across all trials
            mean_PSTH = mean(PSTH_all_trials, 1);
            std_PSTH = std(PSTH_all_trials, 0, 1);

            % Apply Gaussian smoothing for visualization
            smoothed_PSTH = smoothdata(mean_PSTH, 'gaussian', 5);

            % Plot the PSTH
            figure;
            bar(edges_PSTH(1:end-1), smoothed_PSTH);
            xlabel('Time (s) from movement onset');
            ylabel('Spike Count (per 50ms bin)');
            title('Peri-Stimulus Time Histogram (PSTH) - Smoothed');

            % resize and move figure(s) to stay within bounds of display screen
            set(gcf, 'Units', 'normalized', 'OuterPosition', [0.2, 0.2, 0.6, 0.6]);

            % Compute mean and standard deviation firing rate across all trials
            mean_FR = mean(FR_all_trials, 1);
            std_FR = std(FR_all_trials, 0, 1);

            % Plot the average firing rate across trials with error bars
            figure;
            bar(edges_FR(1:end-1), mean_FR);
            xlabel('Time (s) from movement onset');
            ylabel('Firing Rate (spikes/s)');
            title('Average Firing Rate (10ms bins)');

            % resize and move figure(s) to stay within bounds of display screen
            set(gcf, 'Units', 'normalized', 'OuterPosition', [0.2, 0.2, 0.6, 0.6]);

            % Print output
            % fprintf('Movement type (all depths) %s: Mean FR = %f\n', moveType_ids{movT_i}, mean(mean_FR));

        end

end


%% Programatically index STN depth-specific data from All_SpikesPerMove_Tbl

depth_specific_tables = extract_STN_depth_tables(All_SpikesPerMove_Tbl);

% Extract subtables for specific depths and movement types
dorsalSTN_tbl = depth_specific_tables.t;  % 't' corresponds to dorsal STN
centralSTN_tbl = depth_specific_tables.c; % 'c' corresponds to central STN
ventralSTN_tbl = depth_specific_tables.b; % 'b' corresponds to ventral STN

% Test - Access specific movement types for a depth
dorsal_HAND_OC_tbl = depth_specific_tables.t.HANDOC;
central_ARM_EF_tbl = depth_specific_tables.c.ARMEF;
ventral_HAND_PS_tbl = depth_specific_tables.b.HANDPS;

% Iterate over each depth and movement type
for depth_field = fieldnames(depth_specific_tables)'
    depth = depth_field{1}; % Extract depth key ('t', 'c', or 'b')
    depth_tbl = depth_specific_tables.(depth);

    for move_field = fieldnames(depth_tbl)'
        move_type = move_field{1};
        move_tbl = depth_tbl.(move_type);

        fprintf('Processing Depth: %s, Movement Type: %s\n', depth, move_type);

        % Compute FR and PSTH for each movement at each depth
        FR_all_trials = [];
        depth_num3 = cellfun(@(x) x(1), move_tbl.move_trial_ID, 'UniformOutput', false);
        depth_ids = unique(depth_num3);

        for depth_i = 1:length(depth_ids)
            temp_depth_idx = matches(depth_num3, depth_ids{depth_i});
            move_depth_tbl = move_tbl(temp_depth_idx, :);

            for row_i = 1:height(move_depth_tbl)
                temp_depth_row = move_depth_tbl(row_i,:);
                temp_depth_spks = temp_depth_row.C1{1};

                if isempty(temp_depth_spks)
                    warning('No spike data for trial %d at depth %s. Skipping...', row_i, depth_ids{depth_i});
                    continue;
                end

                % Convert spike times to seconds
                temp_depth_seconds = transpose(((temp_depth_spks - temp_depth_row.TTL_spk_idx_Start) / AO_spike_fs) - 0.05);

                % Compute firing rate
                spikeCounts_FR = histcounts(temp_depth_seconds, edges_FR);
                FR_moverep = spikeCounts_FR / binSize_FR;

                % Store FR for movement repetition
                FR_all_trials = [FR_all_trials; FR_moverep];
            end
        end

        % Compute mean and std FR for the current movement type
        mean_FR = mean(FR_all_trials, 1);
        std_FR = std(FR_all_trials, 0, 1);

        % check / compare to old method
        FR_perMoveRep = transpose(mean_FR);
        % Initialize results structure
        if ~exist('results', 'var') || isempty(results)
            results = struct();
        end

        % Check if the field exists before trying to access it
        if isfield(results, depth) && isfield(results.(depth), move_type) && isfield(results.(depth).(move_type), 'FR_perMoveRep') ...
                && ~isempty(results.(depth).(move_type).FR_perMoveRep)
            mean_FR_perDepth = mean(vertcat(results.(depth).(move_type).FR_perMoveRep), 1);
        else
            % warning('Skipping mean_FR_perDepth computation: No valid data for Depth: %s, Movement Type: %s', depth, move_type);
            mean_FR_perDepth = NaN; % Assign NaN to prevent errors
        end
        mean_FR_perMoveRep = mean(mean_FR_perDepth);
        [s_FR, m_FR] = std(FR_perMoveRep);

        % Save results in a struct
        results.(depth).(matlab.lang.makeValidName(move_type)).mean_FR = mean_FR;
        results.(depth).(matlab.lang.makeValidName(move_type)).std_FR = std_FR;
        results.(depth).(matlab.lang.makeValidName(move_type)).FR_perMoveRep = FR_perMoveRep;
        results.(depth).(matlab.lang.makeValidName(move_type)).mean_FR_perMoveRep = mean_FR_perMoveRep;
        results.(depth).(matlab.lang.makeValidName(move_type)).std_FR_perMoveRep = s_FR;

        % Print summary in the command window
        fprintf('Depth: %s, Movement Type: %s -> Mean FR = %.2f spikes/s, Std FR = %.2f spikes/s\n', ...
            depth, move_type, mean(mean_FR), mean(std_FR));

        % % Plot FR histogram
        % figure;
        % bar(edges_FR(1:end-1), mean_FR);
        % xlabel('Time (s) from movement onset');
        % ylabel('Firing Rate (spikes/s)');
        % title(sprintf('FR for %s at %s STN', move_type, depth));
        % resize and move figure(s) to stay within bounds of display screen
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0.2, 0.2, 0.6, 0.6]);

    end
end

%% Incorporate baseline comparisons
% percent change (within subject)
% normalize - z-score (across subject)
    % based on baseline FR avg.
    % average FR across trials per movement type
 % template: https://elifesciences.org/articles/64893
 % Fig 2, c, d
 % Fig 4
 % Fig 7

 
 % 1) Create separate function for sub-function to change ways to extract /
 % align (e.g. to start of move, peak of mov, etc.) - as input

% 2) New function:
% read in AllSPkperMov table, add 2 new cols
% 1. sum # of spikes / rep
% 2. latency b/t last start time of move rep trial and time of 1st spike in
% move rep trial

% save using '-append'


end

%% Sub-function(s)

% Automate STN depth-specific move_subtbl extraction:

function depth_specific_tables = extract_STN_depth_tables(All_SpikesPerMove_Tbl)

% Extract unique STN depths from the 'move_trial_ID' column
depth_ids = unique(cellfun(@(x) x(1), All_SpikesPerMove_Tbl.move_trial_ID, 'UniformOutput', false));

% Extract unique movement types
move_types = unique(All_SpikesPerMove_Tbl.MoveType);

% Initialize a struct to store subtables
depth_specific_tables = struct();

% Iterate over each depth
for depth_idx = 1:length(depth_ids)
    depth = depth_ids{depth_idx};

    % Filter rows for the current depth
    depth_rows = contains(All_SpikesPerMove_Tbl.move_trial_ID, depth);
    depth_tbl = All_SpikesPerMove_Tbl(depth_rows, :);

    % Initialize a nested struct for this depth
    depth_specific_tables.(depth) = struct();

    % Iterate over each movement type
    for move_idx = 1:length(move_types)
        move_type = move_types{move_idx};

        % Filter rows for the current movement type
        move_rows = matches(depth_tbl.MoveType, move_type);
        move_tbl = depth_tbl(move_rows, :);

        % Store the subtable in the nested struct
        depth_specific_tables.(depth).(matlab.lang.makeValidName(move_type)) = move_tbl;
    end
end
end