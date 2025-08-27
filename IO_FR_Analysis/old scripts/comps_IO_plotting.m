function [] = comps_IO_plotting(plot_ID, case_ID, ephys_offset)

% comps_IO_plotting('raster_fig', CaseDate, 1)
plot_ID = 'raster_fig';
ephys_offset = 1;

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


%% Isolate a specific CaseDate / studyID (StudyNum in Subject_AO csv)

% CaseDate = '03_09_2023'; % studyID = 1, ptID 1

% CaseDate = '03_23_2023'; % studyID = 2, ptID 2    *
% CaseDate = '04_05_2023'; % studyID = 3, ptID 2    *

% CaseDate = '04_13_2023_bilateral'; % studyID = 4(L*), 5(R), ptID 3

% CaseDate = '05_11_2023'; % studyID = 6, ptID 4
% CaseDate = '05_18_2023_a'; % studyID = 7, ptID 4

% CaseDate = '05_18_2023_b_bilateral'; % studyID = 8(L*), 9(R), ptID = 5

% CaseDate = '05_31_2023';  % studyID = 10, ptID 6

CaseDate = '06_08_2023_bilateral'; % studyID = 11(L*), 12(R), ptID = 7

% CaseDate = '07_13_2023_bilateral'; % studyID = 15(L), 16(R), ptID = 9

case_ID = CaseDate;

%% Define Case-specific Data Directories

Case_DataDir = [IO_DataDir, filesep, CaseDate];
ephysTbl_Dir = [Case_DataDir, filesep, 'DLC_Ephys'];

%% For Bialteral cases, specify hemisphere + corresponding dir

CaseDate_hem = 'LSTN'; % comment out when N/A
% CaseDate_hem = 'RSTN'; % comment out when N/A

 ephysTbl_Dir = [ephysTbl_Dir, filesep, CaseDate_hem]; % comment out when N/A


%% Define sampling rates (Alpha Omega and Video sampling fs)

TTL_fs = 44000; % Hz
AO_spike_fs = 44000; % Hz
AO_LFP_fs = 1375; % Hz
DLC_fs = 100; % fps


%% Main function

addpath 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\AO_Processing'

% switch case setup            
switch plot_ID
    case 'raster_fig'   % spike response during movement reps
        cd(ephysTbl_Dir)
        Tbl_list = dir('*.mat');
        Tbl_names = {Tbl_list.name};

        if ephys_offset     % offset = 50ms before trial start
            search_index = contains(Tbl_names, 'Spikes') & contains(Tbl_names, 'offset');
        else
            search_index = contains(Tbl_names, 'Spikes') & ~contains(Tbl_names, 'offset');
        end
        spk_case = Tbl_names{search_index};
        load(spk_case, 'All_SpikesPerMove_Tbl');

        % % for 3_23_2023 case - remove duplicates / only plot for primary electrode  %%%% comment out / adjust for other cases
        % All_SpikesPerMove_Tbl = All_SpikesPerMove_Tbl(165:end,1:13);

        % determine count and unique ID(s) of Motor task reps
        moveType_num = numel(unique(All_SpikesPerMove_Tbl.MoveType));
        moveType_ids = unique(All_SpikesPerMove_Tbl.MoveType);

        % plot for all 3 depths
        figure;
        tiledlayout(1,moveType_num)
        for fig_i = 1:moveType_num
            move_subtbl = All_SpikesPerMove_Tbl(matches(All_SpikesPerMove_Tbl.MoveType, moveType_ids{fig_i}),:);

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
            % loop through depths
            for m_row = 1:height(move_subtbl)
                temp_row = move_subtbl(m_row,:);
                temp_spks = temp_row.C1{1};
                depth_idx = find(matches(depth_ids,temp_row.move_trial_ID{1}(1)));
                temp_seconds = transpose(((temp_spks - temp_row.TTL_spk_idx_Start)/AO_spike_fs) - 0.05); % incorporate this in All_SpikesPerMove_Tbl
                m_line = line([temp_seconds; temp_seconds], ...
                    [zeros(1,numel(temp_seconds)) + m_row; ...
                    zeros(1,numel(temp_seconds)) + m_row+1], ...
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
            xlabel('time (s)')
            ylabel('movement repetition')
            xlim([-0.1 0.4])
            % ylim([0 50])

        end


        % case 'PSTH_fig'
        binSize_PSTH = 0.05; % 50 ms bin size for PSTH
        window_PSTH = [-0.05 0.45]; % 50ms before to 450ms after movement onset
        edges_PSTH = window_PSTH(1):binSize_PSTH:window_PSTH(2); % bin edges for 50ms bins
        PSTH_all_trials = [];

        % case 'FR_fig'     % firing rate (spike/s)
        binSize_FR = 0.01; % 10 ms bin size
        window_FR = [-0.05 0.45]; % 50ms before to 450ms after movement onset
        edges_FR = window_FR(1):binSize_FR:window_FR(2); % bin edges for 10ms bins
        FR_all_trials = [];

        
        % loop through movement types
        for fig_i = 1:moveType_num
            move_subtbl = All_SpikesPerMove_Tbl(matches(All_SpikesPerMove_Tbl.MoveType, moveType_ids{fig_i}),:);

            depth_num = cellfun(@(x) x(1), All_SpikesPerMove_Tbl.move_trial_ID, 'UniformOutput', false);
            depth_ids = unique(depth_num);

            % loop through depths
            depth_num2 = cellfun(@(x) x(1), move_subtbl.move_trial_ID, 'UniformOutput', false);

            for depth_i = 1:length(depth_ids)
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

            % Compute mean and standard deviation firing rate across all trials
            mean_FR = mean(FR_all_trials, 1);
            std_FR = std(FR_all_trials, 0, 1);

            % Plot the average firing rate across trials with error bars
            figure;
            bar(edges_FR(1:end-1), mean_FR);
            xlabel('Time (s) from movement onset');
            ylabel('Firing Rate (spikes/s)');
            title('Average Firing Rate (10ms bins)');

        end

end

%% Manually index STN depth-specific data from All_SpikesPerMove_Tbl

%%% automate programatically later %%%

% 03/23/23 case
% dorsalSTN_tbl = All_SpikesPerMove_Tbl(115:end,1:9);
% centralSTN_tbl = All_SpikesPerMove_Tbl(57:114,1:9);
% ventralSTN_tbl = All_SpikesPerMove_Tbl(1:56,1:9);

% % 04/05/23 case   *
% dorsalSTN_tbl = [All_SpikesPerMove_Tbl(17:23,1:9); All_SpikesPerMove_Tbl(40:46,1:9)]; % move_trial_ID contains 't'
% centralSTN_tbl = [All_SpikesPerMove_Tbl(10:16,1:9); All_SpikesPerMove_Tbl(33:39,1:9)]; % move_trial_ID contains 'c'
% ventralSTN_tbl = [All_SpikesPerMove_Tbl(1:9,1:9); All_SpikesPerMove_Tbl(24:32,1:9)]; % move_trial_ID contains 'b'

% % 04/13/23, LSTN case     *
% dorsalSTN_tbl = [All_SpikesPerMove_Tbl(41:end,1:9)]; % move_trial_ID contains 't'
% centralSTN_tbl = [All_SpikesPerMove_Tbl(21:40,1:9)]; % move_trial_ID contains 'c'
% ventralSTN_tbl = [All_SpikesPerMove_Tbl(1:20,1:9)]; % move_trial_ID contains 'b'

% % 05/18/23 (b), LSTN case     *
% dorsalSTN_tbl = [All_SpikesPerMove_Tbl(47:72,1:9)]; % move_trial_ID contains 't'
% centralSTN_tbl = [All_SpikesPerMove_Tbl(25:46,1:9)]; % move_trial_ID contains 'c'
% ventralSTN_tbl = [All_SpikesPerMove_Tbl(1:24,1:9)]; % move_trial_ID contains 'b'

% 06/08/23, LSTN case
dorsalSTN_tbl = All_SpikesPerMove_Tbl(73:end,1:9); % move_trial_ID contains 't'(singleTrial)
centralSTN_tbl = All_SpikesPerMove_Tbl(20:38,1:9); % move_trial_ID contains 'c'(singleTrial)
ventralSTN_tbl = All_SpikesPerMove_Tbl(1:19,1:9); % move_trial_ID contains 'b' (singleTrial)

% % 07_13_2023, LSTN case
% dorsalSTN_tbl = [All_SpikesPerMove_Tbl(38:end,1:9)]; % move_trial_ID contains 't'
% centralSTN_tbl = [All_SpikesPerMove_Tbl(15:37,1:9)]; % move_trial_ID contains 'c'
% ventralSTN_tbl = [All_SpikesPerMove_Tbl(1:14,1:9)]; % move_trial_ID contains 'b'

% % % 07_13_2023, RSTN case     *
% dorsalSTN_tbl = [All_SpikesPerMove_Tbl(17:end,1:9)]; % move_trial_ID contains 't'
% centralSTN_tbl = [All_SpikesPerMove_Tbl(9:16,1:9)]; % move_trial_ID contains 'c'
% ventralSTN_tbl = [All_SpikesPerMove_Tbl(1:8,1:9)]; % move_trial_ID contains 'b'


%% dorsal STN

move_subtbl = dorsalSTN_tbl;

% run code within loop above (starting with depth_num2)
depth_num2 = cellfun(@(x) x(1), move_subtbl.move_trial_ID, 'UniformOutput', false);
for depth_i = 1:length(depth_ids)
    temp_trialID = [move_subtbl.move_trial_ID{depth_i},'_',move_subtbl.MoveType{depth_i},'_',num2str(move_subtbl.MoveN(depth_i))];

    temp_depth_row = move_subtbl(depth_i,:);
    temp_depth_spks = temp_depth_row.C1{1};
    temp_depth_seconds = transpose(((temp_depth_spks - temp_depth_row.TTL_spk_idx_Start)/AO_spike_fs) - 0.05); % incorporate this in All_SpikesPerMove_Tbl

    % Extract spike times for current trial
    spikeTimes = temp_depth_seconds;

    % Create histogram of spike times
    spike_counts_FR = histcounts(spikeTimes, edges_FR);

    % Convert spike counts to FR (spikes per second)
    FR_moverep = spike_counts_FR/binSize_FR;

    % Store FR for move trial
    FR_all_trials = [FR_all_trials; FR_moverep];
end

% Compute mean and standard deviation firing rate across all trials
d_mean_FR = mean(FR_all_trials, 1);
d_std_FR = std(FR_all_trials, 0, 1);

dorsal_FR_perMoveRep = transpose(d_mean_FR);
mean_dorsalFR = mean(dorsal_FR_perMoveRep);
[s_dFR, m_dFR] = std(dorsal_FR_perMoveRep)


%% central STN

move_subtbl = centralSTN_tbl;

% run code within loop above (starting with depth_num2)
depth_num2 = cellfun(@(x) x(1), move_subtbl.move_trial_ID, 'UniformOutput', false);
for depth_i = 1:length(depth_ids)
    temp_trialID = [move_subtbl.move_trial_ID{depth_i},'_',move_subtbl.MoveType{depth_i},'_',num2str(move_subtbl.MoveN(depth_i))];

    temp_depth_row = move_subtbl(depth_i,:);
    temp_depth_spks = temp_depth_row.C1{1};
    temp_depth_seconds = transpose(((temp_depth_spks - temp_depth_row.TTL_spk_idx_Start)/AO_spike_fs) - 0.05); % incorporate this in All_SpikesPerMove_Tbl

    % Extract spike times for current trial
    spikeTimes = temp_depth_seconds;

    % Create histogram of spike times
    spike_counts_FR = histcounts(spikeTimes, edges_FR);

    % Convert spike counts to FR (spikes per second)
    FR_moverep = spike_counts_FR/binSize_FR;

    % Store FR for move trial
    FR_all_trials = [FR_all_trials; FR_moverep];
end

% Compute mean and standard deviation firing rate across all trials
c_mean_FR = mean(FR_all_trials, 1);
c_std_FR = std(FR_all_trials, 0, 1);

central_FR_perMoveRep = transpose(c_mean_FR);
mean_centralFR = mean(central_FR_perMoveRep);
[s_cFR, m_cFR] = std(central_FR_perMoveRep)


%% ventral STN

move_subtbl = ventralSTN_tbl;

% run code within loop above (starting with depth_num2)
depth_num2 = cellfun(@(x) x(1), move_subtbl.move_trial_ID, 'UniformOutput', false);
for depth_i = 1:length(depth_ids)
    temp_trialID = [move_subtbl.move_trial_ID{depth_i},'_',move_subtbl.MoveType{depth_i},'_',num2str(move_subtbl.MoveN(depth_i))];

    temp_depth_row = move_subtbl(depth_i,:);
    temp_depth_spks = temp_depth_row.C1{1};
    temp_depth_seconds = transpose(((temp_depth_spks - temp_depth_row.TTL_spk_idx_Start)/AO_spike_fs) - 0.05); % incorporate this in All_SpikesPerMove_Tbl

    % Extract spike times for current trial
    spikeTimes = temp_depth_seconds;

    % Create histogram of spike times
    spike_counts_FR = histcounts(spikeTimes, edges_FR);

    % Convert spike counts to FR (spikes per second)
    FR_moverep = spike_counts_FR/binSize_FR;

    % Store FR for move trial
    FR_all_trials = [FR_all_trials; FR_moverep];
end

% Compute mean and standard deviation firing rate across all trials
v_mean_FR = mean(FR_all_trials, 1);
v_std_FR = std(FR_all_trials, 0, 1);

ventral_FR_perMoveRep = transpose(v_mean_FR);
mean_ventralFR = mean(ventral_FR_perMoveRep);
[s_vFR, m_vFR] = std(ventral_FR_perMoveRep)


end