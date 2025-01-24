%% Compute FR mean and stdev per STN depth in All_SpikesPerMove_Tbl

addpath 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\AO_Processing'

% Specify directory where case-specific data files are located
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


%% Hardcode Case-specific Data directories

% Isolate specific CaseDate / studyID (StudyNum in Subject_AO csv)
% CaseDate = '03_09_2023'; % studyID = 1, ptID 1

% CaseDate = '03_23_2023'; % studyID = 2, ptID 2    *
% CaseDate = '04_05_2023'; % studyID = 3, ptID 2    *

% CaseDate = '04_13_2023_bilateral'; % studyID = 4(L*), 5(R), ptID 3

% CaseDate = '05_11_2023'; % studyID = 6, ptID 4
% CaseDate = '05_18_2023_a'; % studyID = 7, ptID 4

% CaseDate = '05_18_2023_b_bilateral'; % studyID = 8(L*), 9(R), ptID = 5

% CaseDate = '05_31_2023';  % studyID = 10, ptID 6

% CaseDate = '06_08_2023_bilateral'; % studyID = 11(L*), 12(R), ptID = 7

CaseDate = '07_13_2023_bilateral'; % studyID = 15(L), 16(R), ptID = 9

% % For Bialteral cases, specify hemisphere:
% CaseDate_hem = 'LSTN'; % comment out when N/A
CaseDate_hem = 'RSTN'; % comment out when N/A

case_ID = CaseDate;

% define case-specific data directory
Case_DataDir = [IO_DataDir, filesep, CaseDate];
ephysTbl_Dir = [Case_DataDir, filesep, 'DLC_Ephys'];

% % For Bialteral cases:
ephysTbl_Dir = [ephysTbl_Dir, filesep, CaseDate_hem]; % comment out when N/A

%% define datastream sampling rates (Alpha Omega and Video sampling fs)

TTL_fs = 44000; % Hz
AO_spike_fs = 44000; % Hz
AO_LFP_fs = 1375; % Hz
DLC_fs = 100; % fps

%% Load All Spikes per Motor Trial table

cd(ephysTbl_Dir)
Tbl_list = dir('*.mat');
Tbl_names = {Tbl_list.name};

ephys_offset = 1;
if ephys_offset    % offset = 50ms before trial start
    search_index = contains(Tbl_names, 'Spikes') & contains(Tbl_names, 'offset');
else
    search_index = contains(Tbl_names, 'Spikes') & ~contains(Tbl_names, 'offset');
end
spk_case = Tbl_names{search_index};

load(spk_case, 'All_SpikesPerMove_Tbl');


%% Case-specific modifications

% for 3_23_2023 case - remove duplicates / only plot for primary electrode
% All_SpikesPerMove_Tbl = All_SpikesPerMove_Tbl(158:end,1:13); % for 3_23_2023 case %%% comment out / adjust for other cases


%% Compute mean FR for each trial per Motor Task type across STN depths

% define FR bin params % firing rate (spike/s)
binSize_FR = 0.01; % 10 ms bin size
window_FR = [-0.05 0.45]; % 50ms before to 450ms after movement onset
edges_FR = window_FR(1):binSize_FR:window_FR(2); % bin edges for 10ms bins
FR_all_trials = [];

% determine movement type count and task IDs
moveType_num = numel(unique(All_SpikesPerMove_Tbl.MoveType));
moveType_ids = unique(All_SpikesPerMove_Tbl.MoveType);

% loop through movement types
for movT_i = 1:moveType_num
    move_subtbl = All_SpikesPerMove_Tbl(matches(All_SpikesPerMove_Tbl.MoveType, moveType_ids{movT_i}),:);

    % Subset by depth
    depth_num = cellfun(@(x) x(1), All_SpikesPerMove_Tbl.move_trial_ID, 'UniformOutput', false);
    depth_ids = unique(depth_num);

    % Clear after processing a depth for each movement type
    FR_all_trials = []; 

    % loop through depths
    depth_num3 = cellfun(@(x) x(1), move_subtbl.move_trial_ID, 'UniformOutput', false);
    for depth_i = 1:length(depth_ids)
        temp_trialID = [move_subtbl.move_trial_ID{depth_i},'_',move_subtbl.MoveType{depth_i},'_',num2str(move_subtbl.MoveN(depth_i))];

        temp_depth_idx = matches(depth_num3, depth_ids{depth_i}); % logical idx
        move_depth_tbl = move_subtbl(temp_depth_idx,:);

        temp_depth_row = move_subtbl(depth_i,:);
        temp_depth_spks = temp_depth_row.C1{1};
        % temp_depth_spks = vertcat(move_depth_tbl.C1{:});
        temp_depth_seconds = transpose(((temp_depth_spks - temp_depth_row.TTL_spk_idx_Start)/AO_spike_fs) - 0.05); % incorporate this in All_SpikesPerMove_Tbl
        % temp_depth_seconds = (temp_depth_spks - move_depth_tbl.TTL_spk_idx_Start(1)) / AO_spike_fs - 0.05;

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
    mean_FR = mean(FR_all_trials, 1);
    std_FR = std(FR_all_trials, 0, 1);

    % Debug output
    fprintf('Movement type %s: Mean FR = %f, STD FR = %f\n', moveType_ids{movT_i}, mean(mean_FR), mean(std_FR));

end


%% Compare mean FR for each trial per Motor Task type between STN depths

%%% automate programatically later %%%

% Index STN depth-specific data from All_SpikesPerMove_Tbl

% dorsal, All_SpikesPerMove_Tbl.move_trial_ID contains 't'
% central, All_SpikesPerMove_Tbl.move_trial_ID contains 'c'
% ventral, All_SpikesPerMove_Tbl.move_trial_ID contains 'b'

% switch All_SpikesPerMove_Tbl.move_trial_ID ...
%     case 'dorsal'
%         depthID = 't';
%     case 'central'
%         depthID = 'c';
%     case 'ventral'
%         depthID = 'b';
% end


%% Manually index STN depth-specific data from All_SpikesPerMove_Tbl by CaseDate

%%% automate programatically later %%%

% % 03/09/23 case
% % dorsalSTN_tbl, N/A
% centralSTN_tbl = All_SpikesPerMove_Tbl(55:end,1:9);
% ventralSTN_tbl = All_SpikesPerMove_Tbl(1:54,1:9);

% % 03/23/23 case   *
% dorsalSTN_tbl = All_SpikesPerMove_Tbl(111:end,1:9); % move_trial_ID contains 't'
% centralSTN_tbl = All_SpikesPerMove_Tbl(55:110,1:9); % move_trial_ID contains 'c'
% ventralSTN_tbl = All_SpikesPerMove_Tbl(1:54,1:9); % move_trial_ID contains 'b'

% % 04/05/23 case   *
% dorsal_REST_tbl = [All_SpikesPerMove_Tbl(23,1:9); All_SpikesPerMove_Tbl(46,1:9)];
% dorsalSTN_tbl = [All_SpikesPerMove_Tbl(17:22,1:9); All_SpikesPerMove_Tbl(40:45,1:9)]; % move_trial_ID contains 't'
% central_REST_tbl = [All_SpikesPerMove_Tbl(10,1:9); All_SpikesPerMove_Tbl(33,1:9)];
% centralSTN_tbl = [All_SpikesPerMove_Tbl(11:16,1:9); All_SpikesPerMove_Tbl(34:39,1:9)]; % move_trial_ID contains 'c'
% ventral_REST_tbl = [All_SpikesPerMove_Tbl(1,1:9); All_SpikesPerMove_Tbl(24,1:9)];
% ventralSTN_tbl = [All_SpikesPerMove_Tbl(2:9,1:9); All_SpikesPerMove_Tbl(25:32,1:9)]; % move_trial_ID contains 'b'

% % 04/13/23, LSTN case     *
% dorsal_REST_tbl = [All_SpikesPerMove_Tbl(41,1:9); All_SpikesPerMove_Tbl(48,1:9); All_SpikesPerMove_Tbl(55,1:9)];
% dorsalSTN_tbl = [All_SpikesPerMove_Tbl(42:47,1:9); All_SpikesPerMove_Tbl(49:54,1:9); All_SpikesPerMove_Tbl(56:end,1:9)]; % move_trial_ID contains 't'
% central_REST_tbl = [All_SpikesPerMove_Tbl(21,1:9); All_SpikesPerMove_Tbl(28,1:9)];
% centralSTN_tbl = [All_SpikesPerMove_Tbl(22:27,1:9); All_SpikesPerMove_Tbl(29:40,1:9)]; % move_trial_ID contains 'c'
% ventral_REST_tbl = [All_SpikesPerMove_Tbl(7,1:9); All_SpikesPerMove_Tbl(14,1:9)];
% ventralSTN_tbl = [All_SpikesPerMove_Tbl(1:6,1:9); All_SpikesPerMove_Tbl(8:13,1:9); All_SpikesPerMove_Tbl(15:20,1:9)]; % move_trial_ID contains 'b'

% % 05/11/23 case
% dorsalSTN_tbl = All_SpikesPerMove_Tbl(19:30,1:9); % move_trial_ID% contains 't' ~
% centralSTN_tbl = All_SpikesPerMove_Tbl(31:84,1:9); % move_trial_ID% contains 'c' ~
% ventralSTN_tbl, N/A

% % 05/18/23 (b), LSTN case     *
% dorsal_REST_tbl = [All_SpikesPerMove_Tbl(47,1:9); All_SpikesPerMove_Tbl(64:65,1:9)]; % move_trial_ID contains 't'
% dorsalSTN_tbl = [All_SpikesPerMove_Tbl(48:63,1:9); All_SpikesPerMove_Tbl(66:72,1:9)]; % move_trial_ID contains 't'
% central_REST_tbl = [All_SpikesPerMove_Tbl(39,1:9)]; % move_trial_ID contains 'c'
% centralSTN_tbl = [All_SpikesPerMove_Tbl(25:38,1:9); All_SpikesPerMove_Tbl(40:46,1:9)]; % move_trial_ID contains 'c'
% ventral_REST_tbl = [All_SpikesPerMove_Tbl(1,1:9); All_SpikesPerMove_Tbl(8,1:9); All_SpikesPerMove_Tbl(17,1:9)]; % move_trial_ID contains 'b'
% ventralSTN_tbl = [All_SpikesPerMove_Tbl(2:7,1:9); All_SpikesPerMove_Tbl(9:16,1:9); All_SpikesPerMove_Tbl(18:24,1:9)]; % move_trial_ID contains 'b'

% % 06/08/23, LSTN case
% dorsalSTN_tbl = All_SpikesPerMove_Tbl(71:end,1:9); % move_trial_ID contains 't'(singleTrial)
% centralSTN_tbl = All_SpikesPerMove_Tbl(54:70,1:9); % move_trial_ID contains 'c'(singleTrial)
% ventralSTN_tbl = All_SpikesPerMove_Tbl(1:18,1:9); % move_trial_ID contains 'b' (singleTrial)

% % 07_13_2023, LSTN case
% dorsalSTN_tbl = [All_SpikesPerMove_Tbl(39:51,1:9); All_SpikesPerMove_Tbl(53:59,1:9)]; % move_trial_ID contains 't'
% centralSTN_tbl = [All_SpikesPerMove_Tbl(17:30,1:9); All_SpikesPerMove_Tbl(32:37,1:9)]; % move_trial_ID contains 'c'
% ventralSTN_tbl = [All_SpikesPerMove_Tbl(1:7,1:9); All_SpikesPerMove_Tbl(9:14,1:9)]; % move_trial_ID contains 'b'

% % 07_13_2023, RSTN case     *
dorsal_REST_tbl = [All_SpikesPerMove_Tbl(17,1:9)]; % move_trial_ID contains 't'
dorsalSTN_tbl = [All_SpikesPerMove_Tbl(18:24,1:9)]; % move_trial_ID contains 't'
central_REST_tbl = [All_SpikesPerMove_Tbl(9,1:9)]; % move_trial_ID contains 'c'
centralSTN_tbl = [All_SpikesPerMove_Tbl(10:16,1:9)]; % move_trial_ID contains 'c'
ventral_REST_tbl = [All_SpikesPerMove_Tbl(1,1:9)]; % move_trial_ID contains 'b'
ventralSTN_tbl = [All_SpikesPerMove_Tbl(2:8,1:9)]; % move_trial_ID contains 'b'


%% Define STN depth of interest

%%% automate later %%%

STN_depth = 'dorsal'; % 't'
% STN_depth = 'central'; % 'c'
% STN_depth = 'ventral'; % 'b'

%% Compute FR mean and stdev for each STN depth

% loop through movement types
for movT_i = 1:moveType_num

    switch STN_depth
        case 'dorsal'
            move_subtbl = dorsalSTN_tbl; % 't'
        case 'central'
            move_subtbl = centralSTN_tbl; % 'c'
        case 'ventral'
            move_subtbl = ventralSTN_tbl; % 'b'
    end

    FR_all_trials = []; % Clear after processing a depth for each movement type

    % Subset by depth
    depth_num = cellfun(@(x) x(1), move_subtbl.move_trial_ID, 'UniformOutput', false);
    depth_ids = unique(depth_num);

    % loop through depths
    depth_num3 = cellfun(@(x) x(1), move_subtbl.move_trial_ID, 'UniformOutput', false);
    for depth_i = 1:length(depth_ids)
        temp_trialID = [move_subtbl.move_trial_ID{depth_i},'_',move_subtbl.MoveType{depth_i},'_',num2str(move_subtbl.MoveN(depth_i))];

        temp_depth_idx = matches(depth_num3, depth_ids{depth_i}); % logical idx
        move_depth_tbl = move_subtbl(temp_depth_idx,:);

        temp_depth_row = move_subtbl(depth_i,:);
        temp_depth_spks = temp_depth_row.C1{1};
        temp_depth_seconds = transpose(((temp_depth_spks - temp_depth_row.TTL_spk_idx_Start)/AO_spike_fs) - 0.05); % incorporate this in All_SpikesPerMove_Tbl

        % Extract spike times for current trial
        spikeTimes = temp_depth_seconds;

        % Create histogram of spike times
        spikeCounts_PSTH = histcounts(spikeTimes, edges_PSTH);
        spikeCounts_FR = histcounts(spikeTimes, edges_FR);

        % Store PSTH for this trials
        PSTH_all_trials = [PSTH_all_trials; spikeCounts_PSTH];

        % Convert spike counts to FR (spikes per second)
        FR_moverep = spikeCounts_FR/binSize_FR;

        % Store FR for move trial
        FR_all_trials = [FR_all_trials; FR_moverep];

    end

    % Compute mean and standard deviation firing rate across all trials
    d_mean_FR = mean(FR_all_trials, 1);
    d_std_FR = std(FR_all_trials, 0, 1);

    % Print output
    fprintf('Movement type %s: Mean FR (dorsal) = %f\n', moveType_ids{movT_i}, mean(d_mean_FR));

    % % compute depth-specific mean FR per move rep
    dorsal_FR_perMoveRep = transpose(d_mean_FR);
    mean_dorsalFR = mean(dorsal_FR_perMoveRep);
    [s_dFR, m_dFR] = std(dorsal_FR_perMoveRep)

    % loop through depths
    for depth_i = 1:length(depth_ids)
        temp_trialID = [move_subtbl.move_trial_ID{depth_i},'_',move_subtbl.MoveType{depth_i},'_',num2str(move_subtbl.MoveN(depth_i))];

        temp_depth_idx = matches(depth_num, depth_ids{depth_i}); % logical idx
        move_depth_tbl = move_subtbl(temp_depth_idx,:);

        temp_depth_row = move_subtbl(depth_i,:);
        temp_depth_spks = temp_depth_row.C1{1};
        % temp_depth_spks = vertcat(move_depth_tbl.C1{:});
        temp_depth_seconds = transpose(((temp_depth_spks - temp_depth_row.TTL_spk_idx_Start)/AO_spike_fs) - 0.05); % incorporate this in All_SpikesPerMove_Tbl
        % temp_depth_seconds = (temp_depth_spks - move_depth_tbl.TTL_spk_idx_Start(1)) / AO_spike_fs - 0.05;

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
    mean_FR = mean(FR_all_trials, 1);
    std_FR = std(FR_all_trials, 0, "all");

    % Debug output
    fprintf('Movement type %s: Mean FR = %f, STD FR = %f\n', moveType_ids{movT_i}, mean(mean_FR), mean(std_FR));

end

switch STN_depth
    case 'dorsal'
        dorsal_FR_perMoveRep = transpose(mean_FR);
        [std_dFR, mean_dFR] = std(dorsal_FR_perMoveRep)
    case 'central'
        central_FR_perMoveRep = transpose(mean_FR);
        [std_cFR, mean_cFR] = std(central_FR_perMoveRep)
    case 'ventral'
        ventral_FR_perMoveRep = transpose(mean_FR);
        [std_vFR, mean_vFR] = std(ventral_FR_perMoveRep)
end



