function [] = comps_IO_plotting(plot_ID, case_ID, ephys_offset)

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

%% Hardcode Case-specific Data directories

% isolate a specific CaseDate / studyID (StudyNum in Subject_AO csv)
CaseDate = '03_09_2023'; % studyID = 1
% CaseDate = '03_23_2023'; % studyID = 2
case_ID = CaseDate;

% define case-specific data directory
Case_DataDir = [IO_DataDir, filesep, CaseDate];
ephysTbl_Dir = [Case_DataDir, filesep, 'DLC_MER'];


%% define datastream sampling rates (Alpha Omega and Video sampling fs)

TTL_fs = 44000; % Hz
AO_spike_fs = 44000; % Hz
AO_LFP_fs = 1375; % Hz
DLC_fs = 100; % fps


%% switch case setup

addpath 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\AO_Processing'

% raster                % offset = 50ms before trial start
switch plot_ID
    case 'raster_fig'   % spike response during movement reps
        cd(ephysTbl_Dir)
        Tbl_list = dir('*.mat');
        Tbl_names = {Tbl_list.name};
        if ephys_offset
            search_index = contains(Tbl_names, 'Spikes') & contains(Tbl_names, 'offset');
        else
            search_index = contains(Tbl_names, 'Spikes') & ~contains(Tbl_names, 'offset');
        end
        spk_case = Tbl_names{search_index};
        load(spk_case, 'All_SpikesPerMove_Tbl');
        % plot all 3 depths
        moveType_num = numel(unique(All_SpikesPerMove_Tbl.MoveType));
        moveType_ids = unique(All_SpikesPerMove_Tbl.MoveType);

        tiledlayout(1,moveType_num)
        for fig_i = 1:moveType_num
            move_subtbl = All_SpikesPerMove_Tbl(matches(All_SpikesPerMove_Tbl.MoveType, moveType_ids{fig_i}),:);

            depth_num = cellfun(@(x) x(1), All_SpikesPerMove_Tbl.move_trial_ID, 'UniformOutput', false);
            depth_ids = unique(depth_num);

            % uisetcolor, coloors.com
            color_vec = [0    0.4471    0.7412; ...
                0.7176    0.2745    1.0000; ...
                0.4667    0.6745    0.1882];

            % initialize legend info and color mapping
            legend_info = {};
            depth_color_idx = [];
            legend_handles = []; % gobjects(0) % empty array of graphic objects

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
                if ~ismember(depth_idx, depth_color_idx)
                    legend_info{end+1} = sprintf('depth %s', depth_ids{depth_idx});
                    depth_color_idx(end+1) = depth_idx;
                    depth_color = color_vec(depth_idx,:);
                    legend_handles = [legend_handles; m_line];
                end

            end

            title(moveType_ids{fig_i})
            xlabel('time (s)')
            ylabel('movement repetition')

        end

        % Add legend on top-right, outside the plot
        legend(legend_handles, legend_info, 'Location', 'northeastoutside');

        % case 'FR_fig'     % firing rate (spike/s)

end








end