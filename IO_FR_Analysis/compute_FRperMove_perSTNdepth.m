function [] = compute_FRperMove_perSTNdepth(CaseDate, ephysTbl_Dir, All_SpikesPerMove_Tbl)

% Compute mean FR per Movement context per STN depth, Export full Data Table

addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis'

%% Sampling Rates

AO_spike_fs = 44000; % Hz
DLC_fs = 100; % fps

%% Auto-split STN depths

dorsalSTN_tbl  = All_SpikesPerMove_Tbl(contains(All_SpikesPerMove_Tbl.move_trial_ID,'t'),:);
centralSTN_tbl = All_SpikesPerMove_Tbl(contains(All_SpikesPerMove_Tbl.move_trial_ID,'c'),:);
ventralSTN_tbl = All_SpikesPerMove_Tbl(contains(All_SpikesPerMove_Tbl.move_trial_ID,'b'),:);

depthTables = struct('dorsal', dorsalSTN_tbl, 'central', centralSTN_tbl, 'ventral', ventralSTN_tbl);

%% Define FR parameters

binSize_FR = 0.01; % 10 ms bin size
window_FR  = [-0.05 0.45]; % 50 ms pre, 450 ms post
edges_FR   = window_FR(1):binSize_FR:window_FR(2);

moveType_ids = unique(All_SpikesPerMove_Tbl.MoveType);

%% Compute FR for each depth & move type

% Initialize storage containers
FR_Results = struct();
summary_Data = {};
full_FR_Data  = {};

for depthName = {'dorsal','central','ventral'}
    depth_tbl = depthTables.(depthName{1});
    if isempty(depth_tbl)
        fprintf('\n[INFO] No trials for %s STN\n', depthName{1});
        continue;
    end

    fprintf('\n=== %s STN ===\n', upper(depthName{1}));

    for movT_i = 1:numel(moveType_ids)
        move_subtbl = depth_tbl(strcmp(depth_tbl.MoveType, moveType_ids{movT_i}),:);
        if isempty(move_subtbl)
            continue;
        end
        % if numel(move_subtbl) < 3
        %     continue;
        % end

        FR_all_trials = [];
        spikeTimes_all = {}; % For raster

        for row_i = 1:height(move_subtbl)
            spkTimes = ((move_subtbl.C1{row_i} - move_subtbl.TTL_spk_idx_Start(row_i)) / AO_spike_fs) - 0.05;
            spikeTimes_all{end+1} = spkTimes; 
            spikeCounts = histcounts(spkTimes, edges_FR);
            FR_all_trials = [FR_all_trials; spikeCounts/binSize_FR];
        end

        % Compute FR mean and std across trials
        mean_FR = mean(FR_all_trials, 1);
        std_FR  = std(FR_all_trials, 0, 1);

        fprintf('Move %s: Mean FR = %.2f Hz, STD FR = %.2f Hz\n', ...
            moveType_ids{movT_i}, mean(mean_FR), mean(std_FR));

        % Safe movement name for struct field
        safeMoveName = matlab.lang.makeValidName(moveType_ids{movT_i});  % converts 'ARM EF' -> 'ARM_EF'

        % Store in results struct
        FR_Results.(depthName{1}).(safeMoveName).mean_FR = mean_FR;
        FR_Results.(depthName{1}).(safeMoveName).std_FR  = std_FR;
        FR_Results.(depthName{1}).(safeMoveName).FR_all_trials = FR_all_trials;
        FR_Results.(depthName{1}).(safeMoveName).BinEdges = edges_FR;
        FR_Results.(depthName{1}).(safeMoveName).SpikeTimes = spikeTimes_all;

        % Append to summaryData (for CSV 1)
        summary_Data(end+1,:) = {CaseDate, depthName{1}, moveType_ids{movT_i}, ...
                                mean(mean_FR), mean(std_FR)};

        % Append to fullFRData (for CSV 2)
        full_FR_Data(end+1,:) = {CaseDate, depthName{1}, moveType_ids{movT_i}, ...
                               strjoin(string(edges_FR),','), ...
                               strjoin(string(mean_FR),',')};

    end
end
   
% Append CaseDate to struct
FR_Results.CaseDate = CaseDate;


%% Convert to summary cell arrays to tables for export

FR_Summary_Table = cell2table(summary_Data, ...
    'VariableNames', {'CaseDate','STN_Depth','MoveType','MeanFR_Hz','StdFR_Hz'});

FR_Full_Table = cell2table(full_FR_Data, ...
    'VariableNames', {'CaseDate','STN_Depth','MoveType','BinEdges','FR_BinValues'});

%% Save outputs

cd(ephysTbl_Dir)
save(['FR_Results_' CaseDate '.mat'], 'FR_Results');

summaryFile = ['FR_Summary_Table_' CaseDate '.csv'];
fullFile    = ['FR_FullPSTH_Table_' CaseDate '.csv'];

writetable(FR_Summary_Table, summaryFile);
writetable(FR_Full_Table, fullFile);

fprintf('\n[INFO] Exported:\n - %s\n - %s\n - FR_Results_%s.mat\n', summaryFile, fullFile, CaseDate);

%%