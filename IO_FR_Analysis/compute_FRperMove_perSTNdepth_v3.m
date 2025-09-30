function [FR_perTrialRep_All, FR_perMoveType_perDepth_Summary] = compute_FRperMove_perSTNdepth_v3(CaseDate, All_SpikesPerMove_Tbl, AO_spike_fs, Case_FRKin_Dir)

% Compute mean FR per Movement context per STN depth, Export full Data Table

%% Auto-split STN depths

dorsalSTN_tbl  = All_SpikesPerMove_Tbl(contains(All_SpikesPerMove_Tbl.move_trial_ID,'t'),:);
centralSTN_tbl = All_SpikesPerMove_Tbl(contains(All_SpikesPerMove_Tbl.move_trial_ID,'c'),:);
ventralSTN_tbl = All_SpikesPerMove_Tbl(contains(All_SpikesPerMove_Tbl.move_trial_ID,'b'),:);

depthTables = struct('dorsal', dorsalSTN_tbl, 'central', centralSTN_tbl, 'ventral', ventralSTN_tbl);

moveType_ids = unique(All_SpikesPerMove_Tbl.MoveType);

%% Compute FR for each depth & move type

% window_FR  = [-0.05 0.45]; % 50 ms pre, 450 ms post

% Initialize storage containers
FR_Results = struct();
% summary_Data_tbl = {};
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

        %% run MaxSpkDuration_Raster_PSTH
        % [Max_SpikeDuration_samples, spikesMatrix] = MaxSpkDuration_Raster_PSTH(CaseDate, move_subtbl, AO_spike_fs, Case_FRKin_Dir);

        % Calculate the max spike segment duration (# of samples in longest segment)
        max_SpkDuration_All = zeros(height(move_subtbl), 1);

        for row_i = 1:height(move_subtbl)
            tempSpk_vec = move_subtbl.C1{row_i};
            if isempty(tempSpk_vec) || numel(tempSpk_vec) < 3
                max_SpkDuration_All(row_i) = NaN;
            else
                max_SpkDuration_All(row_i) = max(tempSpk_vec - move_subtbl.TTL_spk_idx_Start(row_i));
            end
        end

        Max_SpikeDuration_samples = max(max_SpkDuration_All); % samples in longest segment
        Max_SpkDur_seconds = Max_SpikeDuration_samples/AO_spike_fs; % seconds
        Max_SpkDus_ms = Max_SpkDur_seconds * 1000; % milliseconds


        spikesMatrix = zeros(height(move_subtbl), Max_SpikeDuration_samples, 'logical');
        for row_i = 1:height(move_subtbl)
            tempSpk_vec = move_subtbl.C1{row_i};
            if isempty(tempSpk_vec) || numel(tempSpk_vec) < 3
                continue;
            else
                spike_indices = tempSpk_vec - move_subtbl.TTL_spk_idx_Start(row_i);
                spike_indices = spike_indices(spike_indices >= 1 & spike_indices <= Max_SpikeDuration_samples);
                spikesMatrix(row_i, spike_indices) = true;
            end

        end


        %%

        % --- Compute FR per trial inside segment window (TTL start→end) ---
        FR_all_trials = nan(height(move_subtbl),1);
        Dur_all_trials_sec = nan(height(move_subtbl),1);
        SpikeCount_all_trials = nan(height(move_subtbl),1);

        for row_i = 1:height(move_subtbl)
            tempSpk_vec = move_subtbl.C1{row_i};
            if isempty(tempSpk_vec) || numel(tempSpk_vec) < 1
                continue;
            end
            t_start = move_subtbl.TTL_spk_idx_Start(row_i);
            t_end   = move_subtbl.TTL_spk_idx_End(row_i);

            % spikes inside [start, end]
            in_window = tempSpk_vec(tempSpk_vec >= t_start & tempSpk_vec <= t_end);
            dur_sec  = (t_end - t_start) / AO_spike_fs; % seconds

            Dur_all_trials_sec(row_i)     = dur_sec;
            SpikeCount_all_trials(row_i) = numel(in_window);
            FR_all_trials(row_i)        = numel(in_window) / dur_sec; % Hz
        end

        % Attach per-trial FR to the move_subtbl (for exporting later)
        move_subtbl.FR_Hz          = FR_all_trials;
        move_subtbl.FR_Duration_s  = Dur_all_trials_sec;
        move_subtbl.FR_SpikeCount  = SpikeCount_all_trials;

        % Keep only relevant columns for export
        keepCols = {'CasDate','move_trial_ID','MoveType','Depth','TTL_spk_idx_Start','TTL_spk_idx_End','FR_Hz','FR_Duration_s','FR_SpikeCount'};
        keepCols = keepCols(ismember(keepCols, move_subtbl.Properties.VariableNames));
        full_FR_Data{end+1,1} = move_subtbl(:, keepCols);

        % Summary stats for this depth × MoveType
        mu_FR  = mean(FR_all_trials, 'omitnan');
        sd_FR  = std(FR_all_trials,  'omitnan');
        n_FR   = sum(~isnan(FR_all_trials));

        fprintf('Move %s: Mean FR = %.2f Hz, STD FR = %.2f Hz\n', ...
            moveType_ids{movT_i}, mean(mu_FR), mean(sd_FR));

        % Safe movement name for struct field
        safeMoveName = matlab.lang.makeValidName(moveType_ids{movT_i});  % converts 'ARM EF' -> 'ARM_EF'

        % Store in results struct
        FR_Results.(depthName{1}).(safeMoveName).mean_FR = mu_FR;
        FR_Results.(depthName{1}).(safeMoveName).std_FR  = sd_FR;
        FR_Results.(depthName{1}).(safeMoveName).n_FR = n_FR;
        FR_Results.(depthName{1}).(safeMoveName).FR_all_trials = FR_all_trials;
        FR_Results.(depthName{1}).(safeMoveName).Duration_all_trials_sec = Dur_all_trials_sec;
        FR_Results.(depthName{1}).(safeMoveName).SpikeCount_all_trials = SpikeCount_all_trials;

        % summary_Data_tbl{end+1,1} = table(...
        %     string(depthName{1}), string(moveType_ids{movT_i}), n_FR, mu_FR, sd_FR, ...
        %     'VariableNames', {'Depth','MoveType','N','Mean_FR_Hz','SD_FR_Hz'});

        % Append to summaryData (for CSV 1)
        summary_Data(end+1,:) = {CaseDate, depthName{1}, moveType_ids{movT_i}, ...
                                n_FR, mu_FR, sd_FR};


    end
end

% Append CaseDate to struct
FR_Results.CaseDate = CaseDate;

%% Concatenate for export

if ~isempty(full_FR_Data)
    FR_perTrialRep_All = vertcat(full_FR_Data{:});
else
    FR_perTrialRep_All = table();
end


%% Convert summaryData cell array to table for export

FR_perMoveType_perDepth_Summary = cell2table(summary_Data, ...
    'VariableNames', {'CaseDate','STN_Depth','MoveType','N','MeanFR_Hz','StdFR_Hz'});

FR_perMoveTperDepth_summary = ['FR_perMoveType_perDepth_Summary_' CaseDate '.csv'];

%% Save outputs

cd(Case_FRKin_Dir)
save(['FR_Results_' CaseDate '.mat'], 'FR_Results');

% Write CSVs
writetable(FR_perTrialRep_All, fullfile(Case_FRKin_Dir, sprintf('%s_FR_perTrialRep_All.csv', CaseDate)));
writetable(FR_perMoveType_perDepth_Summary, FR_perMoveTperDepth_summary);

fprintf('\nSaved per-trial FR and summary CSVs to: %s\n', Case_FRKin_Dir);

%%