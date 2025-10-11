function aggregate_FRKinematic_Correlations(KinematicsDir, Ephys_Kin_Dir, FR_Kin_Dir)

% Aggregates all Correlation_Summary.csv files from subfolders and visualizes results

%% Define output dir

output_agg_FRKinCorr = [Ephys_Kin_Dir, filesep, 'aggregated_FRKinCorr'];
if ~exist(output_agg_FRKinCorr, 'dir')
    mkdir(output_agg_FRKinCorr);
end

% Initialize output storage containers
% agg_Master_CorrTbl = table(); % for all 'Correlation Summary.csv' files
agg_Merged_FRKinTbl = table(); % for all 'Merged_FR_Kin_Summary.csv' files


%% Loop through each CaseDate in FR_Kin_Dir
%
% after running run_FR_KinematicCorr function 
%
% cd(FR_Kin_Dir)
%
% % Navigate to the 'FR_Kinematics' subfolder in each CaseDate folder
% caseFolders = dir(FR_Kin_Dir);
% caseFolders = caseFolders([caseFolders.isdir] & ~startsWith({caseFolders.name}, '.'));
%
% % Loop through each CaseDate in FR_Kin_Dir
% for case_i = 1:numel(caseFolders)
%     caseName = caseFolders(case_i).name;
%     casePath = [FR_Kin_Dir, filesep, caseName, filesep, 'FR_Kinematics'];
%
%     % Skip if not a directory
%     if ~isfolder(casePath), continue; end
%
%     % Extract the 'Correlation Summary.csv' and 'Merged_FR_Kin_Summary.csv' files
%     CorrSummaryPath = fullfile(casePath, 'Correlation_Summary.csv');
%
%     % Aggregate all 'Correlation Summary.csv' files
%     if exist(CorrSummaryPath, 'file')
%         try
%             corrTbl = readtable(CorrSummaryPath);
%             corrTbl.CaseDate = repmat({caseName}, height(corrTbl), 1);
%
%             agg_Master_CorrTbl = [agg_Master_CorrTbl; corrTbl];
%             fprintf('[ADDED] Correlation_Summary.csv from: %s\n', caseName);
%         catch
%             warning('[SKIP] Failed to read: %s', CorrSummaryPath);
%         end
%     else
%         fprintf('[MISSING] No Correlation_Summary.csv in: %s\n', caseName);
%     end
% end

%% Loop through each MoveDir_CaseID in KinematicsDir

cd(KinematicsDir)

% Define MoveDir_CaseIDs in KinematicsDir
% MoveCaseFolders = dir(KinematicsDir);
% MoveCaseFolders = MoveCaseFolders([MoveCaseFolders.isdir] & ~startsWith({MoveCaseFolders.name}, '.'));

MoveCaseFolders = {'IO_03_23_2023_LSTN', 'IO_04_05_2023_RSTN', 'IO_05_18_2023_b_LSTN', 'IO_05_31_2023_LSTN', 'IO_06_08_2023_LSTN', ...
    'IO_07_06_2023_LSTN', 'IO_07_13_2023_LSTN', 'IO_08_23_2023_RSTN'};
% 'IO_03_23_2023_LSTN'; % NER 2025
% 'IO_04_05_2023_RSTN'; % NER 2025
% 'IO_05_18_2023_b_LSTN'; % NER 2025
% 'IO_05_31_2023_LSTN'; % INS 2026
% 'IO_06_08_2023_LSTN'; % NER 2025
% 'IO_07_06_2023_LSTN'; % INS 2026
% 'IO_07_13_2023_LSTN'; % % INS 2026
% 'IO_08_23_2023_RSTN'; % NANS 2026

% Loop through each MoveDir_CaseID in KinematicsDir
for movecase_i = 1:numel(MoveCaseFolders)
    MoveCaseName = MoveCaseFolders{movecase_i};

    % Extract the 'Merged_FR_Kin_Summary.csv' file
    masterFRKinPath = fullfile(KinematicsDir, MoveCaseFolders{movecase_i}, 'Merged_FR_Kin_Summary.csv');

    % Aggregate all 'Master_FR_Kin.csv' files
    % fprintf('[LOOKING] Merged_FR_Kin_Summary.csv in: %s\n', masterFRKinPath);
    % if exist(masterFRKinPath, 'file')
    % try
    master_FRKin_Tbl = readtable(masterFRKinPath);
    if matches('Hemisphere', master_FRKin_Tbl.Properties.VariableNames)
        master_FRKin_Tbl = removevars(master_FRKin_Tbl, 'Hemisphere');
    end
    master_FRKin_Tbl.CaseDate = repmat(MoveCaseFolders(movecase_i), height(master_FRKin_Tbl), 1);
    agg_Merged_FRKinTbl = [agg_Merged_FRKinTbl; master_FRKin_Tbl];
    fprintf('[ADDED] Master_FR_Kin.csv from: %s\n', MoveCaseName);
    % catch
    %     warning('[SKIP] Failed to read: %s', masterFRKinPath);
    % end
    % else
    %     fprintf('[MISSING] No Master_FR_Kin.csv in: %s\n', MoveCaseName);
    % end
end


%% Correlation analysis per MoveType and Depth

movTypes = unique(agg_Merged_FRKinTbl.MoveType);
depths   = unique(agg_Merged_FRKinTbl.Depth);

% Kinematic features to correlate (keep simple for now)
kinFeatures = {'MeanAmp', 'MeanVel'};

% Table to store correlation results
corrResults = table([], [], [], [], [], [], [], [], ...
    'VariableNames', {'MovementType','Depth','KinematicFeature','N','R','p','FisherZ','Method'});

% Table to store correlation results
corrResults2 = table([], [], [], [], [], [], [], ...
    'VariableNames', {'Depth','KinematicFeature','N','R','p','FisherZ','Method'});


% % Convert MoveType and Depth to usable character format
% if iscategorical(agg_Merged_FRKinTbl.MoveType)
%     agg_Merged_FRKinTbl.MoveType = cellstr(agg_Merged_FRKinTbl.MoveType);
% end
% if iscategorical(agg_Merged_FRKinTbl.Depth)
%     agg_Merged_FRKinTbl.Depth = cellstr(agg_Merged_FRKinTbl.Depth);
% elseif isnumeric(agg_Merged_FRKinTbl.Depth)
%     agg_Merged_FRKinTbl.Depth = cellstr(string(agg_Merged_FRKinTbl.Depth));
% end


%% Perform correlation analysis

movTypes = movTypes(~matches(movTypes, 'REST'));
for m = 1:numel(movTypes)
    for d = 1:numel(depths)
        fprintf('  %s × %s: ', movTypes{m}, depths{d});
        sel = matches(agg_Merged_FRKinTbl.MoveType, movTypes{m}) & matches(agg_Merged_FRKinTbl.Depth, depths{d});
        fprintf('%d trials\n', sum(sel));

        for k = 1:numel(kinFeatures)
            kinFeat = kinFeatures{k};

            % Select subset of trials
            sub = agg_Merged_FRKinTbl(sel, :);
            if height(sub) < 4 || sum(isnan(sub.(kinFeat))) > 4
                continue;
            end

            % Ensure any NaNs in a feature are removed from both variable columns
            if any(isnan(sub.(kinFeat)))
                sub = sub(~isnan(sub.(kinFeat)), :);
            end
            % Ensure any 0s in meanFR col are removed from both variable columns
            if any(sub.Mean_FR_Hz == 0)
                sub = sub(sub.Mean_FR_Hz ~= 0, :);
            end
            fprintf('%d trials\n', sum(sel));

            x_corr_FR = sub.Mean_FR_Hz;
            y_corr_kinFts = sub.(kinFeat);


            % Pearson correlation for all trials
            [R, p] = corr(x_corr_FR, y_corr_kinFts, 'Type', 'Pearson');
            method = 'Pearson';  % using Pearson correlation exclusively for now

            % Fisher z-transform for effect size
            eff = atanh(R);

            % Store result
            newRow = {movTypes{m}, depths{d}, kinFeat, height(sub), R, p, eff, method};
            corrResults = [corrResults;
                cell2table(newRow, 'VariableNames', ...
                {'MovementType','Depth','KinematicFeature','N','R','p','FisherZ','Method'})];

            % Save scatter plot
            fig = figure('Visible', 'off');
            scatter(x_corr_FR, y_corr_kinFts, 50, 'filled'); hold on;
            lsline;
            xlabel('Firing Rate (Hz)');
            ylabel(kinFeat);
            title(sprintf('%s | %s | %s: R = %.2f, p = %.3f (%s)', ...
                depths{d}, movTypes{m}, kinFeat, R, p, method));
            saveas(fig, fullfile(output_agg_FRKinCorr, ...
                sprintf('%s_%s_%s_corr.png', movTypes{m}, depths{d}, kinFeat)));
            close(fig);
        end
    end
end


% %% Perform correlation analysis
% 
% for d = 1:numel(depths)
%     sel = matches(agg_Merged_FRKinTbl.Depth, depths{d});
%     fprintf('%d trials\n', sum(sel));
% 
%     kinFeat = 'MeanVel';
% 
%     % Select subset of trials
%     sub = agg_Merged_FRKinTbl(sel, :);
%     if height(sub) < 4 || sum(isnan(sub.(kinFeat))) > 4
%         continue;
%     end
% 
%     % Ensure any NaNs in a feature are removed from both variable columns
%     if any(isnan(sub.(kinFeat)))
%         sub = sub(~isnan(sub.(kinFeat)), :);
%     end
%     % Ensure any 0s in meanFR col are removed from both variable columns
%     if any(sub.Mean_FR_Hz == 0)
%         sub = sub(sub.Mean_FR_Hz ~= 0, :);
%     end
%     fprintf('%d trials\n', sum(sel));
% 
%     x_corr_FR = sub.Mean_FR_Hz;
%     y_corr_kinFts = sub.(kinFeat);
% 
% 
%     % Pearson correlation for all trials
%     [R, p] = corr(x_corr_FR, y_corr_kinFts, 'Type', 'Pearson');
%     method = 'Pearson';  % using Pearson correlation exclusively for now
% 
%     % Fisher z-transform for effect size
%     eff = atanh(R);
% 
%     % Store result
%     newRow = {depths{d}, kinFeat, height(sub), R, p, eff, method};
%     corrResults2 = [corrResults2;
%         cell2table(newRow, 'VariableNames', ...
%         {'Depth','KinematicFeature','N','R','p','FisherZ','Method'})];
% end
% test = 1;

%%
% Run kai square analysis on the R-values for each MoveType
% Binary - negative or not
% Increase in Vel, Decrease in FR
% Mixed for Central and Ventral, All 1s (all neg) for Dorsal


% 6/6 motor trials in dorsal STN revealed a negative correlation (mean R =
% +- stdev R) between FR and movement velocity.
% Motor inhibition and breaking in sensorimotor region

%% Save aggregated

writetable(corrResults,  fullfile(output_agg_FRKinCorr, 'All_Correlation_Summary.csv'));
writetable(agg_Merged_FRKinTbl, fullfile(output_agg_FRKinCorr, 'All_Master_FR_Kin.csv'));

fprintf('[DONE] Aggregated CSV(s) saved to: %s\n', output_agg_FRKinCorr);
cd(output_agg_FRKinCorr)

end

