function kinTbl = run_MovementFeatureAnalysis_IO(CaseDate, CaseDate_hem, MoveDataDir, MoveDir_CaseID, fps, px_to_mm, zscoreKin)

% run_MovementFeatureAnalysis_IO
% Computes per-trial kinematic features for intraoperative motor recordings.
%
% Inputs:
% - CaseDate: e.g. '03_23_2023'
% - CaseDate_hem: 'LSTN' or 'RSTN'
% - MoveDataDir: full path to top-level 'Processed DLC' directory
% - MoveDir_CaseID: subfolder with movement CSVs, e.g., 'IO_03_23_2023_LSTN'
% - fps: video framerate (e.g. 100)
% - px_to_mm: scale conversion (e.g. 2.109 mm/px)
% - zscoreKin: true/false toggle to z-score output features

fprintf('[INFO] Running movement feature analysis for: %s\n', MoveDir_CaseID);

%% Setup
curPC = getenv('COMPUTERNAME');
switch curPC
    case 'DESKTOP-I5CPDO7'
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    otherwise
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
end
output_kinAnalysisDir = fullfile(IO_DataDir, 'Kinematic Analyses', [CaseDate '_' CaseDate_hem]);
if ~exist(output_kinAnalysisDir, 'dir')
    mkdir(output_kinAnalysisDir);
end

% Define subdirs
caseDir = fullfile(MoveDataDir, MoveDir_CaseID);
vidDir  = fullfile(caseDir, 'video folder');
csvDir  = fullfile(caseDir, 'csv folder');

moveFiles = dir(fullfile(vidDir, '*Move*.csv'));
fprintf('[INFO] Found %d Move Index files\n', numel(moveFiles));

kinTbl = table();  % Initialize output

%% Loop over movement trials
for i = 1:numel(moveFiles)
    try
        moveFile = moveFiles(i).name;
        moveIndex = readtable(fullfile(vidDir, moveFile));
        moveIndex = moveIndex(moveIndex.BeginF > 0 & moveIndex.EndF > 0, :);

        % Match DLC CSV file
        baseName = erase(moveFile, '.csv');
        dlcPath = fullfile(csvDir, [baseName, '.csv']);
        if ~isfile(dlcPath)
            warning('[SKIP] DLC CSV not found for %s', baseName); continue;
        end

        rawDLC = readtable(dlcPath);
        if any(contains(rawDLC.Properties.VariableNames, 'likelihood'))
            rawDLC = cleanLowConfidence(rawDLC, 0.6);
        end

        % Use index fingertip marker
        X = rawDLC.index_tip_x;
        Y = rawDLC.index_tip_y;
        dist = sqrt(diff(X).^2 + diff(Y).^2) * px_to_mm;
        time = (1:length(dist)) / fps;
        smoothed = smoothdata(dist, 'gaussian', 15);

        for j = 1:height(moveIndex)
            bF = moveIndex.BeginF(j);
            eF = min(moveIndex.EndF(j), length(smoothed));

            seg = smoothed(bF:eF);
            seg_t = time(bF:eF);

            % Detect peaks
            [pks, locs, widths] = findpeaks(seg, ...
                'MinPeakProminence', std(seg)*0.5, ...
                'MinPeakDistance', round(0.15 * fps));

            if isempty(pks), continue; end

            amp = pks;
            ampMean = mean(amp);
            repDur = diff(locs) / fps;
            interDur = repDur(1:end-1);
            vel = amp ./ (widths(1:length(amp)) / fps); % rough

            % Store results
            newRow = table( ...
                moveIndex.MoveN(j), ...
                string(moveIndex.MoveType{j}), ...
                bF, eF, ...
                ampMean, ...
                mean(repDur, 'omitnan'), ...
                mean(interDur, 'omitnan'), ...
                mean(vel, 'omitnan'), ...
                'VariableNames', {'MoveN','MoveType','BeginF','EndF', ...
                                  'MeanAmp','MeanRepDur','MeanInterRepDur','MeanVel'});
            kinTbl = [kinTbl; newRow];
        end
    catch ME
        warning('[ERROR] %s: %s', moveFile, ME.message); continue;
    end
end

%% Append metadata: move_trial_ID, Hemisphere, Depth
kinTbl.move_trial_ID = strcat(lower(MoveDir_CaseID(15)), string(kinTbl.MoveN), '_', lower(MoveDir_CaseID(15)));
kinTbl.Hemisphere = repmat({CaseDate_hem}, height(kinTbl), 1);
if contains(MoveDir_CaseID, 't')
    kinTbl.Depth = repmat({'dorsal'}, height(kinTbl), 1);
elseif contains(MoveDir_CaseID, 'c')
    kinTbl.Depth = repmat({'central'}, height(kinTbl), 1);
elseif contains(MoveDir_CaseID, 'b')
    kinTbl.Depth = repmat({'ventral'}, height(kinTbl), 1);
else
    kinTbl.Depth = repmat({'unknown'}, height(kinTbl), 1);
end

%% Optional z-scoring
if zscoreKin
    kinCols = {'MeanAmp','MeanRepDur','MeanInterRepDur','MeanVel'};
    kinTbl{:, kinCols} = zscore(kinTbl{:, kinCols});
end

%% Save CSV
outfile = fullfile(output_kinAnalysisDir, 'Kinematic_Features.csv');
writetable(kinTbl, outfile);
fprintf('[SAVED] Kinematic features written to:\n%s\n', outfile);
fprintf('[DONE] %d entries processed.\n', height(kinTbl));

end

%% === Subfunction: confidence filter ===
function out = cleanLowConfidence(tbl, CT)
out = tbl;
cols = tbl.Properties.VariableNames;
for i = 1:3:length(cols)
    if contains(cols{i+2}, 'likelihood')
        bad = tbl{:, cols{i+2}} < CT;
        out{bad, cols{i}}   = NaN;
        out{bad, cols{i+1}} = NaN;
    end
end
end
