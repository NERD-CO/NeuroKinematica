function [kinTbl, kinSummaryTbl] = run_MovementFeatureAnalysis_IO_v2(IO_DataDir, MoveDataDir, MoveDir_CaseID)

% run_MovementFeatureAnalysis_IO
% Extract movement features from DLC-labeled kinematic CSVs using Movement Index files
% Computes per-trial kinematic features for intraoperative motor recordings.

% Inputs:
% - CaseDate: e.g. '03_23_2023'
% - CaseDate_hem: 'LSTN' or 'RSTN'
% - MoveDataDir: full path to top-level 'Processed DLC' directory
% - MoveDir_CaseID: subfolder with movement CSVs, e.g., 'IO_03_23_2023_LSTN'
% - zscoreKin: true/false toggle to z-score output features

% Outputs:
%   kinTbl - trial-wise feature table with MeanAmp, StdAmp, MeanVel, StdVel
%   kinSummaryTbl - summary of MeanRepDur & MeanInterRepDur per MoveType and trial ID

fprintf('[INFO] Running movement feature analysis for: %s\n', MoveDir_CaseID);

%% Define MoveData subdirs

caseDir = fullfile(MoveDataDir, MoveDir_CaseID);
vidDir  = fullfile(caseDir, 'video folder'); % Contains DLC-labeled videos + Movement Index CSVs
csvDir  = fullfile(caseDir, 'csv folder'); % Contains Kinematic Timeseries data

moveFiles = dir(fullfile(vidDir, '*Move*.csv'));
fprintf('[INFO] Found %d Move Index files\n', numel(moveFiles));


%% Define output dir

output_kinAnalysisDir = fullfile(IO_DataDir, 'Kinematic Analyses', MoveDir_CaseID);
if ~exist(output_kinAnalysisDir, 'dir')
    mkdir(output_kinAnalysisDir);
end

% Initialize output storage containers
kinTbl = table();
kinSummaryTbl = table();   % For storing per trial+type repetition durations
% SkippedTbl = table();      % For logging short/failed peak trials
SkippedTbl = table([], [], [], ...
    'VariableNames', {'TrialID','FileName','Reason'});


%% Sampling rates and Calibration

fps = 100; % IO video frame rate

% calibration / conversion metric
px_to_mm = 2.109; % distance conversion factor

% update to calibration-image based method later
% % Load checkerboard image
% I = imread('calibration_frame.png');
% [imagePoints, boardSize] = detectCheckerboardPoints(I);
% squareSize_mm = 25.4;
% pixels_to_mm = squareSize_mm; % assuming 1 square in image = 25.4 mm
%
% % Or estimate from image width:
% avg_square_px = mean(diff(imagePoints(:,1)));
% pixels_to_mm = squareSize_mm / avg_square_px;


%% Main Loop over movement trials

for i = 1:numel(moveFiles)
    moveFile = moveFiles(i).name;

    % Extract base name of moveFiles (Movement Index CSVs in vidDir)
    [~, basePrefix, ~] = fileparts(moveFile);
    % disp(basePrefix);

    % Attempt to extract coreID using flexible patterns:
    % Matches:
    %   20230823_b1_d2p6_session001
    %   20230823_R_b1_d2p6_session001
    %   20230823_LH_b1_d2p6_session001
    core_pat = '(20\d{6}_(?:[LR]H|[LR])?_[bct]\d+_d\d+p\d+_session\d+)';
    tokC = regexp(basePrefix, core_pat, 'tokens', 'once');

    if ~isempty(tokC)
        coreID = tokC{1};
    else
        % Fall back to older forms just in case
        tokA = regexp(basePrefix, '(20\d{6}_[bct]\d+_d\d+p\d+_session\d+)', 'tokens', 'once');
        tokB = regexp(basePrefix, '(20\d{6}_[LR]H_[bct]\d+_d\d+p\d+_session\d+)', 'tokens', 'once');
        if ~isempty(tokA)
            coreID = tokA{1};
        elseif ~isempty(tokB)
            coreID = tokB{1};
        else
            SkippedTbl = [SkippedTbl; {string(basePrefix), moveFile, "CoreID_Extraction_Failed"}];
            continue;
        end
    end



    %  Find matching DLC CSV
    dlcMatch = dir(fullfile(csvDir, ['*' coreID '*DLC*.csv']));
    if isempty(dlcMatch)
        % Fallback: some exports may not include 'DLC' in the filename
        dlcMatch = dir(fullfile(csvDir, ['*' coreID '*.csv']));
    end
    if isempty(dlcMatch)
        SkippedTbl = [SkippedTbl; {string(coreID), moveFile, "DLC_Missing"}];
        continue;
    end


    try
        moveIndex = readtable(fullfile(vidDir, moveFile));
        moveIndex = moveIndex(moveIndex.BeginF > 0 & moveIndex.EndF > 0, :);
        rawDLC = readDLCwithReconstructedHeader(fullfile(csvDir, dlcMatch(1).name));
        rawDLC = cleanLowConfidence(rawDLC, 0.6);

        % Define anatomical label data to extract
        preferredMarkers = {'PalmBase', 'MCP1', 'fTip1','fTip2','fTip3','fTip4','fTip5'}; % update / change to Morgan H's PCs
        % preferredMarkers = {casePCs_MH.ID}; % top 6 PCs
        markerFound = false;

        % Initialize arrays to store X and Y positions for all markers
        allX = [];
        allY = [];

        % Loop through each preferred marker
        for k = 1:length(preferredMarkers)
            xVar = sprintf('%s_x', preferredMarkers{k});
            yVar = sprintf('%s_y', preferredMarkers{k});

            % Check if the marker exists in the data
            if all(ismember({xVar, yVar}, rawDLC.Properties.VariableNames))
                % Add the X and Y data for this marker to the arrays
                allX = [allX, rawDLC.(xVar)];
                allY = [allY, rawDLC.(yVar)];
                markerFound = true;
            end
        end
        % If no markers are found, skip this file
        % if ~markerFound, warning('[SKIP] No valid marker found in: %s', moveFile); continue; end
        if ~markerFound
            % warning('[SKIP] No valid marker found in: %s', moveFile);
            SkippedTbl = [SkippedTbl; {string(coreID), moveFile, "NoValidMarker"}]; continue;
        end


        % Compute the average X and Y positions across all preferred markers
        avgLabels_X = mean(allX, 2);  % Mean across columns (markers)
        avgLabels_Y = mean(allY, 2);  % Mean across columns (markers)

        if all(isnan(avgLabels_X)) || all(isnan(avgLabels_Y))
            warning('[SKIP] All marker data below confidence threshold in: %s', moveFile);
            SkippedTbl = [SkippedTbl; {string(coreID), moveFile, "AllLowConf"}];
            continue;
        end

        % Compute Euclidean distance between consecutive frames average X and Y dlc-label positions (in millimeters).
        dist_mm = sqrt(diff(avgLabels_X).^2 + diff(avgLabels_Y).^2) * px_to_mm;

        % Time vector
        time_sec = (1:length(dist_mm)) / fps;

        % Smooth movement distance data using Gaussian filter (for viz)
        smoothed_data = smoothdata(dist_mm, 'gaussian', 20);

        uniqueMoveTypes = unique(moveIndex.MoveType);
        for m = 1:numel(uniqueMoveTypes)
            moveType = uniqueMoveTypes{m};
            moveSubset = moveIndex(strcmpi(moveIndex.MoveType, moveType), :);

            allReps = table();  % store per-rep kinematic metrics
            repDurAll = []; interDurAll = []; velAll = []; ampAll = [];

            firstFrame = moveSubset.BeginF(1);
            lastFrame = min(moveSubset.EndF(end), length(smoothed_data));
            MovSeg = smoothed_data(firstFrame:lastFrame);

            % Time vector for this segment
            MoveSeg_t = time_sec(firstFrame:lastFrame);

            % Ensure MinPeakDistance is valid
            range_t = range(MoveSeg_t);
            if range_t < 0.1
                % warning('[SKIP] Movement segment too short in: %s - %s', coreID, moveType);
                SkippedTbl = [SkippedTbl; {string(coreID), moveFile, "TooShort"}];
                continue;
            end

            minDist_sec_default = min(0.12, 0.5 * range_t);
            minDist_sec = max(0.01, minDist_sec_default);  % ensure reasonable floor

            %  Set default parameters for findpeaks function
            minProm = std(MovSeg) * 0.4;
            minHeight = mean(MovSeg) + std(MovSeg) * 0.2;
            minDist = minDist_sec;

            normalizedMoveType = upper(strrep(moveType, ' ', ''));

            % Determine findpeaks parameters based on MoveType
            switch upper(strrep(moveType, ' ', ''))  % normalize input
                case {'HANDOC', 'HAND_OPEN_CLOSE'}
                    minProm = std(MovSeg) * 0.4;
                    minHeight = mean(MovSeg) + std(MovSeg) * 0.3;
                    minDist = minDist_sec;

                case {'HANDPS', 'HAND_PRONATE_SUPINATE'}
                    minProm = std(MovSeg) * 0.3;
                    minHeight = mean(MovSeg) + std(MovSeg) * 0.2;
                    minDist = minDist_sec * 0.75;

                case {'ARMEF', 'ARM_EXTENSION_FLEXION'}
                    minProm = std(MovSeg) * 0.2;
                    minHeight = mean(MovSeg) + std(MovSeg) * 0.1;
                    minDist = minDist_sec * 1.4;

                case {'REST'}
                    minProm = std(MovSeg) * 0.5;
                    minHeight = mean(MovSeg) + std(MovSeg) * 0.1;
                    minDist = minDist_sec;

                otherwise
                    warning('Unknown MoveType: %s. Using default thresholds.', moveType);
                    minProm = minProm;
                    minHeight = minHeight;
                    minDist = minDist_sec;
            end


            % Run peak detection: findpeaks
            [pks, locs, widths, proms] = findpeaks(MovSeg, MoveSeg_t, ...
                'MinPeakProminence', minProm, ...
                'MinPeakHeight', minHeight, ...
                'MinPeakDistance', minDist, Annotate ='extents');

            % if isempty(pks), continue; end
            if isempty(pks)
                % warning('[SKIP] No peaks found for: %s - %s', coreID, moveType);
                SkippedTbl = [SkippedTbl; {string(coreID), moveFile, "NoPeaksDetected"}];
                continue;
            end

            ampAll = [ampAll, pks(:)'];

            vel = pks ./ widths;  % mm/s
            velAll = [velAll, vel(:)'];

            % repDur = diff(locs) / fps;
            repDur   = widths;       % seconds, peak 'extent' as repetition duration
            repDurAll = [repDurAll, repDur(:)'];

            % interDur = repDur(1:end-1);
            interDur = diff(locs);   % seconds between successive peaks
            interDurAll = [interDurAll, interDur(:)'];

            % Plot once per MoveType per trial
            fig = figure('Visible','off');
            plot(MoveSeg_t, MovSeg, 'k'); hold on;
            plot(locs, pks, 'ro');
            xlabel('Time (s)'); ylabel('Displacement (mm)');
            title(sprintf('%s - %s', coreID, moveType));
            saveas(fig, fullfile(output_kinAnalysisDir, sprintf('%s_%s_peaks_summary.png', coreID, moveType)));
            close(fig);

            % Save per-repetition metrics
            for j = 1:height(moveSubset)
                beginF = moveSubset.BeginF(j); endF = min(moveSubset.EndF(j), length(smoothed_data));
                MovSeg_j = smoothed_data(beginF:endF);
                MoveSeg_tj = time_sec(beginF:endF);

                % Per-rep peak detection thresholds (same as segment-level)
                range_tj = range(MoveSeg_tj);
                minDist_sec_rep = min(0.10, 0.5 * range_tj);
                minDist_sec_rep = max(0.01, minDist_sec_rep);
                switch upper(moveType)
                    case 'HAND OC'
                        minProm = std(MovSeg_j) * 0.4;
                        % minHeight = mean(MovSeg_j) + std(MovSeg_j) * 0.3;
                        minDist = minDist_sec_rep;

                    case 'HAND PS'
                        minProm = std(MovSeg_j) * 0.3;
                        % minHeight = mean(MovSeg_j) + std(MovSeg_j) * 0.2;
                        minDist = minDist_sec_rep * 0.75;

                    case 'ARM EF'
                        minProm = std(MovSeg_j) * 0.2;
                        % minHeight = mean(MovSeg_j) + std(MovSeg_j) * 0.1;
                        minDist = minDist_sec_rep * 0.9;

                    case 'REST'
                        minProm = std(MovSeg_j) * 0.5;
                        % minHeight = mean(MovSeg_j) + std(MovSeg_j) * 0.4;
                        minDist = minDist_sec_rep;

                    otherwise
                        minProm = std(MovSeg_j) * 0.4;
                        % minHeight = mean(MovSeg_j) + std(MovSeg_j) * 0.2;
                        minDist = minDist_sec_rep;
                end

                % Run peak detection: findpeaks
                [pks_j, ~, widths_j, ~] = findpeaks(MovSeg_j, MoveSeg_tj, ...
                    'MinPeakProminence', minProm, ...
                    'MinPeakDistance', minDist);

                if isempty(pks_j), continue; end
                vel_j = pks_j ./ widths_j;  % mm/s
                valid_vel_j = vel_j(isfinite(vel_j) & ~isnan(vel_j));

                meanAmp = mean(pks_j, 'omitnan');
                stdAmp = std(pks_j, 'omitnan');
                meanVel = mean(valid_vel_j, 'omitnan');
                stdVel = std(valid_vel_j, 'omitnan');

                % Ensure scalar fallback
                if isempty(meanAmp), meanAmp = NaN; end
                if isempty(stdAmp), stdAmp = NaN; end
                if isempty(meanVel), meanVel = NaN; end
                if isempty(stdVel), stdVel = NaN; end

                kinTbl = [kinTbl; table(string(coreID), moveSubset.MoveN(j), string(moveType), ...
                    beginF, endF, meanAmp, stdAmp, meanVel, stdVel, ...
                    'VariableNames', {'TrialID','MoveN','MoveType','BeginF','EndF','MeanAmp','StdAmp','MeanVel','StdVel'})];

            end

            % Store summary metrics per MoveType
            kinSummaryTbl = [kinSummaryTbl; table(string(coreID), string(moveType), ...
                mean(ampAll,'omitnan'), std(ampAll,'omitnan'), ...
                mean(velAll,'omitnan'), std(velAll,'omitnan'), ...
                mean(repDurAll,'omitnan'), std(repDurAll,'omitnan'), ...
                mean(interDurAll,'omitnan'), std(interDurAll,'omitnan'), ...
                'VariableNames', {'TrialID','MoveType','MeanAmp','StdAmp','MeanVel','StdVel', ...
                'MeanRepDur','StdRepDur','MeanInterRepDur','StdInterRepDur'})];
        end

    catch ME
        warning('[ERROR] %s: %s', basePrefix, ME.message);
        continue;
    end
end

% Add Completion Message and Summary Trial Counts
% fprintf('[INFO] Movement feature analysis complete.\n');
fprintf('[INFO] %d movement trials summarized in kinTbl.\n', height(kinTbl));
fprintf('[INFO] %d trial-movement combinations summarized in kinSummaryTbl.\n', height(kinSummaryTbl));

% Add a Warning if Output Tables Are Still Empty
if isempty(kinTbl)
    warning('[WARN] kinTbl is empty. Check DLC marker quality or peak detection thresholds.');
end
if isempty(kinSummaryTbl)
    warning('[WARN] kinSummaryTbl is empty. May indicate all movements were skipped or failed peak detection.');
end

% Write output tables
writetable(kinTbl, fullfile(output_kinAnalysisDir, sprintf('kinTbl_%s.csv', MoveDir_CaseID)));
writetable(kinSummaryTbl, fullfile(output_kinAnalysisDir, sprintf('kinSummaryTbl_%s.csv', MoveDir_CaseID)));

% fprintf('[DEBUG] SkippedTbl has %d entries\n', height(SkippedTbl));
SkippedTbl = sortrows(SkippedTbl, {'Reason', 'TrialID'});
disp(SkippedTbl)
% if ~isempty(SkippedTbl)
%     writetable(SkippedTbl, fullfile(output_kinAnalysisDir, sprintf('SkippedTbl_%s.csv', MoveDir_CaseID)));
% end
SkippedTblPath = fullfile(output_kinAnalysisDir, sprintf('SkippedTbl_%s.csv', MoveDir_CaseID));
writetable(SkippedTbl, SkippedTblPath);
fprintf('[SAVED] SkippedTbl written to: %s\n', SkippedTblPath);


fprintf('[SAVED] All tables written to: %s\n', output_kinAnalysisDir);

end

%% Helper: DLC Header Reconstruction

function dlcTable = readDLCwithReconstructedHeader(dlcFilePath)
raw = readcell(dlcFilePath, 'NumHeaderLines', 0);
markerNames = raw(2,2:end); coordTypes = raw(3,2:end);
newColNames = strcat(markerNames, "_", coordTypes);
newColNames = matlab.lang.makeValidName(newColNames, 'ReplacementStyle','delete');
opts = detectImportOptions(dlcFilePath, 'NumHeaderLines', 3);
opts.VariableNames = [{'frame'}, newColNames];
dlcTable = readtable(dlcFilePath, opts);
end

%% Helper: Remove Low-Confidence Points

function dataOut = cleanLowConfidence(dataIn, threshold)
dataOut = dataIn;
vars = dataIn.Properties.VariableNames;
for i = 1:3:length(vars)
    if i+2 > numel(vars), break; end
    conf = dataIn{:, i+2};
    lowConf = conf < threshold;
    dataOut{lowConf, i} = NaN;
    dataOut{lowConf, i+1} = NaN;
end
end
