function [kinTbl, RepSummaryTbl] = run_MovementFeatureAnalysis_IO(CaseDate, CaseDate_hem, MoveDataDir, MoveDir_CaseID, fps, px_to_mm, zscoreKin)

% run_MovementFeatureAnalysis_IO
% Extract movement features from DLC-labeled kinematic CSVs using Movement Index files
% Computes per-trial kinematic features for intraoperative motor recordings.

% Inputs:
% - CaseDate: e.g. '03_23_2023'
% - CaseDate_hem: 'LSTN' or 'RSTN'
% - MoveDataDir: full path to top-level 'Processed DLC' directory
% - MoveDir_CaseID: subfolder with movement CSVs, e.g., 'IO_03_23_2023_LSTN'
% - fps: video framerate (e.g. 100)
% - px_to_mm: scale conversion (e.g. 2.109 mm/px)
% - zscoreKin: true/false toggle to z-score output features

% Outputs:
%   kinTbl - trial-wise feature table with MeanAmp, StdAmp, MeanVel, StdVel
%   RepSummaryTbl - summary of MeanRepDur & MeanInterRepDur per MoveType and trial ID

% close all; clc;
clearvars -except CaseDate CaseDate_hem ephys_offset; clc;

%% Directory Setup

curPC = getenv('COMPUTERNAME');
switch curPC
    case 'DESKTOP-I5CPDO7'
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    otherwise
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
end

% Define MoveData dirs
MoveDir_CaseID = 'IO_03_23_2023_LSTN';
MoveDataDir = fullfile(IO_DataDir, 'Processed DLC');
fprintf('[INFO] Running movement feature analysis for: %s\n', MoveDir_CaseID);

% Define MoveData subdirs
caseDir = fullfile(MoveDataDir, MoveDir_CaseID);
vidDir  = fullfile(caseDir, 'video folder');
csvDir  = fullfile(caseDir, 'csv folder');

moveFiles = dir(fullfile(vidDir, '*Move*.csv'));
fprintf('[INFO] Found %d Move Index files\n', numel(moveFiles));


%% Define output dir

output_kinAnalysisDir = fullfile(IO_DataDir, 'Kinematic Analyses', MoveDir_CaseID);
if ~exist(output_kinAnalysisDir, 'dir')
    mkdir(output_kinAnalysisDir);
end

% Initialize outputs and storage
kinTbl = table();
RepSummaryTbl = table();   % For storing per trial+type repetition durations
SkippedTbl = table();      % For logging short/failed peak trials

%% Sampling rates and Calibration

AO_spike_fs = 44000; % MER sampling rate
fps = 100; % IO video frame rate

% calibration / conversion metric
px_to_mm      = 2.109; % distance conversion factor

%% update to calibration-image based method later

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
    [~, basePrefix, ~] = fileparts(moveFile);
    tokens = regexp(basePrefix, '(20\d{6}_[bct]\d+_d\d+p\d+_session\d+)', 'tokens');
    if isempty(tokens), continue; end
    coreID = tokens{1}{1};

    dlcMatch = dir(fullfile(csvDir, [coreID, '*DLC*.csv']));
    if isempty(dlcMatch)
        warning('[SKIP] DLC CSV not found for %s', basePrefix);
        continue;
    end

    try
        moveIndex = readtable(fullfile(vidDir, moveFile));
        moveIndex = moveIndex(moveIndex.BeginF > 0 & moveIndex.EndF > 0, :);

        rawDLC = readDLCwithReconstructedHeader(fullfile(csvDir, dlcMatch(1).name));
        rawDLC = cleanLowConfidence(rawDLC, 0.6);

        preferredMarkers = {'PalmBase', 'fTip1','fTip2','fTip3','fTip4','fTip5'};
        markerFound = false;
        for k = 1:length(preferredMarkers)
            xVar = sprintf('%s_x', preferredMarkers{k});
            yVar = sprintf('%s_y', preferredMarkers{k});
            if all(ismember({xVar, yVar}, rawDLC.Properties.VariableNames))
                X = rawDLC.(xVar); Y = rawDLC.(yVar);
                markerUsed = preferredMarkers{k}; markerFound = true; break;
            end
        end

        if ~markerFound
            warning('[SKIP] No valid marker found in: %s', moveFile);
            continue;
        end

        % Compute smoothed_data after confirming marker is found
        dist_mm = sqrt(diff(X).^2 + diff(Y).^2) * px_to_mm;
        time_sec = (1:length(dist_mm)) / fps;
        smoothed_data = smoothdata(dist_mm, 'gaussian', 15);

    catch ME
        warning('[ERROR] Problem processing %s: %s', moveFile, ME.message);
        continue;
    end


    % Per-repetition loop
    repDurAll = [];
    interDurAll = [];
    storedPlot = false;  % flag to save only one peak plot per MoveType

    for j = 1:height(moveIndex)
        beginF = moveIndex.BeginF(j);
        endF = min(moveIndex.EndF(j), length(smoothed_data));
        MovSeg = smoothed_data(beginF:endF);
        MoveSeg_t = time_sec(beginF:endF);

        % Pre-check: Skip short segments before calling findpeaks
        if length(MovSeg) < 3  % adjust threshold as needed
            SkippedTbl = [SkippedTbl; table(string(coreID), moveIndex.MoveN(j), string(moveIndex.MoveType{j}), ...
                sprintf("Segment too short: %d frames", length(MovSeg)), ...
                'VariableNames', {'TrialID','MoveN','MoveType','Reason'})];
            continue;
        end

        % Compute MinPeakDistance in seconds
        minDist = max(1, min(round(0.15 * fps), length(MovSeg)-1));
        minDist_sec = minDist / fps;

        % Check if time span is too short for peak detection
        if MoveSeg_t(end) - MoveSeg_t(1) <= minDist_sec
            SkippedTbl = [SkippedTbl; table(string(coreID), moveIndex.MoveN(j), string(moveIndex.MoveType{j}), ...
                sprintf("MinPeakDistance too large (%.3f sec) for segment time span (%.3f sec)", ...
                minDist_sec, MoveSeg_t(end) - MoveSeg_t(1)), ...
                'VariableNames', {'TrialID','MoveN','MoveType','Reason'})];
            continue;
        end


        % Peak detection: findpeaks
        [pks, locs, widths, proms] = findpeaks(MovSeg, MoveSeg_t, ...
            'MinPeakProminence', std(MovSeg)*0.5, ...
            'MinPeakDistance', minDist_sec);


        % Post-check: Skip segments without peaks
        if isempty(pks)
            SkippedTbl = [SkippedTbl; table(string(coreID), moveIndex.MoveN(j), string(moveIndex.MoveType{j}), ...
                "Too few peaks", 'VariableNames', {'TrialID','MoveN','MoveType','Reason'})];
            continue;
        end


        % Kinematic feature calc
        if isempty(pks)
            ampMean = NaN;
            ampStd = NaN;
        else
            ampMean = mean(pks);
            ampStd = std(pks);
        end

        % Ensure vel is a numeric vector with finite values only
        vel = pks ./ (widths / fps);
        valid_vel = vel(isfinite(vel) & ~isnan(vel));
        if isempty(valid_vel)
            velMean = NaN;
            velStd = NaN;
        else
            velMean = mean(valid_vel);
            velStd = std(valid_vel);
        end
        repDur = diff(locs) / fps;
        interDur = repDur(1:end-1);

        % % Force all variables being inserted into kinTbl to be scalars, with fallback defaults if empty
        % if isempty(ampMean), ampMean = NaN; end
        % if isempty(ampStd), ampStd = NaN; end
        % if isempty(velMean), velMean = NaN; end
        % if isempty(velStd), velStd = NaN; end

        disp([ampMean, ampStd, velMean, velStd])

        % Store trial-wise metrics
        kinTbl = [kinTbl; ...
            table(string(coreID), moveIndex.MoveN(j), string(moveIndex.MoveType{j}), beginF, endF, ...
            ampMean, ampStd, velMean, velStd, ...
            'VariableNames', {'TrialID','MoveN','MoveType','BeginF','EndF','MeanAmp','StdAmp','MeanVel','StdVel'})];

        % Accumulate rep/inter-rep durations
        repDurAll = [repDurAll, repDur(:)'];
        interDurAll = [interDurAll, interDur(:)'];

        % Save a plot only once per MoveType per trial
        if ~storedPlot
            f = figure('Visible','off'); hold on
            plot(MoveSeg_t, MovSeg, 'k');
            plot(locs, pks, 'ro');  % locs already contains time values when x-data is passed to findpeaks
            xlabel('Time (s)'); ylabel('Displacement (mm)');
            title(sprintf('%s - %s', coreID, moveIndex.MoveType{j}));
            saveas(f, fullfile(output_kinAnalysisDir, sprintf('%s_%s_peaks_summary.png', ...
                coreID, moveIndex.MoveType{j})));
            close(f);
            storedPlot = true;
        end
    end

    % After full trial loop, store 1 row per (TrialID, MoveType)
    if ~isempty(repDurAll)
        RepSummaryTbl = [RepSummaryTbl; ...
            table(string(coreID), string(moveIndex.MoveType{1}), ...
            mean(repDurAll, 'omitnan'), mean(interDurAll, 'omitnan'), ...
            'VariableNames', {'TrialID','MoveType','MeanRepDur','MeanInterRepDur'})];
    end

    %% Format outputs
    fprintf('[DONE] Processed %d entries into kinTbl.\n', height(kinTbl));

    % Save outputs
    writetable(kinTbl, fullfile(output_kinAnalysisDir, sprintf('kinTbl_%s.csv', MoveDir_CaseID)));
    writetable(RepSummaryTbl, fullfile(output_kinAnalysisDir, sprintf('RepSummaryTbl_%s.csv', MoveDir_CaseID)));

    if ~isempty(SkippedTbl)
        writetable(SkippedTbl, fullfile(output_kinAnalysisDir, sprintf('SkippedTbl_%s.csv', MoveDir_CaseID)));
    end
    fprintf('[SAVED] All feature tables written to: %s\n', output_kinAnalysisDir);

end

end

% ======= Helper Functions =======

%% Subfunction: Reconstruct DLC Header
function dlcTable = readDLCwithReconstructedHeader(dlcFilePath)
raw = readcell(dlcFilePath, 'NumHeaderLines', 0);
markerNames = raw(2,2:end); coordTypes = raw(3,2:end);
newColNames = strcat(markerNames, "_", coordTypes);
newColNames = matlab.lang.makeValidName(newColNames, 'ReplacementStyle','delete');
opts = detectImportOptions(dlcFilePath, 'NumHeaderLines', 3);
opts.VariableNames = [{'frame'}, newColNames];
dlcTable = readtable(dlcFilePath, opts);
end

%% Subfunction: Clean low-confidence DLC values
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
