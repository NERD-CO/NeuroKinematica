function [] = run_MovementStats_Clin_v2(mainDir, casedate_hem, hemisphere)

%% STATS functions - run on processed movement data

% Goal: Statistically summarize & compare movement timeseries data based on videos that have been anatomically labeled (13pt per frame) and analyzed via a trained DeepLabCut model

%% Directory set-up - Navigate b/t machines
pcname = getenv('COMPUTERNAME');

switch pcname
    case 'DESKTOP-I5CPDO7'   %%% JAT Desktop

        % mainDir = '';

    case 'DSKTP-JTLAB-EMR'   %%% ER Desktop

        mainDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\Kinematic Analyses';

    case 'NSG-M-FQBPFK3'     %%% ER PC

        mainDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\Kinematic Analyses';
end


%% Analyze data isolated by casedate and hemisphere

% Define casedate and hemisphere

 casedate_hem = '09_12_2023_LSTN';
% casedate_hem = '09_12_2023_RSTN';

mainDir2 = [mainDir , filesep , casedate_hem];
cd(mainDir2)

%% Isolate dlc outputs of interest

% Generate list of dlc-video-labeled CSV files
mainCSV = dir('*.csv');
mainCSV2 = {mainCSV.name};

% Generate list of dlc-video-labeled MAT files
mainMat = dir('*.mat');
mainMAT2 = {mainMat.name};

% Generate list of Motor Index CSVs (filters for CSVs that contain 'Move' string)
moveCSV = mainCSV2(contains(mainCSV2,'Move'));

%% Main function

% create an outputs directory
outputDir = [mainDir2 filesep 'movementStats'];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end


% Define framerate of videos (time conversion factor)
fps = 60; % frames per second

% Convert distance units to mm (distance conversion factor)
pixels_to_mm = 2.109; % 232 mm / 110 pxl = 2.1091 mm per pixel
% Anthropometry: vertical distance from the bottom of the chin (menton) to the top of the head: https://upload.wikimedia.org/wikipedia/commons/0/06/AvgHeadSizes.png
% US adult male, 50th percentile: Avg. = 23.2 cm, 9.1 inches
% Subject in video frames: Avg. = 110 pixels

% set up a results structure outside of loop
results = struct();

% initialize variables to store the mean and variability for each video:
mean_amplitudes = zeros(1, length(moveCSV));
std_amplitudes = zeros(1, length(moveCSV));
var_amplitudes = zeros(1, length(moveCSV));

mean_widths = zeros(1, length(moveCSV));
std_widths = zeros(1, length(moveCSV));
var_widths = zeros(1, length(moveCSV));

mean_peakDists = zeros(1, length(moveCSV));
std_peakDists = zeros(1, length(moveCSV));
var_peakDists = zeros(1, length(moveCSV));

% Initialize storage structures per condition
conditionsData = struct('OffOff', struct('amplitudes', [], 'widths', [], 'peakDists', []), ...
    'OffOn', struct('amplitudes', [], 'widths', [], 'peakDists', []), ...
    'OnOff', struct('amplitudes', [], 'widths', [], 'peakDists', []), ...
    'OnOn', struct('amplitudes', [], 'widths', [], 'peakDists', []));

% Initialize condition-based counters
countOffOff = 1;
countOffOn = 1;
% set ramp counter (countRamp = 1)
countOnOff = 1;
countOnOn = 1;


% Loop through CSV files
for csv_i = 1:length(moveCSV)

    tmpCSV = moveCSV{csv_i};

    % Split file names to extract relevant parts (dateID, sessID, and hemID)
    nameParts = split(tmpCSV,'_');
    dateID = nameParts{1};
    sessID = nameParts{3};
    hemID = nameParts{8};
    matName_title = [dateID , '-' , sessID, '-', hemID]

    % Find and load corresponding dlcDAT MAT file
    matTempfind = [dateID , '_' , sessID];
    matInd = contains(mainMAT2 , matTempfind);
    matName = mainMAT2{matInd};
    load(matName)

    % Process dlcDAT MAT file (all points, all frames) per vid first (Split column names of outDATA)
    colNames = outDATA.Properties.VariableNames; % outDATA should be a table containing labeled coordinate data from DeepLabCut
    colNames2 = cellfun(@(x) split(x,'_'), colNames,...
        'UniformOutput',false);
    colNames3 = unique(cellfun(@(x) x{1}, colNames2,...
        'UniformOutput',false));
    colNames4 = colNames3(~matches(colNames3,'frames'));

        % Initialize 'euclidall' to store Euclidean distances between successive points
    euclidall = zeros(height(outDATA)-1,length(colNames4));

    % smooth the position data first, then compute euclidean distances

    % Iterate over each label and compute Euclidean distance for each frame
    for label_i = 1:length(colNames4)

        tmpLabel_x = [colNames4{label_i} , '_x'];
        tmpLabel_y = [colNames4{label_i} , '_y'];

        tmpXdata = outDATA.(tmpLabel_x);
        tmpYdata = outDATA.(tmpLabel_y);

        labelData = [tmpXdata , tmpYdata];

        for frame_i = 1:height(labelData)
            if frame_i ~= height(labelData)
                point1 = labelData(frame_i,:);
                point2 = labelData(frame_i + 1,:);
                euclidall(frame_i , label_i) = pdist2(point1 , point2);
            end
        end
    end

    % Convert distance variables to mm usng conversion factor
    euclidall = euclidall * pixels_to_mm; % converting euclidean distances to mm

    % Filter the computed distances related to fingertip movements
    fTipInds = contains(colNames4,'fTip');
    fTipEuclid = euclidall(:,fTipInds);

    % Average the computed distances related to fingertip movements
    fTipAverage = mean(fTipEuclid,2);

    % Process dlcDAT MAT files using MoveIndex CSV files to select specific portions of the averaged fingertip distances
    moveINDtab = readtable(tmpCSV);
    moveINDtab = moveINDtab(~moveINDtab.BeginF == 0,:); % clean up - filters out rows in moveINDtab where the BeginF field is zero.
    moveINDtab = moveINDtab(~moveINDtab.EndF == 0,:); % clean up - filters out rows in moveINDtab where the EndF field is zero.

    % Align with Euclidean distance frames
    firstBegin = moveINDtab.BeginF(1) - 1; % assigned to value of the first element in the BeginF column of moveINDtab - 1 (because MATLAB uses 1-based indexing)
    lastEnd = moveINDtab.EndF(height(moveINDtab)) - 1; % assigned to value of the last element in the EndF column of moveINDtab - 1

    % Extract and store the average fingertip distance for specified frames (in fTip Average Block)
    fTipAveBlk = fTipAverage(firstBegin:lastEnd); % extracts subset of fTipAverage w/in frame range from firstBegin to lastEnd (represents specific portion of data where specified movement is detected, as indicated in MoveIndex CSV file)

    % Smooth out edges -- smoothdata function w/ 'guassian' method
    window_Width = 5; % set windowWidth as needed
    smoothed_fTipAveBlk = smoothdata(fTipAveBlk, 'gaussian', window_Width); % read documentation, window overlap

    % set thresholds based on mean +/- 3xStd
    MinPeakHeight = mean(smoothed_fTipAveBlk);
    MinPeakProminence = mean(smoothed_fTipAveBlk)*2;

    % Find peak amplitudes and compute widths -- findpeaks function [review documentation]
    [peaks_tmp, locs_tmp, widths_tmp, prominences_tmp] = findpeaks(smoothed_fTipAveBlk, fps, MinPeakHeight=mean(smoothed_fTipAveBlk), MinPeakDistance=0.15, MinPeakProminence=mean(smoothed_fTipAveBlk)*2, Annotate ='extents');
    [peaks, locs, widths, prominences] = findpeaks(smoothed_fTipAveBlk, MinPeakHeight=mean(smoothed_fTipAveBlk), MinPeakDistance=15, MinPeakProminence=mean(smoothed_fTipAveBlk)*2, Annotate ='extents');

    % Convert distance variables to mm usng distance conversion factor
    amplitudes = peaks * pixels_to_mm; % converting amplitudes to mm

    % Compute timepoints from locs (vector of integer indices corresponding to video frame number)
    timepoints = locs / fps; % Convert frame numbers to time (in seconds) using video sampling rate (Fs) conversion factor

    % Compute distances between consecutive peaks
    peakDists_frames = diff(locs); % by frame indice
    peakDists = diff(timepoints); % by timepoint (in seconds)

    % Convert frame-relative variables to seconds using time conversion factor
    widths_fps = widths / fps; % converting widths to seconds
    halfWidths = widths_fps / 2;
    % slope

    % normalize within patient - relative to bl (Off, off)

    % Compute timepoints for all indices
    timepoints__fTipAveBlk = (1:length(smoothed_fTipAveBlk))/fps;


    % Compute the mean and variability for each measurement
    statsData(csv_i).amplitudes = struct('mean', mean(amplitudes), 'std', std(amplitudes), 'var', var(amplitudes));
    statsData(csv_i).widths = struct('mean', mean(widths_fps), 'std', std(widths_fps), 'var', var(widths_fps));
    statsData(csv_i).peakDists = struct('mean', mean(peakDists), 'std', std(peakDists), 'var', var(peakDists));


    % Store the measures to the respective array based on the video index
    if hemisphere == 'L'
        % Off Med conditions
        if ismember(csv_i, [1, 2]) % L: sessions 1 & 3, R: sessions 2 & 4
            conditionsData.OffOff.amplitudes = [conditionsData.OffOff.amplitudes; amplitudes];
            conditionsData.OffOff.widths = [conditionsData.OffOff.widths; widths_fps];
            conditionsData.OffOff.peakDists = [conditionsData.OffOff.peakDists; peakDists];
        elseif ismember(csv_i, 3) % L: session 5, R: session 6 (stim ramping - fingertaps)
            % code for stim ramping - fingertaps
        elseif ismember(csv_i, [4, 5]) % L: sessions 7 & 9, R: sessions 8 & 10
            conditionsData.OffOn.amplitudes = [conditionsData.OffOn.amplitudes; amplitudes];
            conditionsData.OffOn.widths = [conditionsData.OffOn.widths; widths_fps];
            conditionsData.OffOn.peakDists = [conditionsData.OffOn.peakDists; peakDists];
            % On Med conditions
        elseif ismember(csv_i, [6, 7]) % L: sessions 13 & 15, R: sessions 14 & 17
            conditionsData.OnOff.amplitudes = [conditionsData.OnOff.amplitudes; amplitudes];
            conditionsData.OnOff.widths = [conditionsData.OnOff.widths; widths_fps];
            conditionsData.OnOff.peakDists = [conditionsData.OnOff.peakDists; peakDists];
        elseif ismember(csv_i, 8) % L: session 18, R: session 19 (stim ramping - fingertaps)
            % code for stim ramping - fingertaps
        elseif ismember(csv_i, [9, 10]) % L: sessions 20 & 22, R: sessions 21 & 23
            conditionsData.OnOn.amplitudes = [conditionsData.OnOn.amplitudes; amplitudes];
            conditionsData.OnOn.widths = [conditionsData.OnOn.widths; widths_fps];
            conditionsData.OnOn.peakDists = [conditionsData.OnOn.peakDists; peakDists];
        end

    elseif hemisphere == 'R'
        % Off Med conditions
        if ismember(csv_i, [1, 2]) % L: sessions 1 & 3, R: sessions 2 & 4
            conditionsData.OffOff.amplitudes = [conditionsData.OffOff.amplitudes; amplitudes];
            conditionsData.OffOff.widths = [conditionsData.OffOff.widths; widths_fps];
            conditionsData.OffOff.peakDists = [conditionsData.OffOff.peakDists; peakDists];
            %elseif ismember(csv_i, 3) % L: session 5, R: session 6 (stim ramping - fingertaps)
            % code for stim ramping - fingertaps
        elseif ismember(csv_i, [3, 4]) % L: sessions 7 & 9, R: sessions 8 & 10
            conditionsData.OffOn.amplitudes = [conditionsData.OffOn.amplitudes; amplitudes];
            conditionsData.OffOn.widths = [conditionsData.OffOn.widths; widths_fps];
            conditionsData.OffOn.peakDists = [conditionsData.OffOn.peakDists; peakDists];
            % On Med conditions
        elseif ismember(csv_i, [5, 6]) % L: sessions 13 & 15, R: sessions 14 & 17
            conditionsData.OnOff.amplitudes = [conditionsData.OnOff.amplitudes; amplitudes];
            conditionsData.OnOff.widths = [conditionsData.OnOff.widths; widths_fps];
            conditionsData.OnOff.peakDists = [conditionsData.OnOff.peakDists; peakDists];
            %elseif ismember(csv_i, 8) % L: session 18, R: session 19 (stim ramping - fingertaps)
            % code for stim ramping - fingertaps
        elseif ismember(csv_i, [7, 8]) % L: sessions 20 & 22, R: sessions 21 & 23
            conditionsData.OnOn.amplitudes = [conditionsData.OnOn.amplitudes; amplitudes];
            conditionsData.OnOn.widths = [conditionsData.OnOn.widths; widths_fps];
            conditionsData.OnOn.peakDists = [conditionsData.OnOn.peakDists; peakDists];
        end
    end

    % Loop to compute statistics for each condition
    conditionKeys = fieldnames(conditionsData); % Get all condition names
    for k = 1:length(conditionKeys)
        conditionName = conditionKeys{k};
        condition = conditionsData.(conditionName);

        % call compute_ConditionStats function
        conditionsData.(conditionName).stats.amplitudes = compute_ConditionStats(condition.amplitudes);
        conditionsData.(conditionName).stats.widths = compute_ConditionStats(condition.widths);
        conditionsData.(conditionName).stats.peakDists = compute_ConditionStats(condition.peakDists);
    end


    % Store results for each csv_i into a structure
    results(csv_i).fileName = tmpCSV;
    results(csv_i).matName = matName;
    results(csv_i).rawMovement = fTipAveBlk;
    results(csv_i).smoothMovement = smoothed_fTipAveBlk;
    results(csv_i).timepoints = timepoints;
    results(csv_i).locs = locs;
    results(csv_i).peaks = peaks;
    results(csv_i).amplitudes = amplitudes;
    results(csv_i).prominences = prominences;
    results(csv_i).widths = widths_fps;
    results(csv_i).halfWidths = halfWidths;
    results(csv_i).peakDists = peakDists;

    results(csv_i).mean_amplitudes = mean_amplitudes(csv_i);
    results(csv_i).std_amplitudes = std_amplitudes(csv_i);
    results(csv_i).var_amplitudes = var_amplitudes(csv_i);

    results(csv_i).mean_widths = mean_widths(csv_i);
    results(csv_i).std_widths = std_widths(csv_i);
    results(csv_i).var_widths = var_widths(csv_i);

    results(csv_i).mean_peakDists = mean_peakDists(csv_i);
    results(csv_i).std_peakDists = std_peakDists(csv_i);
    results(csv_i).var_peakDists = var_peakDists(csv_i);


    % Define unique name for the summary results output file based on the current CSV name
    descriptiveStats_summary = [outputDir filesep 'summary_output_' tmpCSV(1:end-44) '.csv']; % assumes tmpCSV is a string ending in '.csv'

    % Create a table to store summary results
    T2 = table(mean_amplitudes(csv_i), std_amplitudes(csv_i), var_amplitudes(csv_i), ...
        mean_widths(csv_i), std_widths(csv_i), var_widths(csv_i), ...
        mean_peakDists(csv_i), std_peakDists(csv_i), var_peakDists(csv_i), ...
        'VariableNames', {'MeanAmplitude', 'StdAmplitude', 'VarAmplitude', 'MeanWidth', 'StdWidth', 'VarWidth', 'MeanPeakDist', 'StdPeakDist', 'VarPeakDist'});

    % Write summary results table to a CSV file
    writetable(T2, descriptiveStats_summary);

end

% Save results structure in MAT file
save([outputDir filesep 'descriptiveStats_summary.mat'], 'results');


%% Plotting smooth movement comparisons

if hemisphere == 'L'
    % Concatenate smoothed fTip movement data for left hemisphere conditions
    concatenated_smoothed_OffOff = cat(1, results(1).smoothMovement, results(2).smoothMovement); % [OffMed, OffStim condition]
    concatenated_smoothed_OffOn = cat(1, results(4).smoothMovement, results(5).smoothMovement); % [OffMed, OnStim condition]
    concatenated_smoothed_OnOff = cat(1, results(6).smoothMovement, results(7).smoothMovement); % [OnMed, OffStim condition]
    concatenated_smoothed_OnOn = cat(1, results(9).smoothMovement, results(10).smoothMovement); % [OnMed, OnStim condition]
elseif hemisphere == 'R'
    % Concatenate smoothed fTip movement data for right hemisphere conditions
    concatenated_smoothed_OffOff = cat(1, results(1).smoothMovement, results(2).smoothMovement); % [OffMed, OffStim condition]
    concatenated_smoothed_OffOn = cat(1, results(3).smoothMovement, results(4).smoothMovement); % [OffMed, OnStim condition]
    concatenated_smoothed_OnOff = cat(1, results(5).smoothMovement, results(6).smoothMovement); % [OnMed, OffStim condition]
    concatenated_smoothed_OnOn = cat(1, results(7).smoothMovement, results(8).smoothMovement); % [OnMed, OnStim condition]
end

% Create a new time vector for concatenated data
timepoints_concatenated_OffOff = (1:length(concatenated_smoothed_OffOff))/fps;
timepoints_concatenated_OffOn = (1:length(concatenated_smoothed_OffOn))/fps;
timepoints_concatenated_OnOff = (1:length(concatenated_smoothed_OnOff))/fps;
timepoints_concatenated_OnOn = (1:length(concatenated_smoothed_OnOn))/fps;

% Plot the results
figure;

OffOff_color = [0.4 0.2 0.6]; % indigo/purple
OffOn_color = [0.2 0.7 0.8]; % teal/turquoise
OnOff_color = [0.8 0.3 0.1]; % orange
OnOn_color = [0.5 0.7 0.2]; % green

hold on;

% Plot concatenated smoothed data for videos [OffMed, OffStim condition]
plot(timepoints_concatenated_OffOff, concatenated_smoothed_OffOff, 'Color', OffOff_color);

% Plot concatenated smoothed data for videos [OnMed, OffStim condition]
plot(timepoints_concatenated_OffOn, concatenated_smoothed_OffOn, 'Color', OffOn_color);

% Plot concatenated smoothed data for videos [OffMed, OnStim condition]
plot(timepoints_concatenated_OnOff, concatenated_smoothed_OnOff, 'Color', OnOff_color);

% Plot concatenated smoothed data for videos [OnMed, OnStim condition]
plot(timepoints_concatenated_OnOn, concatenated_smoothed_OnOn, 'Color', OnOn_color);

xlabel('time (s)');
ylabel('finger movement amplitude');
legend('Off Med, Off Stim', 'Off Med, On Stim', 'Off Med, On Stim', 'On Med, On Stim');

hold off;


%% Call summarized descriptive stats function for each condition

OffOff = summarize_ConditionStats('OffOff', conditionsData.OffOff.amplitudes, conditionsData.OffOff.widths, conditionsData.OffOff.peakDists, outputDir);
OffOn = summarize_ConditionStats('OffOn', conditionsData.OffOn.amplitudes, conditionsData.OffOn.widths, conditionsData.OffOn.peakDists, outputDir);
OnOff = summarize_ConditionStats('OnOff', conditionsData.OnOff.amplitudes, conditionsData.OnOff.widths, conditionsData.OnOff.peakDists, outputDir);
OnOn = summarize_ConditionStats('OnOn', conditionsData.OnOn.amplitudes, conditionsData.OnOn.widths, conditionsData.OnOn.peakDists, outputDir);

% Combine the tables vertically
descriptiveStats_perCondition = [OffOff; OffOn; OnOff; OnOn];

% Save the combined table to a CSV file
writetable(descriptiveStats_perCondition, [outputDir filesep 'descriptiveStats-per-Condition_summary.csv'], 'WriteRowNames', true);


%% Computing stat comparisons (between 2 states)

% Compare results from OffMed to OnMed, OffMed to OffStim vs OnStim sessions - Hand OC
ttest_plot_2States_RandTrim('OffMed, OffStim', conditionsData.OffOff, 'OffMed, OnStim', conditionsData.OffOn, outputDir);
ttest_plot_2States_RandTrim('OnMed, OffStim', conditionsData.OnOff, 'OnMed, OnStim', conditionsData.OnOn, outputDir);


%% Computing stat comparisons (between all states) - ANOVA

% Ensure no NaN values
all_amplitudes = [conditionsData.OffOff.amplitudes; conditionsData.OffOn.amplitudes; conditionsData.OnOff.amplitudes; conditionsData.OnOn.amplitudes];
all_widths = [conditionsData.OffOff.widths; conditionsData.OffOn.widths; conditionsData.OnOff.widths; conditionsData.OnOn.widths];
all_peakDists = [conditionsData.OffOff.peakDists; conditionsData.OffOn.peakDists; conditionsData.OnOff.peakDists; conditionsData.OnOn.peakDists];

% remove NaN values or replace them
all_amplitudes = rmmissing(all_amplitudes);
all_widths = rmmissing(all_widths);
all_peakDists = rmmissing(all_peakDists);


% Determine minimum dataset size for amplitudes
minSize_amplitudes = min([length(conditionsData.OffOff.amplitudes), length(conditionsData.OffOn.amplitudes), ...
                    length(conditionsData.OnOff.amplitudes), length(conditionsData.OnOn.amplitudes)]);

% Randomly sample from each group to match the minimum group size
trimmed_OffOff_amplitudes = randsample(conditionsData.OffOff.amplitudes, minSize_amplitudes);
trimmed_OffOn_amplitudes = randsample(conditionsData.OffOn.amplitudes, minSize_amplitudes);
trimmed_OnOff_amplitudes = randsample(conditionsData.OnOff.amplitudes, minSize_amplitudes);
trimmed_OnOn_amplitudes = randsample(conditionsData.OnOn.amplitudes, minSize_amplitudes);

% Combine the trimmed groups
all_amplitudes = [trimmed_OffOff_amplitudes; trimmed_OffOn_amplitudes; ...
                  trimmed_OnOff_amplitudes; trimmed_OnOn_amplitudes];

% Create corresponding group labels (ensure data and group_labels are same length)
group_labels_amplitudes = [repmat({'OffOff'}, minSize_amplitudes, 1); ...
                           repmat({'OffOn'}, minSize_amplitudes, 1); ...
                           repmat({'OnOff'}, minSize_amplitudes, 1); ...
                           repmat({'OnOn'}, minSize_amplitudes, 1)];


% Determine minimum dataset size for widths
minSize_widths = min([length(conditionsData.OffOff.widths), length(conditionsData.OffOn.widths), ...
                    length(conditionsData.OnOff.widths), length(conditionsData.OnOn.widths)]);

% Randomly sample from each group to match the minimum group size
trimmed_OffOff_widths = randsample(conditionsData.OffOff.widths, minSize_widths);
trimmed_OffOn_widths = randsample(conditionsData.OffOn.widths, minSize_widths);
trimmed_OnOff_widths = randsample(conditionsData.OnOff.widths, minSize_widths);
trimmed_OnOn_widths = randsample(conditionsData.OnOn.widths, minSize_widths);

% Combine the trimmed groups
all_widths = [trimmed_OffOff_widths; trimmed_OffOn_widths; ...
              trimmed_OnOff_widths; trimmed_OnOn_widths];

% Create corresponding group labels (ensure data and group_labels are same length)
group_labels_widths = [repmat({'OffOff'}, minSize_widths, 1); ...
                        repmat({'OffOn'}, minSize_widths, 1); ...
                        repmat({'OnOff'}, minSize_widths, 1); ...
                        repmat({'OnOn'}, minSize_widths, 1)];


% Determine minimum dataset size for peakDists
minSize_peakDists = min([length(conditionsData.OffOff.peakDists), length(conditionsData.OffOn.peakDists), ...
                    length(conditionsData.OnOff.peakDists), length(conditionsData.OnOn.peakDists)]);

% Randomly sample from each group to match the minimum group size
trimmed_OffOff_peakDists = randsample(conditionsData.OffOff.peakDists, minSize_peakDists);
trimmed_OffOn_peakDists = randsample(conditionsData.OffOn.peakDists, minSize_peakDists);
trimmed_OnOff_peakDists = randsample(conditionsData.OnOff.peakDists, minSize_peakDists);
trimmed_OnOn_peakDists = randsample(conditionsData.OnOn.peakDists, minSize_peakDists);

% Combine the trimmed groups
all_peakDists = [trimmed_OffOff_peakDists; trimmed_OffOn_peakDists; ...
                 trimmed_OnOff_peakDists; trimmed_OnOn_peakDists];

% Create corresponding group labels (ensure data and group_labels are same length)
group_labels_peakDists = [repmat({'OffOff'}, minSize_peakDists, 1); ...
                           repmat({'OffOn'}, minSize_peakDists, 1); ...
                           repmat({'OnOff'}, minSize_peakDists, 1); ...
                           repmat({'OnOn'}, minSize_peakDists, 1)];


%% Call ANOVA_plot_allStates function for each measure

ANOVA_plot_allStates(all_amplitudes, group_labels_amplitudes, 'Amplitudes', 'Amplitude (mm)', 'Amplitude Comparison Across Conditions');
ANOVA_plot_allStates(all_widths, group_labels_widths, 'Intra-movement Durations', 'Intra-movement Durations (s)', 'Intra-movement Duration Comparison Across Conditions');
ANOVA_plot_allStates(all_peakDists, group_labels_peakDists, 'Inter-movement Durations', 'Inter-movement Durations (s)', 'Inter-movement Duration Comparison Across Conditions');

%%

cd(outputDir)

end
