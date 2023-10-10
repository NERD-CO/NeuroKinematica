%% fTip tracking script (v2)

% Goal: analyze  fingertip movement timeseries data based on videos that have been anatomically labeled (13pt per frame) and analyzed via a trained DeepLabCut model

%% Directory set-up - Navigate b/t machines
pcname = getenv('COMPUTERNAME');

switch pcname
    case 'DESKTOP-I5CPDO7'   %%% JAT

        mainDir = 'W:\RadcliffeE\INS_2024';

    case 'DSKTP-JTLAB-EMR'   %%% ER

        mainDir = 'Z:\RadcliffeE\INS_2024';

    case 'NSG-M-FQBPFK3'

        mainDir = 'Z:\RadcliffeE\INS_2024';
end


%% Analyze data isolated by hemisphere

hemisphere = 'L';

switch hemisphere
    case 'L'

        mainDir2 = [mainDir , filesep , 'LSTN'];

    case 'R'

        mainDir2 = [mainDir , filesep , 'LSTN'];
end

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
outputDir = [mainDir2 filesep 'fTipTracking_outputs_mm'];
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

% initialize arrays to store concatenated data for selected videos/conditions:
OffOff_amplitudes = [];
OffOff_widths = [];
OffOff_peakDists = [];

OffOn_amplitudes = [];
OffOn_widths = [];
OffOn_peakDists = [];

% initialize a figure outside of loop
figure;

% Loop through CSV files
for csv_i = 1:length(moveCSV)

    tmpCSV = moveCSV{csv_i};

    % Split file names to extract relevant parts (dateID, sessID, and hemID)
    nameParts = split(tmpCSV,'_');
    dateID = nameParts{1};
    sessID = nameParts{3};
    hemID = nameParts{8};
    matName_title = [dateID , '-' , sessID, '-', hemID];

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
    smoothed_fTipAveBlk = smoothdata(fTipAveBlk, 'gaussian', window_Width);

    % Find peak amplitudes and compute half-widths -- findpeaks function [review documentation]
    [peaks, locs, widths, prominences] = findpeaks(smoothed_fTipAveBlk, MinPeakHeight=20, MinPeakDistance=20, MinPeakProminence=10, Annotate ='extents'); 

    % Convert distance variables to mm usng conversion factor
    amplitudes = peaks * pixels_to_mm; % converting amplitudes to mm
    widths = widths * pixels_to_mm; % converting widths to mm
    halfWidths = widths / 2;

    % Compute timepoints from locs (vector of integer indices corresponding to video frame number)
    timepoints = locs / fps; % Convert frame numbers to time (in seconds) using video sampling rate (Fs) conversion factor

    % Compute distances between consecutive peaks
    peakDists = diff(locs); % by frame indice
    peakDists = diff(timepoints); % by timepoint

    % Compute timepoints for all indices
    timepoints__fTipAveBlk = (1:length(smoothed_fTipAveBlk))/fps;

    % Plot smooth movement for each CSV iteration
    subplot(length(moveCSV), 1, csv_i);
    hold on
    % Adjust the parameters in findpeaks
    findpeaks(smoothed_fTipAveBlk, timepoints__fTipAveBlk, MinPeakHeight=20, MinPeakDistance=0.20, MinPeakProminence=10, Annotate ='extents');
    % define axes labels and subplot titles
    xlabel('time (s)');
    ylabel('amplitude');
    hold off
    title(['Smooth Hand O/C Movement, ', num2str(matName_title)])


    % Compute the mean and variability for each measurement
    mean_amplitudes(csv_i) = mean(amplitudes);
    std_amplitudes(csv_i) = std(amplitudes);
    var_amplitudes(csv_i) = var(amplitudes);

    mean_widths(csv_i) = mean(widths);
    std_widths(csv_i) = std(widths);
    var_widths(csv_i) = var(widths);

    mean_peakDists(csv_i) = mean(peakDists);
    std_peakDists(csv_i) = std(peakDists);
    var_peakDists(csv_i) = var(peakDists);

    % Store the measures to the respective array based on the video index
    if ismember(csv_i, [1, 2]) % sessions 1 & 3
        OffOff_amplitudes = [OffOff_amplitudes; amplitudes];
        OffOff_widths = [OffOff_widths; widths];
        OffOff_peakDists = [OffOff_peakDists; peakDists];
    elseif ismember(csv_i, [4, 5]) % sessions 7 & 9
        OffOn_amplitudes = [OffOn_amplitudes; amplitudes];
        OffOn_widths = [OffOn_widths; widths];
        OffOn_peakDists = [OffOn_peakDists; peakDists];
    end

    % Compute the mean and variability for each measurement in the OffMed, OffStim condition
    mean_OffOff_amplitudes_i(csv_i) = mean(OffOff_amplitudes);
    std_OffOff_amplitudes_i(csv_i) = std(OffOff_amplitudes);
    var_OffOff_amplitudes_i(csv_i) = var(OffOff_amplitudes);

    mean_OffOff_widths_i(csv_i) = mean(OffOff_widths);
    std_OffOff_widths_i(csv_i) = std(OffOff_widths);
    var_OffOff_widths_i(csv_i) = var(OffOff_widths);

    mean_OffOff_peakDists_i(csv_i) = mean(OffOff_peakDists);
    std_OffOff_peakDists_i(csv_i) = std(OffOff_peakDists);
    var_OffOff_peakDists_i(csv_i) = var(OffOff_peakDists);

    % Compute the mean and variability for each measurement in the OffMed, OnStim condition
    mean_OffOn_amplitudes_i(csv_i) = mean(OffOn_amplitudes);
    std_OffOn_amplitudes_i(csv_i) = std(OffOn_amplitudes);
    var_OffOn_amplitudes_i(csv_i) = var(OffOn_amplitudes);

    mean_OffOn_widths_i(csv_i) = mean(OffOn_widths);
    std_OffOn_widths_i(csv_i) = std(OffOn_widths);
    var_OffOn_widths_i(csv_i) = var(OffOn_widths);

    mean_OffOn_peakDists_i(csv_i) = mean(OffOn_peakDists);
    std_OffOn_peakDists_i(csv_i) = std(OffOn_peakDists);
    var_OffOn_peakDists_i(csv_i) = var(OffOn_peakDists);


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
    results(csv_i).widths = widths;
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


    % Define unique name for the results output file based on the current CSV name
    fTipTracking_results_mm = [outputDir filesep 'output_' tmpCSV(1:end-44) '.csv']; % assumes tmpCSV is a string ending in '.csv'

    % Create a table to store results based on computed variables
    T1 = table(timepoints, locs, peaks, amplitudes, prominences, widths, halfWidths, 'VariableNames', {'Timepoints', 'Locations', 'Peaks', 'Amplitudes', 'Prominences', 'Widths', 'HalfWidths'});

    % Write results table to a CSV file
    writetable(T1, fTipTracking_results_mm);


    % Define unique name for the summary results output file based on the current CSV name
    fTipTracking_results_summary_mm = [outputDir filesep 'summary_output_' tmpCSV(1:end-44) '.csv']; % assumes tmpCSV is a string ending in '.csv'

    % Create a table to store summary results
    T2 = table(mean_amplitudes(csv_i), std_amplitudes(csv_i), var_amplitudes(csv_i), ...
        mean_widths(csv_i), std_widths(csv_i), var_widths(csv_i), ...
        mean_peakDists(csv_i), std_peakDists(csv_i), var_peakDists(csv_i), ...
        'VariableNames', {'MeanAmplitude', 'StdAmplitude', 'VarAmplitude', 'MeanWidth', 'StdWidth', 'VarWidth', 'MeanPeakDist', 'StdPeakDist', 'VarPeakDist'});

    % Write summary results table to a CSV file
    writetable(T2, fTipTracking_results_summary_mm);


end

% Save results structure in MAT file
save([outputDir filesep 'fTipTracking_results_mm.mat'], 'results');


%% Compute and save concatenated summary results for analyses

% Compute overall mean and variability for OffMed, OffStim condition (sessions 1 & 3)
mean_OffOff_amplitudes = mean(OffOff_amplitudes);
std_OffOff_amplitudes = std(OffOff_amplitudes);
var_OffOff_amplitudes = var(OffOff_amplitudes);

mean_OffOff_widths = mean(OffOff_widths);
std_OffOff_widths = std(OffOff_widths);
var_OffOff_widths = var(OffOff_widths);

mean_OffOff_peakDists = mean(OffOff_peakDists);
std_OffOff_peakDists = std(OffOff_peakDists);
var_OffOff_peakDists = var(OffOff_peakDists);

% Create a table to store overall summary results for OffMed, OffStim condition
T3 = table(mean_OffOff_amplitudes, std_OffOff_amplitudes, var_OffOff_amplitudes, ...
    mean_OffOff_widths, std_OffOff_widths, var_OffOff_widths, ...
    mean_OffOff_peakDists, std_OffOff_peakDists, var_OffOff_peakDists, ...
    'VariableNames', {'MeanAmplitude', 'StdAmplitude', 'VarAmplitude', 'MeanWidth', 'StdWidth', 'VarWidth', 'MeanPeakDist', 'StdPeakDist', 'VarPeakDist'});

% Write overall summary results table to a CSV file for OffMed, OffStim condition
writetable(T3, [outputDir filesep 'fTipTracking_results_OffOff_summary_mm.csv']);


% Compute overall mean and variability for OffMed, OnStim condition (sessions 7 & 9)
mean_OffOn_amplitudes = mean(OffOn_amplitudes);
std_OffOn_amplitudes = std(OffOn_amplitudes);
var_OffOn_amplitudes = var(OffOn_amplitudes);

mean_OffOn_widths = mean(OffOn_widths);
std_OffOn_widths = std(OffOn_widths);
var_OffOn_widths = var(OffOn_widths);

mean_OffOn_peakDists = mean(OffOn_peakDists);
std_OffOn_peakDists = std(OffOn_peakDists);
var_OffOn_peakDists = var(OffOn_peakDists);

% Create a table to store overall summary results for OffMed, OffStim condition
T4 = table(mean_OffOn_amplitudes, std_OffOn_amplitudes, var_OffOn_amplitudes, ...
    mean_OffOn_widths, std_OffOn_widths, var_OffOn_widths, ...
    mean_OffOn_peakDists, std_OffOn_peakDists, var_OffOn_peakDists, ...
    'VariableNames', {'MeanAmplitude', 'StdAmplitude', 'VarAmplitude', 'MeanWidth', 'StdWidth', 'VarWidth', 'MeanPeakDist', 'StdPeakDist', 'VarPeakDist'});

% Write overall summary results table to a CSV file for OffMed, OffStim condition
writetable(T4, [outputDir filesep 'fTipTracking_results_OffOn_summary_mm.csv']);


% Assign row names to T3 (OffOff) and T4 (OffOn) summary results tables
T3.Properties.RowNames = {'OffOff_Condition'};
T4.Properties.RowNames = {'OffOn_Condition'};

% Combine the tables vertically
T5 = [T3; T4];

% Save the combined table to a CSV file
writetable(T5, [outputDir filesep 'fTipTracking_results-per-condition_summary_mm.csv'], 'WriteRowNames', true);


%% Plotting smooth movement comparisions

% Concatenate smoothed fTip movement data for videos 1 and 2 [OffMed, OffStim condition (sessions 1 & 3)]
concatenated_smoothed_OffOff = cat(1, results(1).smoothMovement, results(2).smoothMovement);

% Concatenate smoothed fTip movement data for videos 4 and 5 [OffMed, OnStim condition (sessions 7 & 9)]
concatenated_smoothed_OffOn = cat(1, results(4).smoothMovement, results(5).smoothMovement);

% Create a new time vector for concatenated data
timepoints_concatenated_OffOff = (1:length(concatenated_smoothed_OffOff))/fps;
timepoints_concatenated_OffOn = (1:length(concatenated_smoothed_OffOn))/fps;

% Plot the results
figure;

OffOffcolor = [0.4 0.2 0.6]; % indigo/purple
OffOncolor = [0.2 0.7 0.8]; % teal/turquoise

hold on;

% Plot concatenated smoothed data for videos 1 and 2
plot(timepoints_concatenated_OffOff, concatenated_smoothed_OffOff, 'Color', OffOffcolor);

% Plot concatenated smoothed data for videos 4 and 5
plot(timepoints_concatenated_OffOn, concatenated_smoothed_OffOn, 'Color', OffOncolor);

xlabel('time (s)');
ylabel('finger movement amplitude');
legend('Off Med, Off Stim', 'Off Med, On Stim');

hold off;

%% Computing stat comparisons

% compare results from sessions 1 & 3 (Off Med, Off Stim) to sessions 7 & 9 (Off Med, On Stim) - Hand OC

% Perform t-tests
[~, p_value_amplitude] = ttest2(OffOff_amplitudes, OffOn_amplitudes);
[~, p_value_width] = ttest2(OffOff_widths, OffOn_widths);
[~, p_value_peakDists] = ttest2(OffOff_peakDists, OffOn_peakDists);

% Display p-values
disp(['p-value for amplitude: ', num2str(p_value_amplitude)]);
disp(['p-value for width: ', num2str(p_value_width)]);
disp(['p-value for peak distance: ', num2str(p_value_peakDists)]);

% Create p-value table
p_value_table = table(p_value_amplitude, p_value_width, p_value_peakDists, 'VariableNames', {'p_val_Amplitude', 'p_val_Width', 'p_val_PeakDistance'});

% Save the table to a CSV file
writetable(p_value_table, [outputDir filesep 'p_values_mm.csv']);

%% Plotting stat comparisons

% Input data
OffOff_means = [mean_OffOff_amplitudes, mean_OffOff_widths, mean_OffOff_peakDists];
OffOn_means = [mean_OffOn_amplitudes, mean_OffOn_widths, mean_OffOn_peakDists];

OffOff_sd = [std_OffOff_amplitudes, std_OffOff_widths, std_OffOff_peakDists];
OffOn_sd = [std_OffOn_amplitudes, std_OffOn_widths, std_OffOn_peakDists];

OffOff_var = [var_OffOff_amplitudes, var_OffOff_widths, var_OffOff_peakDists];
OffOn_var = [var_OffOn_amplitudes, var_OffOn_widths, var_OffOn_peakDists];

% Data organization for plotting
means_matrix = [OffOff_means; OffOn_means];
sd_matrix = [OffOff_sd; OffOn_sd];
var_matrix = [OffOff_var; OffOn_var];

% Labels
group_labels = {'Off Med, Off Stim', 'Off Med, On Stim'};
metric_labels = {'Peak Amplitudes', 'Peak Widths', 'Peak Distances'};

% Plotting
figure;
colors = [0.4 0.2 0.6; 0.2 0.7 0.8]; % RGB triplet per condition bar

for i = 1:3
    subplot(1, 3, i);
    b = bar(means_matrix(:, i), 'FaceColor', 'flat');
    hold on;

    % Assign colors per condition bar
    for cond = 1:2
        b.CData(cond,:) = colors(cond, :);
    end

    % Add text annotations
    for cond = 1:2
        x = b.XEndPoints(cond); % get the x-coordinate of each bar
        y = b.YEndPoints(cond); % get the y-coordinate of each bar
        text(x, y + 0.05, sprintf('%.2f (SD=%.2f)', means_matrix(cond,i), sd_matrix(cond,i)),...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end

    % Use either standard deviation (sd_matrix) or variance (var_matrix) for the error bars
    errorbar(1:2, means_matrix(:, i), sd_matrix(:, i), '.k', 'LineWidth', 1.5);

    set(gca, 'XTickLabel', group_labels);
    title(['Comparison of ', metric_labels{i}]);
    ylabel(['Mean ', metric_labels{i}]);
    xlabel('Condition');

    hold off;
end


%% compare set results (6) within session 5 (Off Med, Ramp Stim) - Finger Taps

