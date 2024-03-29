%% STATS functions - run on processed movement data

% Goal: Statistically assess movement timeseries data based on videos that have been anatomically labeled (13pt per frame) and analyzed via a trained DeepLabCut model

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

% Define switch case inputs
casedate = '09_12_2023';
hemisphere = 'R';

switch casedate
    case '09_12_2023'

        mainDir2 = [mainDir , filesep , '09_12_2023'];

    case '[insert relevant casedate]'

        mainDir2 = [mainDir , filesep , 'relevant casedate'];
end


switch hemisphere
    case 'L'

        mainDir3 = [mainDir2 , filesep , 'LSTN'];

    case 'R'

        mainDir3 = [mainDir2 , filesep , 'RSTN'];
end

cd(mainDir3)


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
outputDir = [mainDir3 filesep 'movementStats'];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end


%% Call summarized descriptive stats function for each condition

T3 = summarize_ConditionStats('OffOff', conditionsData.OffOff.amplitudes, conditionsData.OffOff.widths, conditionsData.OffOff.peakDists, outputDir);
T4 = summarize_ConditionStats('OffOn', conditionsData.OffOn.amplitudes, conditionsData.OffOn.widths, conditionsData.OffOn.peakDists, outputDir);
T5 = summarize_ConditionStats('OnOff', conditionsData.OnOff.amplitudes, conditionsData.OnOff.widths, conditionsData.OnOff.peakDists, outputDir);
T6 = summarize_ConditionStats('OnOn', conditionsData.OnOn.amplitudes, conditionsData.OnOn.widths, conditionsData.OnOn.peakDists, outputDir);

% Combine the tables vertically
T7 = [T3; T4; T5; T6];

% Save the combined table to a CSV file
writetable(T7, [outputDir filesep 'fTipTracking_results-per-condition_summary_mm_v2.csv'], 'WriteRowNames', true);


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


%% Call plotAnovaResults function for each measure

ANOVA_plot_allStates(all_amplitudes, group_labels_amplitudes, 'Amplitudes', 'Amplitude (mm)', 'Amplitude Comparison Across Conditions');
ANOVA_plot_allStates(all_widths, group_labels_widths, 'Intra-movement Durations', 'Intra-movement Durations (s)', 'Intra-movement Duration Comparison Across Conditions');
ANOVA_plot_allStates(all_peakDists, group_labels_peakDists, 'Inter-movement Durations', 'Inter-movement Durations (s)', 'Inter-movement Duration Comparison Across Conditions');

