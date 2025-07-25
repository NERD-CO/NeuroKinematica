function IO_plotting_v4(CaseDate, CaseDate_hem, ephys_offset)

% IO_plotting_v4: Generates Raster, PSTH, and FR plots. and summary CSV exports per depth & movement type

clear; clc;
addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis'


%% Inputs: 

CaseDate = '05_18_2023_b_bilateral'; %% adjust as needed
% '03_23_2023'; % unilateral example
% '05_18_2023_b_bilateral'; bilateral example


ephys_offset = 1;

% Example:
% IO_plotting_v4('03_23_2023', '', 1);
% IO_plotting_v4('05_18_2023_b_bilateral', 'LSTN', 1);


%% Directory Setup

curPCname = getenv('COMPUTERNAME');
switch curPCname
    case 'DESKTOP-I5CPDO7'
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    otherwise
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
end


%% CONFIG - script parameters, AO_fs, bin sizes

useSEM = true; % true = use SEM, false = use STD
AO_spike_fs = 44000;

% FR bins (10 ms)
binSize_FR = 0.01; 
window_FR  = [-0.05 0.45];
edges_FR   = window_FR(1):binSize_FR:window_FR(2);
time_FR    = edges_FR(1:end-1)+binSize_FR/2;

% PSTH bins (50 ms)
binSize_PSTH = 0.05;
edges_PSTH   = window_FR(1):binSize_PSTH:window_FR(2);
time_PSTH    = edges_PSTH(1:end-1)+binSize_PSTH/2;

depth_ids = {'t','c','b'}; % dorsal, central, ventral

% dynamic error label
if useSEM
    errLabel = 'SEM';
else
    errLabel = 'STD';
end


%% Define Case-specific Parameters

Case_DataDir = fullfile(IO_DataDir, CaseDate);
ephysTbl_Dir = fullfile(Case_DataDir, 'DLC_Ephys'); % Base ephys directory

% Handle bilateral cases and hemisphere selection
isBilateral = contains(CaseDate, 'bilateral', 'IgnoreCase', true);

if isBilateral
    fprintf('[INFO] Bilateral case detected: %s\n', CaseDate);
    
    % Prompt user for hemisphere choice (LSTN or RSTN)
    CaseDate_hem = input('Enter hemisphere (LSTN or RSTN): ', 's');
    
    % Validate input
    validHems = {'LSTN','RSTN'};
    if ~ismember(CaseDate_hem, validHems)
        error('Invalid input. Please enter LSTN or RSTN.');
    end
    
    % Append hemisphere-specific folder
    ephysTbl_Dir = fullfile(ephysTbl_Dir, CaseDate_hem);
    fprintf('[INFO] Hemisphere-specific case directory: %s\n', ephysTbl_Dir);
else
    CaseDate_hem = ''; % No hemisphere for unilateral cases
    fprintf('[INFO] Unilateral case directory: %s\n', ephysTbl_Dir);
end

%% Load All_SpikesPerMove_Tbl

Tbl_list = dir('*Spikes*.mat');
Tbl_names = {Tbl_list.name};

if ephys_offset
    spk_case = Tbl_names{contains(Tbl_names, 'Spikes') & contains(Tbl_names, 'offset')};
else
    spk_case = Tbl_names{contains(Tbl_names, 'Spikes') & ~contains(Tbl_names, 'offset')};
end

if isempty(spk_case)
    error('No matching spike table found in %s', ephysTbl_Dir);
end

load(spk_case, 'All_SpikesPerMove_Tbl');
move_types = unique(All_SpikesPerMove_Tbl.MoveType);


%% Output Directory

figBaseDir = fullfile(Case_DataDir, 'DLC_Ephys');
if ~isempty(CaseDate_hem)
    figBaseDir = fullfile(figBaseDir, CaseDate_hem, 'Figures', 'Raster_PSTH');
else 
    figBaseDir = fullfile(figBaseDir, 'Figures', 'Raster_PSTH');
end

if ~exist(figBaseDir, 'dir')
    mkdir(figBaseDir);
end
fprintf('[INFO] Figures will be saved in: %s\n', figBaseDir);

csvDir = fullfile(Case_DataDir, 'DLC_Ephys');
if ~isempty(CaseDate_hem)
    csvDir = fullfile(csvDir, CaseDate_hem, 'CSV_Outputs');
else
    csvDir = fullfile(csvDir, 'CSV_Outputs');
end
if ~exist(csvDir, 'dir'), mkdir(csvDir); end
fprintf('[INFO] CSVs will be saved in: %s\n', csvDir);


%% Loop Through Depths and Movement Types

% Initialize storage for CSV
summary_Data = {};
full_Data = {};

for d = 1:numel(depth_ids)
    depth_code = depth_ids{d};
    depth_name = switchDepth(depth_code);

    depth_tbl = All_SpikesPerMove_Tbl(contains(All_SpikesPerMove_Tbl.move_trial_ID, depth_code), :);
    if isempty(depth_tbl)
        continue;
    end
    
    fprintf('\n[INFO] Processing depth: %s (%s)\n', depth_name, depth_code);
    figDepthDir = fullfile(figBaseDir, depth_name);
    if ~exist(figDepthDir, 'dir'), mkdir(figDepthDir); end
    
    for m = 1:numel(move_types)
        move_type = move_types{m};
        move_tbl = depth_tbl(strcmp(depth_tbl.MoveType, move_type), :);
        if isempty(move_tbl)
            continue;
        end
        
        fprintf('  [INFO] Movement: %s\n', move_type);
        
        % Collect spike times for raster and compute PSTH/FR
        spikeTimes_all = {};
        FR_all_trials = [];
        PSTH_all_trials = [];
        
        for row = 1:height(move_tbl)
            spkTimes = (move_tbl.C1{row} - move_tbl.TTL_spk_idx_Start(row)) / AO_spike_fs - 0.05;
            spikeTimes_all{end+1} = spkTimes;
            
            % FR (10 ms bins)
            spikeCounts_FR = histcounts(spkTimes, edges_FR);
            FR_all_trials = [FR_all_trials; spikeCounts_FR/binSize_FR];
            
            % PSTH (50 ms bins)
            spikeCounts_PSTH = histcounts(spkTimes, edges_PSTH);
            PSTH_all_trials = [PSTH_all_trials; spikeCounts_PSTH];
        end
        
        % Compute means and errors
        mean_FR = mean(FR_all_trials, 1);
        mean_PSTH = mean(PSTH_all_trials, 1);
        if useSEM
            err_FR = std(FR_all_trials, 0, 1)/sqrt(size(FR_all_trials,1));
            err_PSTH = std(PSTH_all_trials, 0, 1)/sqrt(size(PSTH_all_trials,1));
        else
            err_FR = std(FR_all_trials, 0, 1);
            err_PSTH = std(PSTH_all_trials, 0, 1);
        end
        
        % ---- Store Summary ----
        summary_Data(end+1,:) = {CaseDate, CaseDate_hem, depth_name, move_type, size(FR_all_trials,1), ...
                                mean(mean_FR), mean(err_FR), mean(mean_PSTH), mean(err_PSTH)};
        
        % ---- Store Full FR/PSTH Data ----
        full_Data(end+1,:) = {CaseDate, CaseDate_hem, depth_name, move_type, ...
                              strjoin(string(time_FR),','), strjoin(string(mean_FR),','), strjoin(string(err_FR),','), ...
                              strjoin(string(time_PSTH),','), strjoin(string(mean_PSTH),','), strjoin(string(err_PSTH),',')};
        
        
        % Plot Raster + PSTH + FR
        figName = sprintf('%s_%s.png', depth_name, matlab.lang.makeValidName(move_type));
        plotRasterPSTH(spikeTimes_all, time_FR, mean_FR, err_FR, time_PSTH, mean_PSTH, err_PSTH, move_type, depth_name, ['±' errLabel]);
        
        % Save Figures
        saveas(gcf, fullfile(figDepthDir, figName));
        savefig(fullfile(figDepthDir, strrep(figName,'.png','.fig')));
        close;
    end
end

fprintf('[INFO] All figures saved in: %s\n', figBaseDir);

% Convert to Tables and Export
summary_Table = cell2table(summary_Data, 'VariableNames', ...
    {'CaseDate','Hemisphere','Depth','MoveType','N_Trials','Mean_FR_Hz',['Err_FR_' errLabel],'Mean_PSTH',['Err_PSTH_' errLabel]});
writetable(summary_Table, fullfile(csvDir, ['FR_PSTH_Summary_' CaseDate '_' errLabel '.csv']));

full_Table = cell2table(full_Data, 'VariableNames', ...
    {'CaseDate','Hemisphere','Depth','MoveType','TimeVector_FR','Mean_FR','Err_FR','TimeVector_PSTH','Mean_PSTH','Err_PSTH'});
writetable(full_Table, fullfile(csvDir, ['FR_PSTH_FullBins_' CaseDate '_' errLabel '.csv']));

fprintf('[INFO] Exported CSVs to: %s\n', csvDir);


end

%% Helper Functions

function depth_name = switchDepth(code)
switch code
    case 't', depth_name = 'dorsal';
    case 'c', depth_name = 'central';
    case 'b', depth_name = 'ventral';
end
end

function plotRasterPSTH(spikeTimes_all, time_FR, mean_FR, err_FR, time_PSTH, mean_PSTH, err_PSTH, move_type, depth_name, errLabel)
figure('Visible','off','Position',[100 100 800 600]);
tiledlayout(3,1,'TileSpacing','compact');

% Raster Plot
nexttile;
hold on;
for t = 1:numel(spikeTimes_all)
    st = spikeTimes_all{t};
    plot(st, t*ones(size(st)), 'k.', 'MarkerSize', 8);
end
xlim([-0.05 0.45]);
ylabel('Trial');
title(sprintf('Raster | %s | %s', depth_name, move_type));

% PSTH Plot
nexttile;
shadedArea(time_PSTH, mean_PSTH, err_PSTH, [0.6 0.6 0.9]);
ylabel(sprintf('Spike Count (50 ms bin %s)', errLabel));
title('PSTH');

% FR Plot
nexttile;
shadedArea(time_FR, mean_FR, err_FR, [0.2 0.6 1]);
xlabel('Time (s)');
ylabel(sprintf('Firing Rate (Hz %s)', errLabel));
title('FR (10 ms bins)');
end

function shadedArea(x, y, err, color)
fill([x fliplr(x)], [y+err fliplr(y-err)], color, 'FaceAlpha',0.3,'EdgeColor','none');
hold on;
plot(x, y, 'Color', color*0.7, 'LineWidth', 2);
end
