function IO_plotting_v5(CaseDate, CaseDate_hem, ephys_offset)

% IO_plotting_v5:
% Generates:
%  - Raster plots per movement type & STN depth (stacked & polished)
%  - PSTH & FR plots per movement type
%  - Exports summary and full data as CSV

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

Case_DataDir = fullfile(IO_DataDir, CaseDate);
ephysTbl_Dir = fullfile(Case_DataDir, 'DLC_Ephys');

% Handle bilateral
isBilateral = contains(CaseDate, 'bilateral', 'IgnoreCase', true);
if isBilateral
    fprintf('[INFO] Bilateral case: %s\n', CaseDate);
    CaseDate_hem = input('Enter hemisphere (LSTN or RSTN): ', 's');
    if ~ismember(CaseDate_hem, {'LSTN','RSTN'})
        error('Invalid hemisphere input.');
    end
    ephysTbl_Dir = fullfile(ephysTbl_Dir, CaseDate_hem);
end
fprintf('[INFO] Data directory: %s\n', ephysTbl_Dir);


%% CONFIG - script parameters, AO_fs, bin sizes

useSEM = true; % true = SEM, false = STD
AO_spike_fs = 44000;

% Bin settings
binSize_FR = 0.01;  % 10 ms
window_FR  = [-0.05 0.45];
edges_FR   = window_FR(1):binSize_FR:window_FR(2);
time_FR    = edges_FR(1:end-1)+binSize_FR/2;

binSize_PSTH = 0.05; % 50 ms
edges_PSTH   = window_FR(1):binSize_PSTH:window_FR(2);
time_PSTH    = edges_PSTH(1:end-1)+binSize_PSTH/2;

depth_ids = {'t','c','b'}; % dorsal, central, ventral
depth_labels = {'dorsal STN','central STN','ventral STN'};
colors = [0 0.6 0; 0.6 0 0.8; 0 0.4 0.8]; % green, purple, blue

errLabel = ternary(useSEM, 'SEM', 'STD');

%% Load Spikes Table

cd(ephysTbl_Dir);
Tbl_list = dir('*Spikes*.mat');
Tbl_names = {Tbl_list.name};

if ephys_offset
    spk_case = Tbl_names{contains(Tbl_names,'Spikes') & contains(Tbl_names,'offset')};
else
    spk_case = Tbl_names{contains(Tbl_names,'Spikes') & ~contains(Tbl_names,'offset')};
end
if isempty(spk_case), error('No spike table found.'); end
load(spk_case,'All_SpikesPerMove_Tbl');
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


%% Loop through movement types

% Initialize CSV storage
summary_Data = {};
full_Data = {};

for m = 1:numel(move_types)
    move_type = move_types{m};
    fprintf('\n[INFO] Processing Movement: %s\n', move_type);
    
    % Collect for stacked raster
    spikeTimes_byDepth = {{} {} {}}; 
    
    for d = 1:numel(depth_ids)
        depth_code = depth_ids{d};
        depth_name = switchDepth(depth_code);
        depth_tbl = All_SpikesPerMove_Tbl(contains(All_SpikesPerMove_Tbl.move_trial_ID, depth_code),:);
        move_tbl = depth_tbl(strcmp(depth_tbl.MoveType,move_type),:);
        
        if isempty(move_tbl), continue; end
        
        fprintf('  Depth: %s (%d trials)\n', depth_name, height(move_tbl));
        
        spikeTimes_list = {};
        FR_all_trials = [];
        PSTH_all_trials = [];
        
        for row = 1:height(move_tbl)
            spkTimes = (move_tbl.C1{row}-move_tbl.TTL_spk_idx_Start(row))/AO_spike_fs - 0.05;
            spikeTimes_list{end+1} = spkTimes;
            
            spikeCounts_FR = histcounts(spkTimes,edges_FR);
            FR_all_trials = [FR_all_trials; spikeCounts_FR/binSize_FR];
            
            spikeCounts_PSTH = histcounts(spkTimes,edges_PSTH);
            PSTH_all_trials = [PSTH_all_trials; spikeCounts_PSTH];
        end
        
        spikeTimes_byDepth{d} = spikeTimes_list;
        
        % Compute means & error
        mean_FR = mean(FR_all_trials,1);
        mean_PSTH = mean(PSTH_all_trials,1);
        if useSEM
            err_FR = std(FR_all_trials,0,1)/sqrt(size(FR_all_trials,1));
            err_PSTH = std(PSTH_all_trials,0,1)/sqrt(size(PSTH_all_trials,1));
        else
            err_FR = std(FR_all_trials,0,1);
            err_PSTH = std(PSTH_all_trials,0,1);
        end
        
        % Store summary for CSV
        summary_Data(end+1,:) = {CaseDate, CaseDate_hem, depth_name, move_type, size(FR_all_trials,1), ...
                                mean(mean_FR), mean(err_FR), mean(mean_PSTH), mean(err_PSTH)};
        
        % Full bin data
        full_Data(end+1,:) = {CaseDate, CaseDate_hem, depth_name, move_type, ...
                              strjoin(string(time_FR),','), strjoin(string(mean_FR),','), strjoin(string(err_FR),','), ...
                              strjoin(string(time_PSTH),','), strjoin(string(mean_PSTH),','), strjoin(string(err_PSTH),',')};
        
        % Plot individual PSTH + FR
        figName = sprintf('%s_%s.png', depth_name, matlab.lang.makeValidName(move_type));
        plotRasterPSTH(spikeTimes_list, time_FR, mean_FR, err_FR, time_PSTH, mean_PSTH, err_PSTH, move_type, depth_name, ['±' errLabel]);
        saveas(gcf, fullfile(figBaseDir, figName));
        savefig(fullfile(figBaseDir, strrep(figName,'.png','.fig')));
        close;
    end
    
    % Plot stacked raster for this movement
    stackedName = fullfile(figBaseDir, sprintf('StackedRaster_%s.png', matlab.lang.makeValidName(move_type)));
    plotStackedRaster(spikeTimes_byDepth, move_type, depth_labels, colors, stackedName);
end

%% Export CSVs
summary_Table = cell2table(summary_Data,'VariableNames',...
    {'CaseDate','Hemisphere','Depth','MoveType','N_Trials','Mean_FR_Hz',['Err_FR_' errLabel],'Mean_PSTH',['Err_PSTH_' errLabel]});
writetable(summary_Table, fullfile(csvDir,['FR_PSTH_Summary_' CaseDate '_' errLabel '.csv']));

full_Table = cell2table(full_Data,'VariableNames',...
    {'CaseDate','Hemisphere','Depth','MoveType','TimeVector_FR','Mean_FR','Err_FR','TimeVector_PSTH','Mean_PSTH','Err_PSTH'});
writetable(full_Table, fullfile(csvDir,['FR_PSTH_FullBins_' CaseDate '_' errLabel '.csv']));

fprintf('[INFO] CSVs exported: %s\n', csvDir);

end

%% Helper Functions
function depth_name = switchDepth(code)
switch code
    case 't', depth_name = 'dorsal';
    case 'c', depth_name = 'central';
    case 'b', depth_name = 'ventral';
end
end

function res = ternary(cond,a,b)
if cond, res=a; else, res=b; end
end

function plotRasterPSTH(spikeTimes_all, time_FR, mean_FR, err_FR, time_PSTH, mean_PSTH, err_PSTH, move_type, depth_name, errLabel)
figure('Visible','off','Position',[100 100 800 600]);
tiledlayout(3,1,'TileSpacing','compact');

nexttile; hold on;
for t = 1:numel(spikeTimes_all)
    st = spikeTimes_all{t};
    plot(st, t*ones(size(st)), 'k.', 'MarkerSize', 8);
end
xlim([-0.05 0.45]);
ylabel('Trial');
title(sprintf('Raster | %s | %s', depth_name, move_type));

nexttile;
shadedArea(time_PSTH, mean_PSTH, err_PSTH, [0.6 0.6 0.9]);
ylabel(sprintf('Spike Count (50ms %s)', errLabel));
title('PSTH');

nexttile;
shadedArea(time_FR, mean_FR, err_FR, [0.2 0.6 1]);
xlabel('Time (s)');
ylabel(sprintf('Firing Rate (Hz %s)', errLabel));
title('FR (10 ms bins)');
end

function shadedArea(x,y,err,color)
fill([x fliplr(x)],[y+err fliplr(y-err)],color,'FaceAlpha',0.3,'EdgeColor','none');
hold on; plot(x,y,'Color',color*0.7,'LineWidth',2);
end

function plotStackedRaster(spikeTimes_byDepth, move_type, depth_labels, colors, savePath)
figure('Visible','off','Position',[100 100 900 700]);
tiledlayout(3,1,'TileSpacing','compact');

% Ensure dorsal-central-ventral order
panel_order = [1 2 3]; % dorsal, central, ventral
for p = 1:3
    d = panel_order(p);
    trials = spikeTimes_byDepth{d};
    nexttile; hold on;
    
    for t = 1:numel(trials)
        st = trials{t};
        plot(st, t*ones(size(st)), '.', 'Color', colors(d,:), 'MarkerSize', 8);
    end
    
    xlim([-0.1 0.4]);
    ylim([0 numel(trials)+1]);
    if p < 3
        set(gca,'XTickLabel',[]); % Hide X labels for top/mid panels
    else
        xlabel('Time (s)');
    end
    
    ylabel('Repetition');
    title(depth_labels{d},'FontWeight','bold');
end

sgtitle(sprintf('STN Spike Times | %s', move_type),'FontSize',14,'FontWeight','bold');

% Legend (depth color map) in a separate invisible axis
lgd = legend(depth_labels,'Orientation','horizontal','Box','off');
lgd.Layout.Tile = 'south';

saveas(gcf, savePath);
savefig(strrep(savePath,'.png','.fig'));
close;
end
