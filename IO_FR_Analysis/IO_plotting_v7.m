function IO_plotting_v7(CaseDate, CaseDate_hem, ephys_offset)

% IO_plotting_v7:
% Generates:
%   - Multi-Movement Raster (Hand OC, Hand PS, Arm EF)
%   - REST vs Movement Rasters (with FR stats, significance stars)
%   - CSV exports: FR summary, full FR bins, REST vs Move stats (incl. p-value + Cohen's d)
%   - Aggregated MoveN Summary for each STN depth

clear; clc;
addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis'

%% Inputs:

CaseDate = '03_23_2023'; %% adjust as needed
% '03_23_2023'; % unilateral example
% '05_18_2023_b_bilateral'; bilateral example

ephys_offset = 1;

% Example:
% IO_plotting_v7('03_23_2023', '', 1);
% IO_plotting_v7('05_18_2023_b_bilateral', 'LSTN', 1);


%% Config - script parameters

AO_spike_fs = 44000;

binSize_FR = 0.01;
window_FR = [-0.05 0.45];
edges_FR = window_FR(1):binSize_FR:window_FR(2);
time_FR = edges_FR(1:end-1)+binSize_FR/2;

plotPSTH = false; % Toggle PSTH/FR bin CSV outputs
binSize_PSTH = 0.05;
edges_PSTH = window_FR(1):binSize_PSTH:window_FR(2);
time_PSTH = edges_PSTH(1:end-1)+binSize_PSTH/2;

depth_ids = {'t','c','b'};
depth_labels = {'dorsal STN','central STN','ventral STN'};
depth_colors = [0 0.6 0; 0.6 0 0.8; 0 0.4 0.8];
rest_color = [0.5 0.5 0.5];
sorted_movements = {'HAND OC','HAND PS','ARM EF','REST'};
active_movements = {'HAND OC','HAND PS','ARM EF'};

fontTitle = 14; fontLabel = 12;
useSEM = false;  % true=SEM, false=STD
errLabel = ternary(useSEM,'SEM','STD');


%% Directory Setup

curPCname = getenv('COMPUTERNAME');
switch curPCname
    case 'DESKTOP-I5CPDO7'
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    otherwise
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
end

Case_DataDir = fullfile(IO_DataDir,CaseDate);
ephysTbl_Dir = fullfile(Case_DataDir,'DLC_Ephys');

%% Handle bilateral cases

isBilateral = contains(CaseDate,'bilateral','IgnoreCase',true);
if isBilateral
    fprintf('[INFO] Bilateral case detected: %s\n',CaseDate);
    CaseDate_hem = input('Enter hemisphere (LSTN or RSTN): ','s');
    if ~ismember(CaseDate_hem,{'LSTN','RSTN'}), error('Invalid hemisphere'); end
    ephysTbl_Dir = fullfile(ephysTbl_Dir,CaseDate_hem);
else 
    CaseDate_hem = '';
end
fprintf('[INFO] Loading spike data from: %s\n', ephysTbl_Dir);

%% Output Directories

figDir = fullfile(ephysTbl_Dir,'Figures'); if ~exist(figDir,'dir'), mkdir(figDir); end
figMultiDir = fullfile(figDir,'Raster_MultiMovement_and_RestComparison'); if ~exist(figMultiDir,'dir'), mkdir(figMultiDir); end
csvDir = fullfile(ephysTbl_Dir,'CSV_Outputs'); if ~exist(csvDir,'dir'), mkdir(csvDir); end

%% Load Spike Table

cd(ephysTbl_Dir);
Tbl_list = dir('*Spikes*.mat'); Tbl_names = {Tbl_list.name};
if ephys_offset
    spk_case = Tbl_names{contains(Tbl_names,'Spikes') & contains(Tbl_names,'offset')};
else
    spk_case = Tbl_names{contains(Tbl_names,'Spikes') & ~contains(Tbl_names,'offset')};
end
if isempty(spk_case), error('No spike table found'); end
load(spk_case,'All_SpikesPerMove_Tbl');

%% Identify move types

move_types = intersect(sorted_movements, unique(All_SpikesPerMove_Tbl.MoveType),'stable');
hasREST = any(strcmp(move_types,'REST'));

%% Initialize storage

allDepthSpikes = cell(3, numel(active_movements));
restSpikes = cell(3,1);
summary_FR = {}; full_FR_bins = {}; rest_vs_move_stats = {};

%% Loop through original table (no aggregation of frontCam, sideCam redundancy)

for m = 1:numel(move_types)
    move_type = move_types{m};
    for d = 1:numel(depth_ids)
        depth_code = depth_ids{d};
        depth_name = depth_labels{d};
        
        move_tbl = All_SpikesPerMove_Tbl(strcmp(All_SpikesPerMove_Tbl.MoveType, move_type) & ...
                                         contains(All_SpikesPerMove_Tbl.move_trial_ID, depth_code), :);
        if isempty(move_tbl), continue; end
        
        % Collect spike times
        spike_list = {};
        FR_trials = [];
        for row = 1:height(move_tbl)
            spkTimes = (move_tbl.C1{row} - move_tbl.TTL_spk_idx_Start(row)) / AO_spike_fs - 0.05;
            spike_list{end+1} = spkTimes;
            FR_trials = [FR_trials; histcounts(spkTimes, edges_FR)/binSize_FR];
        end
        
        % Compute FR stats
        mean_FR_val = mean(FR_trials,'all'); 
        std_FR_val = std(mean(FR_trials,2));
        summary_FR(end+1,:) = {CaseDate, CaseDate_hem, depth_name, move_type, size(FR_trials,1), round(mean_FR_val,1), round(std_FR_val,1)};
        
        % Store for rasters
        if strcmp(move_type,'REST')
            restSpikes{d} = spike_list;
        elseif ismember(move_type,active_movements)
            idx = find(strcmp(active_movements,move_type));
            allDepthSpikes{d,idx} = spike_list;
        end
    end
end

%% Plot Multi-Movement Raster

multiRasterName = fullfile(figMultiDir, sprintf('MultiMovement_Raster_%s%s.png', CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], '')));
plotMultiMovementRaster(allDepthSpikes, active_movements, depth_labels, depth_colors, multiRasterName, fontTitle, fontLabel);

%% REST vs Movement Comparison

if hasREST
    for m = 1:numel(active_movements)
        move_type = active_movements{m};
        figName = fullfile(figMultiDir, sprintf('REST_vs_%s_%s%s.png', matlab.lang.makeValidName(move_type), CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], '')));
        [pVals, cohenD_vals] = plotRestVsMoveRaster(restSpikes, allDepthSpikes(:,m), move_type, depth_labels, depth_colors, rest_color, figName, fontTitle, fontLabel, edges_FR, binSize_FR);
        
        % Add stats for CSV
        for d = 1:3
            rest_vs_move_stats(end+1,:) = {CaseDate, CaseDate_hem, depth_labels{d}, move_type, ...
                length(allDepthSpikes{d,m}), length(restSpikes{d}), pVals(d), cohenD_vals(d)};
        end
    end
end

%% Export CSVs

summaryTable = cell2table(summary_FR, 'VariableNames', {'CaseDate','Hemisphere','Depth','MoveType','N_Trials','Mean_FR_Hz','SD_FR_Hz'});
writetable(summaryTable, fullfile(csvDir, sprintf('FR_Summary_%s%s.csv', CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], ''))));

if ~isempty(rest_vs_move_stats)
    statsTable = cell2table(rest_vs_move_stats, 'VariableNames', {'CaseDate','Hemisphere','Depth','MoveType','N_MoveTrials','N_RESTTrials','p_value','Cohens_d'});
    writetable(statsTable, fullfile(csvDir, sprintf('RESTvsMove_Stats_%s%s.csv', CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], ''))));
end

fprintf('[INFO] CSV export complete: %s\n', csvDir);
end

%% ---------------------- Helper Functions ---------------------- %%

function res = ternary(cond,a,b), if cond,res=a;else,res=b;end,end

function plotMultiMovementRaster(allDepthSpikes, movements, depth_labels, colors, savePath, fontTitle, fontLabel)
figure('Visible','off','Position',[100 100 1300 900]);
tiledlayout(3, numel(movements), 'TileSpacing','compact');
for d = 1:3
    for m = 1:numel(movements)
        nexttile; hold on;
        trials = allDepthSpikes{d, m};
        if isempty(trials)
            text(0.2,0.5,'No Data','HorizontalAlignment','center','Color',[0.7 0.7 0.7]);
            continue; 
        end
        for t2 = 1:numel(trials)
            plot(trials{t2}, t2*ones(size(trials{t2})), '.', 'Color', colors(d,:), 'MarkerSize', 8);
        end
        xlim([-0.05 0.45]); ylim([0 numel(trials)+1]);
        if d == 3, xlabel('Time (s)', 'FontSize', fontLabel); end
        if m == 1, ylabel(depth_labels{d}, 'FontSize', fontLabel); end
        title(sprintf('%s (n=%d)', movements{m}, numel(trials)), 'FontSize', fontTitle);
    end
end
sgtitle('Multi-Movement Raster', 'FontSize', fontTitle, 'FontWeight', 'bold');
saveas(gcf, savePath); savefig(strrep(savePath, '.png', '.fig')); close;
end

function [pVals, cohenD_vals] = plotRestVsMoveRaster(restSpikes, moveSpikes, move_type, depth_labels, colors, restColor, savePath, fontTitle, fontLabel, edges_FR, binSize_FR)
figure('Visible','off','Position',[100 100 900 800]);
tiledlayout(3,2,'TileSpacing','compact');
pVals = nan(1,3); cohenD_vals = nan(1,3);
for d=1:3
    % REST raster
    nexttile; hold on;
    trialsR = restSpikes{d};
    for t=1:numel(trialsR)
        plot(trialsR{t}, t*ones(size(trialsR{t})), '.', 'Color', [0.5 0.5 0.5], 'MarkerSize',8);
    end
    xlim([-0.05 0.45]); ylim([0 numel(trialsR)+1]);
    if d==3, xlabel('Time (s)','FontSize',fontLabel); end
    ylabel('Rep','FontSize',fontLabel);
    title(sprintf('%s | REST', depth_labels{d}), 'FontSize', fontTitle);
    
    % Movement raster
    nexttile; hold on;
    trialsM = moveSpikes{d};
    for t=1:numel(trialsM)
        plot(trialsM{t}, t*ones(size(trialsM{t})), '.', 'Color', colors(d,:), 'MarkerSize',8);
    end
    xlim([-0.05 0.45]); ylim([0 numel(trialsM)+1]);
    if d==3, xlabel('Time (s)','FontSize',fontLabel); end
    title(sprintf('%s | %s', depth_labels{d}, move_type), 'FontSize', fontTitle);
    
    % Compute FR stats
    if ~isempty(trialsM) && ~isempty(trialsR)
        FR_move = cellfun(@(x) numel(x)/(binSize_FR*(length(edges_FR)-1)), trialsM);
        FR_rest = cellfun(@(x) numel(x)/(binSize_FR*(length(edges_FR)-1)), trialsR);
        mu = mean(FR_move); sd = std(FR_move);
        text(0.15, numel(trialsM)+1, sprintf('\\mu = %.1f \\pm %.1f Hz', mu, sd), 'FontSize', 10, 'Color', colors(d,:));
        
        [~,p] = ttest2(FR_move, FR_rest);
        pVals(d) = p;
        cohenD_vals(d) = computeCohensD(FR_move, FR_rest);
        
        stars = getSigStars(p);
        text(0.3, numel(trialsM)+1, stars, 'FontSize', 14, 'FontWeight', 'bold', 'Color', colors(d,:));
    end
end
sgtitle(sprintf('REST vs %s', move_type), 'FontSize', fontTitle,'FontWeight','bold');
saveas(gcf, savePath); savefig(strrep(savePath,'.png','.fig')); close;
end

function stars = getSigStars(p)
if p < 0.001, stars = '***';
elseif p < 0.01, stars = '**';
elseif p < 0.05, stars = '*';
else, stars = 'ns'; end
end

function d = computeCohensD(group1, group2)
n1 = length(group1); n2 = length(group2);
s1 = var(group1); s2 = var(group2);
pooledSD = sqrt(((n1-1)*s1 + (n2-1)*s2)/(n1+n2-2));
d = (mean(group1) - mean(group2)) / pooledSD;
end