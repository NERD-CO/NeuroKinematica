function IO_plotting_vCombined(CaseDate, CaseDate_hem, ephysTbl_Dir, ephys_offset)

% Combined script for original + rescaled raster plots with mean ± SD FR in titles
% Adds p-values and significance stars in subtitles for REST vs Movement plots
% Includes Option A (actual duration) and Option B (fixed window) toggle for FR


%% Directory Setup

curPCname = getenv('COMPUTERNAME');
switch curPCname
    case 'DESKTOP-I5CPDO7'
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    otherwise
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
end

%% Inputs:

CaseDate = '03_23_2023'; % Adjust as needed
% '03_23_2023'
% '05_18_2023_b_bilateral'

ephys_offset = 1;

Case_DataDir = fullfile(IO_DataDir,CaseDate);
ephysTbl_Dir = fullfile(Case_DataDir,'DLC_Ephys');

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


% Example:
% IO_plotting_vCombined('03_23_2023', '', 1);                 % unilateral
% IO_plotting_vCombined('05_18_2023_b_bilateral', 'LSTN', 1); % bilateral


%% CONFIG:

% FR Calculation Logic
useActualDuration = false; % true = Option A, false = Option B
window_FR = [-0.05 0.45];  % Fixed window for Option B

% Raster
AO_spike_fs = 44000;
binSize_FR = 0.01;
% edges_FR = window_FR(1):binSize_FR:window_FR(2);

depth_ids = {'t','c','b'};
depth_labels = {'dorsal STN','central STN','ventral STN'};
depth_colors = [0 0.6 0; 0.6 0 0.8; 0 0.4 0.8];
rest_color = [0.5 0.5 0.5];

sorted_movements = {'HAND OC','HAND PS','ARM EF','REST'};
active_movements = {'HAND OC','HAND PS','ARM EF'};

fontTitle = 14; fontLabel = 12;


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
summary_FR = {}; rest_vs_move_stats = {};

%% Loop through trials

for m = 1:numel(move_types)
    move_type = move_types{m};
    for d = 1:numel(depth_ids)
        depth_code = depth_ids{d};
        depth_name = depth_labels{d};

        move_tbl = All_SpikesPerMove_Tbl(strcmp(All_SpikesPerMove_Tbl.MoveType, move_type) & ...
            contains(All_SpikesPerMove_Tbl.move_trial_ID, depth_code), :);
        if isempty(move_tbl), continue; end

        spike_list = {}; FR_trials = [];
        for row = 1:height(move_tbl)
            spkTimes = (move_tbl.C1{row} - move_tbl.TTL_spk_idx_Start(row)) / AO_spike_fs - 0.05;
            spike_list{end+1} = spkTimes;
            FR_trials(end+1) = computeFR(spkTimes, useActualDuration, window_FR);
        end

        mu_FR = mean(FR_trials); sd_FR = std(FR_trials);
        summary_FR(end+1,:) = {CaseDate, CaseDate_hem, depth_name, move_type, numel(spike_list), round(mu_FR,1), round(sd_FR,1)};

        if strcmp(move_type,'REST')
            restSpikes{d} = spike_list;
        elseif ismember(move_type,active_movements)
            idx = find(strcmp(active_movements,move_type));
            allDepthSpikes{d,idx} = spike_list;
        end
    end
end

%% Plot Original + Rescaled Multi-Movement Raster

multiRasterName = fullfile(figMultiDir, sprintf('MultiMovement_Raster_%s%s.png', CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], '')));
plotMultiMovementRaster(allDepthSpikes, active_movements, depth_labels, depth_colors, multiRasterName, fontTitle, fontLabel, useActualDuration, window_FR);
plotMultiMovementRasterRescaled(allDepthSpikes, active_movements, depth_labels, depth_colors, multiRasterName, fontTitle, fontLabel, useActualDuration, window_FR);

%% REST vs Movement Comparison

if hasREST
    for m = 1:numel(active_movements)
        move_type = active_movements{m};
        figName = fullfile(figMultiDir, sprintf('REST_vs_%s_%s%s.png', matlab.lang.makeValidName(move_type), CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], '')));
        [pVals, cohenD_vals] = plotRestVsMoveRaster(restSpikes, allDepthSpikes(:,m), move_type, depth_labels, depth_colors, rest_color, figName, fontTitle, fontLabel, useActualDuration, window_FR);
        plotRestVsMoveRasterRescaled(restSpikes, allDepthSpikes(:,m), move_type, depth_labels, depth_colors, rest_color, figName, fontTitle, fontLabel, useActualDuration, window_FR);

        for d = 1:3
            rest_vs_move_stats(end+1,:) = {CaseDate, CaseDate_hem, depth_labels{d}, move_type, length(allDepthSpikes{d,m}), length(restSpikes{d}), pVals(d), cohenD_vals(d)};
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



%% === Helper Functions ===

function res = ternary(cond,a,b), if cond,res=a;else,res=b;end,end

%% === Compute FR ===
function FR = computeFR(spikeTimes, useActualDuration, window)
% Compute number of spikes within the window
spikeMask = (spikeTimes >= window(1)) & (spikeTimes <= window(2));
numSpikes = sum(spikeMask);

% Determine correct duration
if useActualDuration
    if any(spikeMask)
        dur = max(spikeTimes(spikeMask)) - min(spikeTimes(spikeMask));
        % Prevent divide-by-zero if only 1 spike
        dur = max(dur, 1e-3);
    else
        dur = window(2) - window(1); % fallback to fixed window
    end
else
    dur = window(2) - window(1); % fixed duration window (e.g., 0.5s)
end

% Compute firing rate
FR = numSpikes / dur;  % in Hz
end

%% === Check Normality ===
function isNormal = checkNormality(data)
% Shapiro-Wilk is not available natively, so we'll use kstest
% If sample size < 4, we’ll assume normal for stability
if numel(data) < 4
    isNormal = true;
    return;
end
data = data(~isnan(data));
data = (data - mean(data)) / std(data); % z-score
[h, ~] = kstest(data); % h=0 means normal
isNormal = (h == 0);
end


%% === Plotting: MultiMovement (Original) ===
function plotMultiMovementRaster(allDepthSpikes, movements, depth_labels, colors, savePath, fontTitle, fontLabel, useActualDuration, window_FR)
figure('Visible','off','Position',[100 100 1300 900]);
tiledlayout(3,numel(movements),'TileSpacing','compact');
for d=1:3
    for m=1:numel(movements)
        nexttile; hold on;
        trials = allDepthSpikes{d,m};
        if isempty(trials)
            text(0.5,0.5,'No Data','HorizontalAlignment','center'); continue;
        end
        FR_vals = cellfun(@(x) computeFR(x,useActualDuration,window_FR),trials);
        mu=mean(FR_vals); sd=std(FR_vals);
        for t2=1:numel(trials)
            plot(trials{t2},t2*ones(size(trials{t2})),'.','Color',colors(d,:),'MarkerSize',8);
        end
        xlim([-0.05 0.45]); ylim([0 numel(trials)+1]);
        if d==3, xlabel('Time (s)','FontSize',fontLabel); end
        if m==1, ylabel(depth_labels{d},'FontSize',fontLabel); end
        title(sprintf('%s (%.1f ± %.1f Hz)',movements{m},mu,sd),'FontSize',fontTitle);
    end
end
sgtitle('Multi-Movement Raster (Original)','FontSize',fontTitle,'FontWeight','bold');
saveas(gcf,savePath); savefig(strrep(savePath,'.png','.fig')); close;
end

% === Plotting: MultiMovement (Rescaled) ===
function plotMultiMovementRasterRescaled(allDepthSpikes,movements,depth_labels,colors,savePath,fontTitle,fontLabel,useActualDuration,window_FR)
figure('Visible','off','Position',[100 100 1300 900]);
tiledlayout(3,numel(movements),'TileSpacing','compact');
scaleFactor=0.5;
for d=1:3
    for m=1:numel(movements)
        nexttile; hold on;
        trials = allDepthSpikes{d,m};
        if isempty(trials)
            text(0.5,0.5,'No Data','HorizontalAlignment','center'); continue;
        end
        FR_vals = cellfun(@(x) computeFR(x,useActualDuration,window_FR),trials);
        mu=mean(FR_vals); sd=std(FR_vals);
        for t2=1:numel(trials)
            yScaled=t2*scaleFactor;
            plot(trials{t2},yScaled*ones(size(trials{t2})),'.','Color',colors(d,:),'MarkerSize',8);
        end
        xlim([-0.05 0.45]); ylim([0 max(1,ceil(numel(trials)*scaleFactor))]);
        if d==3,xlabel('Time (s)','FontSize',fontLabel);end
        if m==1,ylabel(depth_labels{d},'FontSize',fontLabel);end
        title(sprintf('%s (%.1f ± %.1f Hz)',movements{m},mu,sd),'FontSize',fontTitle);
    end
end
sgtitle('Multi-Movement Raster (Rescaled)','FontSize',fontTitle,'FontWeight','bold');
saveas(gcf,strrep(savePath,'.png','_rescaled.png')); savefig(strrep(savePath,'.png','_rescaled.fig')); close;
end


%% === Plotting: REST vs Movement (Original) ===
function [pVals,cohenD_vals]=plotRestVsMoveRaster(restSpikes,moveSpikes,move_type,depth_labels,colors,restColor,savePath,fontTitle,fontLabel,useActualDuration,window_FR)
figure('Visible','off','Position',[100 100 900 800]);
tiledlayout(3,2,'TileSpacing','compact');
pVals=nan(1,3); cohenD_vals=nan(1,3);
for d=1:3
    % REST
    ax1=nexttile; hold on;
    trialsR=restSpikes{d};
    FR_rest=cellfun(@(x)computeFR(x,useActualDuration,window_FR),trialsR);
    muR=mean(FR_rest); sdR=std(FR_rest);
    for t=1:numel(trialsR)
        plot(trialsR{t},t*ones(size(trialsR{t})),'.','Color',restColor,'MarkerSize',8);
    end
    title(ax1,sprintf('%s | REST (%.1f ± %.1f)',depth_labels{d},muR,sdR),'FontSize',fontTitle);
    xlim([-0.05 0.45]); ylim([0 numel(trialsR)+1]);

    % Movement
    ax2=nexttile; hold on;
    trialsM=moveSpikes{d};
    FR_move=cellfun(@(x)computeFR(x,useActualDuration,window_FR),trialsM);
    muM=mean(FR_move); sdM=std(FR_move);
    for t=1:numel(trialsM)
        plot(trialsM{t},t*ones(size(trialsM{t})),'.','Color',colors(d,:),'MarkerSize',8);
    end
    title(ax2,sprintf('%s | %s (%.1f ± %.1f)',depth_labels{d},move_type,muM,sdM),'FontSize',fontTitle);

    if ~isempty(trialsM)&&~isempty(trialsR)
        % Check for normality
        isNormal1 = checkNormality(FR_rest);
        isNormal2 = checkNormality(FR_move);

        if isNormal1 && isNormal2
            [~, p] = ttest2(FR_move, FR_rest);
            testType = 't-test';
        else
            p = ranksum(FR_move, FR_rest);
            testType = 'ranksum';
        end

        pVals(d) = p;
        cohenD_vals(d) = computeCohensD(FR_move, FR_rest);  % still valid to show
        stars = getSigStars(p);

        subtitle(ax2, sprintf('p=%.3f (%s, %s)\nd=%.2f', p, stars, testType, cohenD_vals(d)));
    end
    xlim([-0.05 0.45]); ylim([0 numel(trialsM)+1]);
end
sgtitle(sprintf('REST vs %s (Original)',move_type),'FontSize',fontTitle,'FontWeight','bold');
saveas(gcf,savePath); savefig(strrep(savePath,'.png','.fig')); close;
end


% === Plotting: REST vs Movement (Rescaled) ===
function plotRestVsMoveRasterRescaled(restSpikes,moveSpikes,move_type,depth_labels,colors,restColor,savePath,fontTitle,fontLabel,useActualDuration,window_FR)
figure('Visible','off','Position',[100 100 900 800]);
tiledlayout(3,2,'TileSpacing','compact');
scaleFactor=0.5;
for d=1:3
    % REST
    ax1=nexttile; hold on;
    trialsR=restSpikes{d};
    FR_rest=cellfun(@(x)computeFR(x,useActualDuration,window_FR),trialsR);
    muR=mean(FR_rest); sdR=std(FR_rest);
    for t=1:numel(trialsR)
        yScaled=t*scaleFactor;
        plot(trialsR{t},yScaled*ones(size(trialsR{t})),'.','Color',restColor,'MarkerSize',8);
    end
    title(ax1,sprintf('%s | REST (%.1f ± %.1f)',depth_labels{d},muR,sdR),'FontSize',fontTitle);
    ylim([0 max(1,ceil(numel(trialsR)*scaleFactor))]); xlim([-0.05 0.45]);

    % Movement
    ax2=nexttile; hold on;
    trialsM=moveSpikes{d};
    FR_move=cellfun(@(x)computeFR(x,useActualDuration,window_FR),trialsM);
    muM=mean(FR_move); sdM=std(FR_move);
    for t=1:numel(trialsM)
        yScaled=t*scaleFactor;
        plot(trialsM{t},yScaled*ones(size(trialsM{t})),'.','Color',colors(d,:),'MarkerSize',8);
    end
    title(ax2,sprintf('%s | %s (%.1f ± %.1f)',depth_labels{d},move_type,muM,sdM),'FontSize',fontTitle);

    if ~isempty(trialsM)&&~isempty(trialsR)
        % Check for normality
        isNormal1 = checkNormality(FR_rest);
        isNormal2 = checkNormality(FR_move);

        if isNormal1 && isNormal2
            [~, p] = ttest2(FR_move, FR_rest);
            testType = 't-test';
        else
            p = ranksum(FR_move, FR_rest);
            testType = 'ranksum';
        end

        pVals(d) = p;
        cohenD_vals(d) = computeCohensD(FR_move, FR_rest);  % still valid to show
        stars = getSigStars(p);

        subtitle(ax2, sprintf('p=%.3f (%s, %s)\nd=%.2f', p, stars, testType, cohenD_vals(d)));
    end
    ylim([0 max(1,ceil(numel(trialsM)*scaleFactor))]); xlim([-0.05 0.45]);

end
sgtitle(sprintf('REST vs %s (Rescaled)',move_type),'FontSize',fontTitle,'FontWeight','bold');
saveas(gcf,strrep(savePath,'.png','_rescaled.png')); savefig(strrep(savePath,'.png','_rescaled.fig')); close;
end


%% Stats Helpers
function stars=getSigStars(p)
if p<0.001,stars='***';elseif p<0.01,stars='**';elseif p<0.05,stars='*';else,stars='ns';end
end

function d=computeCohensD(group1,group2)
n1=length(group1); n2=length(group2);
s1=var(group1); s2=var(group2);
pooledSD=sqrt(((n1-1)*s1+(n2-1)*s2)/(n1+n2-2));
d=(mean(group1)-mean(group2))/pooledSD;
end
