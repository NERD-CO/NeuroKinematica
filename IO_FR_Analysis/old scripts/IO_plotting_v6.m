function IO_plotting_v6(CaseDate, CaseDate_hem, ephys_offset)

% IO_plotting_v6:
% Generates:
%   - Raster, PSTH, FR per movement/depth
%   - Stacked rasters per movement
%   - Multi-movement stacked raster
%   - REST vs Movement comparison (raster-only + raster+PSTH/FR)
%   - CSV summaries: FR/PSTH and baseline comparisons
%
% Example:
% IO_plotting_v6('05_18_2023_b_bilateral','LSTN',1);

clear; clc;
addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis'

%% Inputs: 

CaseDate = '03_23_2023'; %% adjust as needed
% '03_23_2023'; % unilateral example
% '05_18_2023_b_bilateral'; bilateral example


ephys_offset = 1;

% Example:
% IO_plotting_v4('03_23_2023', '', 1);
% IO_plotting_v4('05_18_2023_b_bilateral', 'LSTN', 1);


%% Config

smoothPlots = true;      % Gaussian smoothing toggle
useSEM = true;           % true=SEM, false=STD
AO_spike_fs = 44000;
binSize_FR = 0.01; window_FR = [-0.05 0.45];
edges_FR = window_FR(1):binSize_FR:window_FR(2); time_FR = edges_FR(1:end-1)+binSize_FR/2;
binSize_PSTH = 0.05;
edges_PSTH = window_FR(1):binSize_PSTH:window_FR(2); time_PSTH = edges_PSTH(1:end-1)+binSize_PSTH/2;
depth_ids = {'t','c','b'}; depth_labels = {'dorsal STN','central STN','ventral STN'};
depth_colors = [0 0.6 0; 0.6 0 0.8; 0 0.4 0.8]; rest_color = [0.5 0.5 0.5];
sorted_movements = {'HAND OC','HAND PS','ARM EF','REST'};
errLabel = ternary(useSEM,'SEM','STD'); fontTitle = 14; fontLabel = 12;

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
end
fprintf('[INFO] Loading spike data from: %s\n',ephysTbl_Dir);

if ~isBilateral
    CaseDate_hem = '';  % Ensure variable exists for unilateral cases
end


%% Output dirs

figBaseDir = fullfile(ephysTbl_Dir,'Figures');
figRasterDir = fullfile(figBaseDir,'Raster_PSTH');
figMultiDir  = fullfile(figBaseDir,'Raster_MultiMovement_and_RestComparison');
if ~exist(figRasterDir,'dir'), mkdir(figRasterDir); end
if ~exist(figMultiDir,'dir'), mkdir(figMultiDir); end
csvDir = fullfile(ephysTbl_Dir,'CSV_Outputs'); if ~exist(csvDir,'dir'), mkdir(csvDir); end

%% Load Spikes Table
cd(ephysTbl_Dir);
Tbl_list = dir('*Spikes*.mat'); Tbl_names = {Tbl_list.name};
if ephys_offset
    spk_case = Tbl_names{contains(Tbl_names,'Spikes') & contains(Tbl_names,'offset')};
else
    spk_case = Tbl_names{contains(Tbl_names,'Spikes') & ~contains(Tbl_names,'offset')};
end
if isempty(spk_case), error('No spike table found'); end
load(spk_case,'All_SpikesPerMove_Tbl');

move_types = intersect(sorted_movements,unique(All_SpikesPerMove_Tbl.MoveType),'stable');
hasREST = any(strcmp(move_types,'REST'));

%% CSV storage
summary_Data = {}; full_Data = {}; baseline_summary = {};

%% ============ MAIN LOOP ============
allDepthSpikes = cell(3, numel(move_types)-1); % For multi-movement raster
for m = 1:numel(move_types)
    move_type = move_types{m};
    fprintf('\n[INFO] Processing Movement: %s\n',move_type);
    spikeTimes_byDepth = {{} {} {}}; % For stacked raster
    for d = 1:numel(depth_ids)
        depth_code = depth_ids{d}; depth_name = depth_labels{d};
        depth_tbl = All_SpikesPerMove_Tbl(contains(All_SpikesPerMove_Tbl.move_trial_ID,depth_code),:);
        move_tbl = depth_tbl(strcmp(depth_tbl.MoveType,move_type),:);
        if isempty(move_tbl), continue; end
        
        spike_list = {}; FR_trials=[]; PSTH_trials=[];
        for row=1:height(move_tbl)
            spkTimes=(move_tbl.C1{row}-move_tbl.TTL_spk_idx_Start(row))/AO_spike_fs-0.05;
            spike_list{end+1}=spkTimes;
            FR_trials=[FR_trials;histcounts(spkTimes,edges_FR)/binSize_FR];
            PSTH_trials=[PSTH_trials;histcounts(spkTimes,edges_PSTH)];
        end
        spikeTimes_byDepth{d}=spike_list;
        if m<=3, allDepthSpikes{d,m}=spike_list; end
        
        % Compute curves
        mean_FR=mean(FR_trials,1); mean_PSTH=mean(PSTH_trials,1);
        if useSEM, err_FR=std(FR_trials,0,1)/sqrt(size(FR_trials,1)); err_PSTH=std(PSTH_trials,0,1)/sqrt(size(PSTH_trials,1));
        else, err_FR=std(FR_trials,0,1); err_PSTH=std(PSTH_trials,0,1); end
        if smoothPlots, mean_FR=smoothdata(mean_FR,'gaussian',5); mean_PSTH=smoothdata(mean_PSTH,'gaussian',3); end
        
        summary_Data(end+1,:)={CaseDate,CaseDate_hem,depth_name,move_type,size(FR_trials,1),mean(mean_FR),mean(err_FR),mean(mean_PSTH),mean(err_PSTH)};
        full_Data(end+1,:)={CaseDate,CaseDate_hem,depth_name,move_type,strjoin(string(time_FR),','),strjoin(string(mean_FR),','),strjoin(string(err_FR),','),strjoin(string(time_PSTH),','),strjoin(string(mean_PSTH),','),strjoin(string(err_PSTH),',')};
        
        % Save individual raster+PSTH+FR
        figName = sprintf('%s%s_%s.png', CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], ''), matlab.lang.makeValidName(move_type));
        plotRasterPSTH(spike_list,time_FR,mean_FR,err_FR,time_PSTH,mean_PSTH,err_PSTH,move_type,depth_name,['±' errLabel],fontTitle,fontLabel);
        saveas(gcf,fullfile(figRasterDir,figName)); 
        savefig(fullfile(figRasterDir,strrep(figName,'.png','.fig'))); close;
    end
    
    % Save stacked raster for this movement
    figStackName = fullfile(figMultiDir, sprintf('StackedRaster_%s%s_%s.png', CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], ''), matlab.lang.makeValidName(move_type)));
    plotStackedRaster(spikeTimes_byDepth,move_type,depth_labels,depth_colors,figStackName,fontTitle,fontLabel);
end

%% Prepare Multi-Movement Raster Data 

if hasREST
    restSpikes = cell(3,1);
    for d = 1:3
        depth_code = depth_ids{d};
        rest_tbl_depth = All_SpikesPerMove_Tbl(strcmp(All_SpikesPerMove_Tbl.MoveType,'REST') & ...
                                  contains(All_SpikesPerMove_Tbl.move_trial_ID,depth_code),:);
        trials_list = {};
        for r = 1:height(rest_tbl_depth)
            spkTimes = (rest_tbl_depth.C1{r} - rest_tbl_depth.TTL_spk_idx_Start(r))/AO_spike_fs - 0.05;
            trials_list{end+1} = spkTimes;
        end
        restSpikes{d,1} = trials_list;
    end
    % Insert REST as first column
    allDepthSpikes = [restSpikes, allDepthSpikes];
    multiMovements = [{'REST'}, sorted_movements(1:3)];
else
    multiMovements = sorted_movements(1:3);
end

%% Plot Multi-Movement Raster

multiRasterName = fullfile(figMultiDir, sprintf('MultiMovement_Raster_%s%s.png', CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], '')));
plotMultiMovementRaster(allDepthSpikes, multiMovements, depth_labels, depth_colors, rest_color, multiRasterName, fontTitle, fontLabel);

%% REST vs Movement comparisons

if hasREST
    fprintf('[INFO] Generating REST vs Movement comparisons...\n');
    rest_tbl=All_SpikesPerMove_Tbl(strcmp(All_SpikesPerMove_Tbl.MoveType,'REST'),:);
    for m=1:3
        move_type=sorted_movements{m};
        for d=1:3
            depth_code=depth_ids{d}; depth_name=depth_labels{d};
            move_tbl=All_SpikesPerMove_Tbl(strcmp(All_SpikesPerMove_Tbl.MoveType,move_type)&contains(All_SpikesPerMove_Tbl.move_trial_ID,depth_code),:);
            rest_depth_tbl=rest_tbl(contains(rest_tbl.move_trial_ID,depth_code),:);
            if isempty(move_tbl)||isempty(rest_depth_tbl), continue; end
            
            % Compute PSTH/FR for move & REST
            [mFR,eFR,mPSTH,ePSTH]=computeCurves(move_tbl,edges_FR,edges_PSTH,binSize_FR,smoothPlots,useSEM);
            [mFRr,eFRr,mPSTHr,ePSTHr]=computeCurves(rest_depth_tbl,edges_FR,edges_PSTH,binSize_FR,smoothPlots,useSEM);
            pctChange=((mean(mFR)-mean(mFRr))/mean(mFRr))*100;
            baseline_summary(end+1,:)={CaseDate,CaseDate_hem,depth_name,move_type,height(move_tbl),height(rest_depth_tbl),pctChange};
            
            % Raster+PSTH+FR fig
            compName = sprintf('RESTvs%s_PSTH_FR_%s%s.png', matlab.lang.makeValidName(move_type), CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], ''));
            plotBaselinePSTH_FR(time_FR,time_PSTH,mFR,eFR,mFRr,eFRr,mPSTH,ePSTH,mPSTHr,ePSTHr,move_type,depth_name,depth_colors(d,:),rest_color,fontTitle,fontLabel,fullfile(figMultiDir,compName));
        end
    end
else
    fprintf('[INFO] No REST found. Skipping baseline comparisons.\n');
end

%% Export CSVs
summary_Table=cell2table(summary_Data,'VariableNames',{'CaseDate','Hemisphere','Depth','MoveType','N_Trials','Mean_FR_Hz',['Err_FR_' errLabel],'Mean_PSTH',['Err_PSTH_' errLabel]});

writetable(summary_Table, fullfile(csvDir, sprintf('FR_PSTH_Summary_%s%s.csv', CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], ''))));
full_Table=cell2table(full_Data,'VariableNames',{'CaseDate','Hemisphere','Depth','MoveType','TimeVector_FR','Mean_FR','Err_FR','TimeVector_PSTH','Mean_PSTH','Err_PSTH'});
writetable(full_Table, fullfile(csvDir, sprintf('FR_PSTH_FullBins_%s%s.csv', CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], ''))));
if ~isempty(baseline_summary)
    baseline_Table=cell2table(baseline_summary,'VariableNames',{'CaseDate','Hemisphere','Depth','MoveType','N_MoveTrials','N_RESTTrials','PctChange_FR'});
    writetable(baseline_Table, fullfile(csvDir, sprintf('Baseline_Comparison_Summary_%s%s.csv', CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], ''))));
end
fprintf('[INFO] CSV export complete: %s\n',csvDir);
end

%% ===== Helper Functions =====

function res=ternary(cond,a,b), if cond,res=a;else,res=b;end,end

function plotRasterPSTH(spikeTimes_all,time_FR,mean_FR,err_FR,time_PSTH,mean_PSTH,err_PSTH,move_type,depth_name,errLabel,fontTitle,fontLabel)
figure('Visible','off','Position',[100 100 900 700]); tiledlayout(3,1,'TileSpacing','compact');
nexttile; hold on; for t=1:numel(spikeTimes_all), plot(spikeTimes_all{t}, t*ones(size(spikeTimes_all{t})), 'k.', 'MarkerSize', 8); end
xlim([-0.05 0.45]); ylabel('Trial'); title(sprintf('Raster | %s | %s',depth_name,move_type),'FontSize',fontTitle);
nexttile; shadedArea(time_PSTH,mean_PSTH,err_PSTH,[0.6 0.6 0.9]); ylabel(sprintf('Spike Count (%s)',errLabel),'FontSize',fontLabel); title('PSTH','FontSize',fontTitle);
nexttile; shadedArea(time_FR,mean_FR,err_FR,[0.2 0.6 1]); xlabel('Time (s)','FontSize',fontLabel); ylabel(sprintf('FR (Hz %s)',errLabel),'FontSize',fontLabel); title('FR','FontSize',fontTitle);
end

function shadedArea(x,y,err,color)
fill([x fliplr(x)],[y+err fliplr(y-err)],color,'FaceAlpha',0.3,'EdgeColor','none'); hold on; plot(x,y,'Color',color*0.7,'LineWidth',2);
end

function plotStackedRaster(spikeTimes_byDepth,move_type,depth_labels,colors,savePath,fontTitle,fontLabel)
figure('Visible','off','Position',[100 100 900 700]); tiledlayout(3,1,'TileSpacing','compact');
for p=1:3, trials=spikeTimes_byDepth{p}; nexttile; hold on;
for t=1:numel(trials), plot(trials{t}, t*ones(size(trials{t})), '.', 'Color', colors(p,:), 'MarkerSize', 8); end
xlim([-0.1 0.4]); ylim([0 numel(trials)+1]); if p<3, set(gca,'XTickLabel',[]); else, xlabel('Time (s)','FontSize',fontLabel); end
ylabel('Rep','FontSize',fontLabel); title(depth_labels{p},'FontSize',fontTitle);
end
sgtitle(sprintf('Stacked Raster | %s',move_type),'FontSize',fontTitle,'FontWeight','bold'); saveas(gcf,savePath); savefig(strrep(savePath,'.png','.fig')); close;
end

function plotMultiMovementRaster(allDepthSpikes, movements, depth_labels, colors, restColor, savePath, fontTitle, fontLabel)
numMovements = numel(movements);
cols = numMovements; % Will be 4 if REST present
figure('Visible','off','Position',[100 100 1300 900]);
tiledlayout(3, cols, 'TileSpacing','compact');

for d = 1:3 % Rows = depths
    for m = 1:numMovements
        nexttile;
        trials = allDepthSpikes{d, m};
        if isempty(trials), continue; end
        hold on;

        % Determine color (REST = gray)
        if strcmp(movements{m}, 'REST')
            plotColor = restColor;
        else
            plotColor = colors(d,:);
        end

        for t2 = 1:numel(trials)
            plot(trials{t2}, t2*ones(size(trials{t2})), '.', 'Color', plotColor, 'MarkerSize', 8);
        end
        xlim([-0.1 0.4]); ylim([0 numel(trials)+1]);
        if d == 3, xlabel('Time (s)', 'FontSize', fontLabel); end
        if m == 1, ylabel(depth_labels{d}, 'FontSize', fontLabel); end
        title(movements{m}, 'FontSize', fontTitle);
    end
end

% Add Legend
lgd = legend({'REST','dorsal STN','central STN','ventral STN'}, ...
    'Orientation','horizontal','Location','southoutside','FontSize',fontLabel);
lgd.Box = 'off';

sgtitle('Multi-Movement Raster', 'FontSize', fontTitle, 'FontWeight', 'bold');

% Save figure
saveas(gcf, savePath);
savefig(strrep(savePath, '.png', '.fig'));
close;
end


function [mean_FR,err_FR,mean_PSTH,err_PSTH]=computeCurves(tbl,edges_FR,edges_PSTH,binSize_FR,smoothPlots,useSEM)
FR_all=[]; PSTH_all=[];
for i=1:height(tbl)
    spkTimes=(tbl.C1{i}-tbl.TTL_spk_idx_Start(i))/44000-0.05;
    FR_all=[FR_all; histcounts(spkTimes,edges_FR)/binSize_FR];
    PSTH_all=[PSTH_all; histcounts(spkTimes,edges_PSTH)];
end
mean_FR=mean(FR_all,1); mean_PSTH=mean(PSTH_all,1);
if useSEM, err_FR=std(FR_all,0,1)/sqrt(size(FR_all,1)); err_PSTH=std(PSTH_all,0,1)/sqrt(size(PSTH_all,1));
else, err_FR=std(FR_all,0,1); err_PSTH=std(PSTH_all,0,1); end
if smoothPlots, mean_FR=smoothdata(mean_FR,'gaussian',5); mean_PSTH=smoothdata(mean_PSTH,'gaussian',3); end
end

function plotBaselinePSTH_FR(time_FR,time_PSTH,mFR,eFR,mFRr,eFRr,mPSTH,ePSTH,mPSTHr,ePSTHr,...
    move_type,depth_name,moveColor,restColor,fontTitle,fontLabel,savePath)

figure('Visible','off','Position',[100 100 900 800]); 
tiledlayout(3,1,'TileSpacing','compact');

% ---- Raster Placeholder ----
nexttile;
text(0.5,0.5,'Raster Omitted','HorizontalAlignment','center','FontSize',12);
axis off;

% ---- PSTH Comparison ----
nexttile; hold on;
plot(time_PSTH,mPSTH,'Color',moveColor,'LineWidth',2);
plot(time_PSTH,mPSTHr,'Color',restColor,'LineWidth',2);
fillArea(time_PSTH,mPSTH,ePSTH,moveColor);
fillArea(time_PSTH,mPSTHr,ePSTHr,restColor);
ylabel('Spike Count','FontSize',fontLabel);
title('PSTH','FontSize',fontTitle);
legend({move_type,'REST'},'Location','northeast');

% ---- FR Comparison ----
nexttile; hold on;
plot(time_FR,mFR,'Color',moveColor,'LineWidth',2);
plot(time_FR,mFRr,'Color',restColor,'LineWidth',2);
fillArea(time_FR,mFR,eFR,moveColor);
fillArea(time_FR,mFRr,eFRr,restColor);
xlabel('Time (s)','FontSize',fontLabel);
ylabel('FR (Hz)','FontSize',fontLabel);
title('Firing Rate','FontSize',fontTitle);
legend({move_type,'REST'},'Location','northeast');

% ---- Global Title ----
sgtitle(sprintf('%s vs REST | %s',move_type,depth_name),'FontSize',fontTitle,'FontWeight','bold');

% ---- Save ----
saveas(gcf,savePath);
savefig(strrep(savePath,'.png','.fig'));
close;
end

function fillArea(x, y, err, color)
% fillArea: Creates shaded error bands
fill([x fliplr(x)], [y+err fliplr(y-err)], color, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
