function [FR_SummaryTbl] = run_IO_FR_Analysis_and_Plotting(CaseDate, CaseDate_hem, ephysTbl_Dir, ephys_offset, FR_Kin_Dir)

% Combined script for original + rescaled raster plots with mean ± SD FR in titles
% Adds p-values and significance stars in subtitles for REST vs Movement plots
% Includes Option A (actual duration) and Option B (fixed window) toggle for FR


%% Inputs

% CaseDate
% CaseDate_hem,
% ephys_offset

% Example:
% IO_plotting_vCombined('03_23_2023', '', 1);                 % unilateral
% IO_plotting_vCombined('05_18_2023_b_bilateral', 'LSTN', 1); % bilateral


%% CONFIG:

AO_spike_fs = 44000; % MER sampling rate

% FR Calculation Logic
useFullDuration = false; % true = Option A, false = Option B
window_FR = [-0.05 0.45];  % Fixed window for Option B

% Raster
binSize_FR = 0.01;
% edges_FR = window_FR(1):binSize_FR:window_FR(2);

depth_ids = {'t','c','b'};
depth_labels = {'dorsal STN','central STN','ventral STN'};
depth_colors = [0 0.6 0; 0.6 0 0.8; 0 0.4 0.8];
rest_color = [0.5 0.5 0.5];

sorted_moveTypes = {'HAND OC','HAND PS','ARM EF','REST'};
active_movements = {'HAND OC','HAND PS','ARM EF'};

fontTitle = 14; fontLabel = 12;


%% Output Directories

figDir = fullfile(FR_Kin_Dir, filesep, CaseDate, filesep, 'Figures'); if ~exist(figDir,'dir'), mkdir(figDir); end
figMultiDir = fullfile(figDir,'Raster_MultiMovement_and_RestComparison'); if ~exist(figMultiDir,'dir'), mkdir(figMultiDir); end
csvDir = fullfile(FR_Kin_Dir, filesep, CaseDate, filesep,'CSV_Outputs'); if ~exist(csvDir,'dir'), mkdir(csvDir); end


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


%% OPTIONAL: Case-specific cleaning (remove duplicates if needed)

% Example: All_SpikesPerMove_Tbl(158:end,:) for case 03_23_2023

% for 3_23_2023 case - remove duplicates / only plot for primary electrode
if strcmp(char(CaseDate), '03_23_2023')
    % All_SpikesPerMove_Tbl = All_SpikesPerMove_Tbl(158:end,1:13); % Comment or adjust as needed
    All_SpikesPerMove_Tbl = All_SpikesPerMove_Tbl(168:end,1:13);
elseif strcmp(char(CaseDate), '04_25_2023')
    All_SpikesPerMove_Tbl = [All_SpikesPerMove_Tbl(1:68,1:11); All_SpikesPerMove_Tbl(133:197,1:11); All_SpikesPerMove_Tbl(266:329,1:11)];
end


%% Identify move types

move_types = intersect(sorted_moveTypes, unique(All_SpikesPerMove_Tbl.MoveType),'stable');
hasREST = any(strcmp(move_types,'REST'));


%% Initialize storage

FR_SummaryTbl = {}; % FR summary storage (similar to kinSummaryTbl)

allDepthSpikes = cell(3, numel(active_movements));
restSpikes = cell(3,1);
rest_vs_move_stats = {};


%% Loop through trials

% Loop through trials and compute FR summary statistics
for m = 1:numel(move_types)
    move_type = move_types{m};
    for depth_i = 1:numel(depth_ids)
        depth_code = depth_ids{depth_i};
        depth_name = depth_labels{depth_i};

        % Extract the relevant trials for this MoveType and depth
        move_tbl = All_SpikesPerMove_Tbl(strcmp(All_SpikesPerMove_Tbl.MoveType, move_type) & ...
            contains(All_SpikesPerMove_Tbl.move_trial_ID, depth_code), :);
        if isempty(move_tbl), continue; end

        % Initialize a container to store FR for each MoveTrialID
        FR_per_MoveTrialID = containers.Map;

        % Initialize container to store spkTimes list
        spike_list = {};

        % Loop through the trials in move_tbl (each trial is identified by move_trial_ID)
        for move_row = 1:height(move_tbl)
            % Get the MoveTrialID for the current trial (e.g., 't1', 'c1', etc.)
            moveTrialID = move_tbl.move_trial_ID{move_row};

            % Get spike times for this trial (convert spike times as needed)
            spkTimes = (move_tbl.C1{move_row} - move_tbl.TTL_spk_idx_Start(move_row)) / AO_spike_fs - 0.05;
            spike_list{end+1} = spkTimes;

            % Compute the firing rate for this trial (using the helper function `computeFR`)
            FR = computeFR(spkTimes, useFullDuration, window_FR);

            % Store the computed FR for this trial under its MoveTrialID
            if isKey(FR_per_MoveTrialID, moveTrialID)
                FR_per_MoveTrialID(moveTrialID) = [FR_per_MoveTrialID(moveTrialID), FR];  % Append FR for this MoveTrialID
            else
                FR_per_MoveTrialID(moveTrialID) = FR;  % Initialize with FR for this MoveTrialID
            end

            % For REST vs. Move Stats
            if strcmp(move_type,'REST')
                restSpikes{depth_i} = spike_list;
            elseif ismember(move_type, active_movements)
                idx = find(strcmp(active_movements, move_type));
                allDepthSpikes{depth_i,idx} = spike_list;
            end
        end

        % Compute the mean and std of FR for each unique MoveTrialID
        for moveTrialID = keys(FR_per_MoveTrialID)
            % Get the FRs for this MoveTrialID across all repetitions of each MoveType
            FRs = FR_per_MoveTrialID(moveTrialID{1});

            % Compute the mean and standard deviation of FR for this MoveTrialID
            mu_FR = mean(FRs);
            sd_FR = std(FRs);

            % Append this trial's FR summary to the FRSummaryTbl (including MoveTrialID)
            FR_SummaryTbl(end+1, :) = {CaseDate, CaseDate_hem, depth_name, moveTrialID{1}, move_type, mu_FR, sd_FR};
        end
    end
end

%% Plot Original + Rescaled Multi-Movement Raster

multiRasterName = fullfile(figMultiDir, sprintf('MultiMovement_Raster_%s%s.png', CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], '')));
plotMultiMovementRaster(allDepthSpikes, active_movements, depth_labels, depth_colors, multiRasterName, fontTitle, fontLabel, useFullDuration, window_FR);
plotMultiMovementRasterRescaled(allDepthSpikes, active_movements, depth_labels, depth_colors, multiRasterName, fontTitle, fontLabel, useFullDuration, window_FR);
fprintf('[INFO] Figure export complete: %s\n', figDir);


%% REST vs Movement Comparison

if hasREST
    for m = 1:numel(active_movements)
        move_type = active_movements{m};
        figName = fullfile(figMultiDir, sprintf('REST_vs_%s_%s%s.png', matlab.lang.makeValidName(move_type), CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], '')));

        % Compute p-values and Cohen's d for REST vs Movement comparison
        [pVals, cohenD_vals] = plotRestVsMoveRaster(restSpikes, allDepthSpikes(:,m), move_type, depth_labels, depth_colors, rest_color, figName, fontTitle, fontLabel, useFullDuration, window_FR);
        plotRestVsMoveRasterRescaled(restSpikes, allDepthSpikes(:,m), move_type, depth_labels, depth_colors, rest_color, figName, fontTitle, fontLabel, useFullDuration, window_FR);

        % Store comparison stats (N_MoveTrials, N_RESTTrials, p-value, Cohen's d)
        for depth_i = 1:3
            rest_vs_move_stats(end+1,:) = {CaseDate, CaseDate_hem, depth_labels{depth_i}, move_type, length(allDepthSpikes{depth_i,m}), length(restSpikes{depth_i}), pVals(depth_i), cohenD_vals(depth_i)};
        end
    end
end


%% Export CSVs

% Convert FR_SummaryTbl to a table with appropriate column names
FR_SummaryTbl = cell2table(FR_SummaryTbl, 'VariableNames', ...
    {'CaseDate', 'Hemisphere', 'Depth', 'MoveTrialID', 'MoveType', 'Mean_FR_Hz', 'SD_FR_Hz'});

% Export FR summary table to CSV
writetable(FR_SummaryTbl, fullfile(csvDir, sprintf('FR_Summary_%s%s.csv', CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], ''))));


%% Convert REST vs. Move Stats to a table with appropriate column names

if ~isempty(rest_vs_move_stats)
    statsTable = cell2table(rest_vs_move_stats, 'VariableNames', {'CaseDate','Hemisphere','Depth','MoveType','N_MoveTrials','N_RESTTrials','p_value','Cohens_d'});

    % Export statsTable to CSV
    writetable(statsTable, fullfile(csvDir, sprintf('RESTvsMove_Stats_%s%s.csv', CaseDate, ternary(~isempty(CaseDate_hem), ['_' CaseDate_hem], ''))));
end

fprintf('[INFO] CSV export complete: %s\n', csvDir);

end




%% === Helper Functions ===

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


function res = ternary(cond,a,b), if cond,res=a;else,res=b;end,end

%% === Check Normality ===

function isNormal = checkNormality(data)
% Shapiro-Wilk is not available natively, so we'll use kstest
% If sample size < 4, we’ll assume normal for stability

data = data(~isnan(data));  % Remove NaNs 
    if numel(data) < 4
        warning('[checkNormality] Not enough valid data points (n = %d). Assuming non-normal.', numel(data));
        isNormal = false;
        return;
    end

    % Prevent division by 0 during z-scoring
    sigma = std(data);
    if sigma == 0 || isnan(sigma)
        warning('[checkNormality] Standard deviation is zero or NaN. Cannot z-score. Assuming non-normal.');
        isNormal = false;
        return;
    end

    data = (data - mean(data)) / sigma;  % z-score
    if isempty(data) || all(isnan(data))
        warning('[checkNormality] Data became empty or NaN after z-scoring. Assuming non-normal.');
        isNormal = false;
        return;
    end

    [h, ~] = kstest(data);  % h=0 → normal
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


%% === Plotting: MultiMovement (Rescaled) ===
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
        if all(isnan(FR_rest)) || all(isnan(FR_move))
            warning('[plotRestVsMoveRaster] All FRs are NaN for MoveType %s, Depth %s. Skipping stats.', ...
                move_type, depth_labels{d});
            continue;
        end

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


%% === Plotting: REST vs Movement (Rescaled) ===

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
        if all(isnan(FR_rest)) || all(isnan(FR_move))
            warning('[plotRestVsMoveRaster] All FRs are NaN for MoveType %s, Depth %s. Skipping stats.', ...
                move_type, depth_labels{d});
            continue;
        end

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
