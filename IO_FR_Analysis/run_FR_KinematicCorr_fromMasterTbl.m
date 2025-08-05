function run_FR_KinematicCorr_fromMasterTbl(agg_FRKinCorr_dir, output_agg_FRKinCorr_dir, masterTbl)

% Accepts full merged FR-Kinematic table instead of computing it internally
% Computes correlation between FR and kinematic features from merged table

% INPUTS:
%   masterTbl : table with columns: MoveType, Depth, mean_FR_Hz, kinematic features...
%   outputDir : directory to save correlation plots and CSV summary

agg_FRKinCorr_dir = [Ephys_Kin_Dir, filesep, 'aggregated_FRKinCorr'];
if ~exist(agg_FRKinCorr_dir, 'dir')
    mkdir(agg_FRKinCorr_dir);
end

% Create output directory if needed
output_agg_FRKinCorr_dir = [Ephys_Kin_Dir, filesep, 'aggregated_FRKinCorr', filesep, 'v1'];
if ~exist(output_agg_FRKinCorr_dir, 'dir')
    mkdir(output_agg_FRKinCorr_dir);
end

cd(agg_FRKinCorr_dir)
masterTblPath = fullfile(agg_FRKinCorr_dir, 'master_FRKin_Tbl_MANUAL.csv');
masterTbl = readtable(masterTblPath);

%% TOGGLE OPTIONS 

outlierStrategy   = 'remove';   % 'remove' | 'flag'
logTransformFR    = false;
logTransformKin   = false;
zscoreFR          = false;
zscoreKin         = false;
useSpearman       = false;      % fallback not used currently


%% Transforms and Outlier Handling

numericVars = varfun(@isnumeric, masterTbl, 'OutputFormat', 'uniform');
numericCols = masterTbl(:, numericVars);

if strcmpi(outlierStrategy, 'remove')
    ok = all(~isoutlier(numericCols{:,:}), 2);
    masterTbl = masterTbl(ok, :);
end

if logTransformFR
    masterTbl.Mean_FR_Hz = log(masterTbl.Mean_FR_Hz + eps);
end
if logTransformKin
    masterTbl.MeanAmp = log(masterTbl.MeanAmp + eps);
    masterTbl.MeanVel = log(masterTbl.MeanVel + eps);
end
if zscoreFR
    masterTbl.Mean_FR_Hz = zscore(masterTbl.Mean_FR_Hz);
end
if zscoreKin
    masterTbl.MeanAmp = zscore(masterTbl.MeanAmp);
    masterTbl.MeanVel = zscore(masterTbl.MeanVel);
end


%% Correlation analysis

movTypes = unique(masterTbl.MoveType);
depths   = unique(masterTbl.Depth);
kinFeatures = {'MeanAmp', 'MeanVel'};

AggCorrResults = table([], [], [], [], [], [], [], [], ...
    'VariableNames', {'MovementType','Depth','KinematicFeature','N','R','p','FisherZ','Method'});

for m = 1:numel(movTypes)
    for d = 1:numel(depths)
        sel = strcmp(masterTbl.MoveType, movTypes{m}) & strcmp(masterTbl.Depth, depths{d});
        sub = masterTbl(sel, :);
        if height(sub) < 3
            continue;
        end
        for k = 1:numel(kinFeatures)
            kinFeat = kinFeatures{k};
            x = sub.Mean_FR_Hz;
            y = sub.(kinFeat);
            if all(isnan(y)), continue; end

            [R, p] = corr(x, y, 'Type', 'Pearson');
            eff = atanh(R); % Fisher z
            method = 'Pearson';

            newRow = {movTypes{m}, depths{d}, kinFeat, height(sub), R, p, eff, method};
            AggCorrResults = [AggCorrResults;
                cell2table(newRow, 'VariableNames', ...
                {'MovementType','Depth','KinematicFeature','N','R','p','FisherZ','Method'})];

            % Save scatter plot
            fig = figure('Visible','off');
            scatter(x, y, 50, 'filled'); hold on;
            lsline;
            xlabel('Firing Rate (Hz)');
            ylabel(kinFeat);
            title(sprintf('%s | %s | %s: R = %.2f, p = %.3f', ...
                depths{d}, movTypes{m}, kinFeat, R, p));
            saveas(fig, fullfile(output_agg_FRKinCorr_dir, ...
                sprintf('%s_%s_%s_corr.png', movTypes{m}, depths{d}, kinFeat)));
            close(fig);
        end
    end
end

%% Save output

writetable(AggCorrResults, fullfile(output_agg_FRKinCorr_dir, 'Aggregate_Correlation_Summary.csv'));

fprintf('[DONE] Aggregate Correlation summary written to: %s\n', output_agg_FRKinCorr_dir);

end
