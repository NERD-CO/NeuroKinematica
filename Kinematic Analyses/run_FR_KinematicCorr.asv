function run_FR_KinematicCorr(CaseDate, CaseDate_hem, ephys_offset, MoveDir_CaseID)

% run_FR_KinematicCorr  Correlate per‐trial firing rates with kinematic features.
%
%   run_FR_KinematicCorr(CaseDate,CaseDate_hem,ephys_offset) loads your
%   spike‐per‐move table (All_SpikesPerMove_Tbl) alongside all per‐trial
%   kinematic summary CSVs in the "Kinematic Analyses\CaseDate_hem" folder,
%   merges them, applies outlier handling and optional transforms, computes
%   correlations per movement type & depth, plots scatter+fits, and writes
%   a master CSV with all trial‐level data plus correlation metrics.

clearvars -except CaseDate CaseDate_hem ephys_offset; clc;

MoveDir_CaseID = 'IO_03_23_2023_LSTN';

%% ==== TOGGLE OPTIONS ====

useActualDuration = false;      % false = fixed window_FR Option B
window_FR         = [-0.05 0.45]; % FR window (s) for Option B
AO_spike_fs       = 44000;      % sampling rate

outlierStrategy   = 'remove';   % 'remove' | 'flag'
logTransformFR    = false;      % log‐transform FR if heavy‐tailed
logTransformKin   = false;      % log‐transform kin features
zscoreFR          = false;      % z‐score FR before correlation
zscoreKin         = false;      % z‐score kin features
useSpearman       = false;      % fallback to Spearman if non‐normal

% update to calibration-image based method later
pixels_to_mm      = 2.109;      % distance conversion factor

% % Load checkerboard image
% I = imread('calibration_frame.png');
% [imagePoints, boardSize] = detectCheckerboardPoints(I);
% squareSize_mm = 25.4;
% pixels_to_mm = squareSize_mm; % assuming 1 square in image = 25.4 mm
% 
% % Or estimate from image width:
% avg_square_px = mean(diff(imagePoints(:,1)));
% pixels_to_mm = squareSize_mm / avg_square_px;


%% Directory Setup

curPC = getenv('COMPUTERNAME');
switch curPC
    case 'DESKTOP-I5CPDO7'
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    otherwise
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
end

Case_DataDir       = fullfile(IO_DataDir, CaseDate);
ephysTbl_Dir       = fullfile(Case_DataDir,'DLC_Ephys');
MoveDataDir        = fullfile(IO_DataDir, 'Processed DLC', MoveDir_CaseID);

isBilateral = contains(CaseDate,'bilateral','IgnoreCase',true);
if isBilateral
    fprintf('[INFO] Bilateral case detected: %s\n',CaseDate);
    CaseDate_hem = input('Enter hemisphere (LSTN or RSTN): ','s');
    if ~ismember(CaseDate_hem,{'LSTN','RSTN'}), error('Invalid hemisphere'); end
    ephysTbl_Dir = fullfile(ephysTbl_Dir, CaseDate_hem);
    % MoveDataDir  = fullfile(MoveDataDir, CaseDate_hem);
else
    CaseDate_hem = '';
end
fprintf('[INFO] Loading spike data from: %s\n', ephysTbl_Dir);
fprintf('[INFO] Loading DLC data from: %s\n', MoveDataDir);

output_kinAnalysisDir     = fullfile(IO_DataDir,'Kinematic Analyses',[CaseDate '_' CaseDate_hem]);
output_CorrDir      = fullfile(ephysTbl_Dir,'FR_Kinematics');
if ~exist(output_CorrDir,'dir'), mkdir(output_CorrDir); end

%% Load spike‐per‐trial table

cd(ephysTbl_Dir);
mats = dir('*Spikes*.mat');
names = {mats.name};
if ephys_offset
    matfile = names{contains(names,'offset') & contains(names,'Spikes')};
else
    matfile = names{~contains(names,'offset') & contains(names,'Spikes')};
end
load(matfile,'All_SpikesPerMove_Tbl');


%% Compute per‐trial FR (one row per move_trial_ID)

% Option B: fixed window_FR for all trials
spkTbl = All_SpikesPerMove_Tbl;
[uniqMovT,uniqMovT_indices] = unique(spkTbl.move_trial_ID,'stable');
frTbl = table(uniqMovT, 'VariableNames',{'move_trial_ID'});
fr = nan(size(uniqMovT));
for uniqMovT_indices=1:numel(uniqMovT)
    rowIdx = find(strcmp(spkTbl.move_trial_ID,uniqMovT{uniqMovT_indices}));
    % collect all spikes for that trial (C1 times)
    spikes = cell2mat(cellfun(@(c) c - spkTbl.TTL_spk_idx_Start(rowIdx(1)), ...
        spkTbl.C1(rowIdx),'uni',0)) / AO_spike_fs - window_FR(1);
    fr(uniqMovT_indices) = numel(spikes)/diff(window_FR);
end
frTbl.FR_Hz = fr;

% check frTbl columns
disp('frTbl columns:');
disp(frTbl.Properties.VariableNames);

%% Load all kinematic‐summary CSVs into one table

% MoveDataDir subfolders:
Move_csvs = [MoveDataDir, filesep, 'csv folder']; % contains kinematic timeseries data for each marker (raw DLC output)
Move_mats = [MoveDataDir, filesep, 'mat folder']; % contains csv folder data files converted to .mat
Move_vids = [MoveDataDir, filesep, 'video folder']; % contains DLC-labeled motor trial videos and Movement Index CSVs

kinFiles = dir(fullfile(Move_vids,'*.csv'));
kinTbls = cellfun(@(f) readtable(fullfile(Move_vids,f)), ...
    {kinFiles.name}, 'Uni',false);
kinTbl  = vertcat(kinTbls{:});  % assumes each has move_trial_ID + numeric columns

% check if a table
assert(istable(kinTbl), 'kinTbl is not a table');

% Reconstruct move_trial_ID in kinTbl
% Example format: 't3_c' = trial 3 at central depth
kinTbl.move_trial_ID = cell(height(kinTbl),1);

% Extract depth indicator (e.g., 't', 'c', or 'b') from Move_CaseID
if contains(MoveDir_CaseID,'_LSTN')
    hemPrefix = 'L';
elseif contains(MoveDir_CaseID,'_RSTN')
    hemPrefix = 'R';
else
    hemPrefix = '';
end

% Infer depth char from Move_CaseID
depth_map = {'t','c','b'}; % dorsal, central, ventral
depth_char = '';
for i = 1:numel(depth_map)
    if contains(MoveDir_CaseID, depth_map{i})
        depth_char = depth_map{i};
        break
    end
end

% Construct move_trial_ID
for i = 1:height(kinTbl)
    trialNum = kinTbl.MoveN(i);
    kinTbl.move_trial_ID{i} = sprintf('%s%d_%s', depth_char, trialNum, depth_char);
end

% check kinTbl columns
disp('kinTbl columns:');
disp(kinTbl.Properties.VariableNames);

%% Merge into frTbl and kinTbl into master table

masterTbl = innerjoin(frTbl, kinTbl, 'Keys','move_trial_ID');

% check masterTbl columns
disp('masterTbl columns:');
disp(masterTbl.Properties.VariableNames);

%% Optional transforms & outlier handling

% --- Outliers ---
if strcmpi(outlierStrategy,'remove')
    ok = all(~isoutlier(masterTbl{:,2:end}),2);
    masterTbl = masterTbl(ok,:);
end
% --- Log‐transform if flagged ---
if logTransformFR
    masterTbl.FR_Hz = log(masterTbl.FR_Hz + eps);
end
if logTransformKin
    kinCols = setdiff(masterTbl.Properties.VariableNames,{'move_trial_ID','FR_Hz'});
    masterTbl{:,kinCols} = log(masterTbl{:,kinCols}+eps);
end
% --- Z‐score if flagged ---
if zscoreFR
    masterTbl.FR_Hz = zscore(masterTbl.FR_Hz);
end
if zscoreKin
    kinCols = setdiff(masterTbl.Properties.VariableNames,{'move_trial_ID','FR_Hz'});
    masterTbl{:,kinCols} = zscore(masterTbl{:,kinCols});
end

%% Correlation by MovementType × Depth

% Assume kinTbl has MovementType & Depth variables
movTypes = unique(masterTbl.MovementType);
depths   = unique(masterTbl.Depth);

corrResults = table;
for m = 1:numel(movTypes)
    for d = 1:numel(depths)
        sel = strcmp(masterTbl.MovementType,movTypes{m}) & strcmp(masterTbl.Depth,depths{d});
        sub = masterTbl(sel,:);
        if height(sub)<3, continue; end

        % pick one kinematic feature, e.g. 'MeanDisplacement'
        kinFeat = 'MeanDisplacement';
        x = sub.FR_Hz;
        y = sub.(kinFeat);

        % test normality
        if ~useSpearman && (adtest(x) || adtest(y))
            [R,p] = corr(x,y,'Type','Spearman');
            method='Spearman';
        else
            [R,p] = corr(x,y,'Type','Pearson');
            method='Pearson';
        end

        % effect size = Fisher z
        eff = atanh(R);

        % save
        newRow = {movTypes{m}, depths{d}, height(sub),R,p,eff};
        corrResults = [corrResults;
            cell2table(newRow,'VariableNames',{'MovementType','Depth','N','R','p','FisherZ'})];

        % scatter + fit plot
        fig = figure('Visible','off');
        scatter(x,y,50,'filled'); hold on
        lsline;
        xlabel('Firing rate (Hz)'); ylabel(kinFeat);
        title(sprintf('%s | %s: R=%.2f, p=%.3f (%s)', ...
            depths{d},movTypes{m},R,p,method));
        saveas(fig, fullfile(output_CorrDir, ...
            sprintf('%s_%s_corr.png',movTypes{m},depths{d})));
        close(fig);
    end
end

%% Write master CSV and correlation table
writetable(masterTbl, fullfile(output_CorrDir,'Master_FR_Kin.csv'));
writetable(corrResults, fullfile(output_CorrDir,'Correlation_Summary.csv'));

fprintf('Done! Master CSV and correlation results in:\n%s\n',output_CorrDir);
end
