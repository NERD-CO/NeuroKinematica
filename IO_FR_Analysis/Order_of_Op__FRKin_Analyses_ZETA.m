%% Order of Operations - FR-Kinematic Correlation Script

clear; close all; clc;


%% Functions in FR Analysis and FR-Kinematic Correlation Pipeline

% MaxSpkDuration_Raster_PSTH

% Zeta Testing functions
% 1) zetatest.m
% % > Helper:
% % > Wrapper: 
% 2) getIFR.m

% plot_Raster_PSTH

% compute_FRperMove_perSTNdepth_v3

% run_IO_FR_Analysis_and_Plotting

% run_MovementFeatureAnalysis_IO_v2

% merge_FRKin_SummaryTbls

% aggregate_FRKinematic_Correlations


%% Directory Setup

curPCname = getenv('COMPUTERNAME');
switch curPCname
    case 'DESKTOP-I5CPDO7'  % PC_1
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'DSKTP-JTLAB-EMR'  % Lab Desktop
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\DLC_Processing'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\zetatest'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\Kinematic Analyses'
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_LFP_Analysis'
    case 'NSG-M-H8J3X34'    % PC_2
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        addpath 'C:\GitHub\NeuroKinematica\DLC_Processing'
        addpath 'C:\GitHub\NeuroKinematica\IO_FR_Analysis'
        addpath 'C:\GitHub\NeuroKinematica\zetatest'
        addpath 'C:\GitHub\NeuroKinematica\Kinematic Analyses'
        addpath 'C:\GitHub\NeuroKinematica\IO_LFP_Analysis'
end
cd(IO_DataDir)
Subject_Hem_CaseMap = readtable('Subject_Hem_MetaSummary.xlsx'); % Subject_Hem_MetaSummary.xlsx


%% Config - Define datastream sampling rates (Alpha Omega and Video sampling fs)

TTL_fs = 44000;         % Hz, Alpha Omega TTL clock
AO_spike_fs = 44000;    % Hz, Alpha Omega spike sampling rate
AO_LFP_fs = 1375;       % Hz, Alpha Omega LFP sampling rate
DLC_fs = 100;           % fps, Video/DLC frame rate


%% Config - Define pre-trial offset duration for useOffset_spikes function

pre_offset_ms = 50; % milliseconds
offset_seconds = pre_offset_ms / 1000; % seconds

% Calculate number of TTL samples
offset_TTLs = round(TTL_fs * offset_seconds); % ensure value is integer

% Calculate number TTL samples in AO_LFP sample domain by downsampling TTL_fs
offset_spikes = round((offset_TTLs/TTL_fs)*AO_spike_fs); % ensure value is integer

useOffset = true;
% If useOffset == true or pre_offset_ms>0, useOffset_spikes function returns the #
% of samples to pre-pad in spike sample domain based on a set pre-trial offset time
% (and a meta struct for reference).
% If useOffset == false or pre_offset_ms<=0, useOffset_spike function returns 0.


%% Inputs:

% Ephys data folder:
CaseDate = '11_30_2023_bilateral'; % Adjust as needed

% 1: '03_23_2023';             % NER 2025
% 1: '04_05_2023';             % NER 2025

% 2: '04_13_2023_bilateral';   % GRC 2026      %% errors - Unrecognized field name "dblZetaT".

% 3: '05_18_2023_b_bilateral'; % NER 2025      %% errors using new pipeline --> check MIs; 
%   %   % Other error: Unrecognized field name "dblZetaT".

% 4: '05_31_2023';             % INS 2026

% 5: '06_08_2023_bilateral';   % NER 2025      %% errors using new pipeline --> check MIs; 
%   %   % Other error: Error using findpeaks>parse_inputs (line 225): Data set must contain at least 3 samples.

% 6: '07_06_2023_bilateral';   % INS 2026
% 7: '07_13_2023_bilateral';   % INS 2026
% 10: '08_23_2023';             % INS 2026
% 12: '11_30_2023_bilateral';   % GRC 2026


% Kinematic data folder:
MoveDir_CaseID = 'IO_11_30_2023_RSTN'; % Adjust as needed

% 'IO_03_23_2023_LSTN';   % NER 2025
% 'IO_04_05_2023_RSTN';   % NER 2025
% 'IO_04_13_2023_LSTN';   % GRC 2026
% 'IO_05_18_2023_b_LSTN'; % NER 2025
% 'IO_05_31_2023_LSTN';   % INS 2026
% 'IO_06_08_2023_LSTN';   % NER 2025
% 'IO_07_06_2023_LSTN';   % INS 2026
% 'IO_07_13_2023_LSTN';   % INS 2026
% 'IO_07_13_2023_RSTN';   % GRC 2026
% 'IO_08_23_2023_RSTN';   % INS 2026
% 'IO_11_30_2023_LSTN';   % GRC 2026
% 'IO_11_30_2023_RSTN';   % GRC 2026


%% Data folder paths

Case_DataDir = fullfile(IO_DataDir, CaseDate);
% ephysTbl_Dir = fullfile(Case_DataDir,'DLC_Ephys');
MoveDataDir = fullfile(IO_DataDir, 'Processed DLC');
KinematicsDir = fullfile(IO_DataDir, 'Kinematic Analyses');
Ephys_Kin_Dir = fullfile(IO_DataDir, 'Ephys_Kinematics');
FR_Kin_Dir = fullfile(Ephys_Kin_Dir, 'FR_Kinematic_Analyses');
Case_FRKin_Dir = fullfile(FR_Kin_Dir, CaseDate); % input dir (new ephysTbleDir) and output (results) dir


%% Handle Bilateral Cases

% Specify hemisphere in command window if bilateral
isBilateral = contains(CaseDate,'bilateral','IgnoreCase',true);
if isBilateral
    fprintf('[INFO] Bilateral case detected: %s\n',CaseDate);
    CaseDate_hem = input('Enter hemisphere (LSTN or RSTN): ','s');
    if ~ismember(CaseDate_hem,{'LSTN','RSTN'})
        error('Invalid hemisphere');
    end
    % Append hemisphere-specific folder
    Case_FRKin_Dir = fullfile(Case_FRKin_Dir, CaseDate_hem);
else
    CaseDate_hem = '';
end
fprintf('[INFO] Loading spike data from: %s\n', Case_FRKin_Dir);


%% Load All_SpikesPerMove_Tbl

cd(Case_FRKin_Dir)
Tbl_list = dir('*Spikes*.mat');
Tbl_names = {Tbl_list.name};
spk_case = Tbl_names{contains(Tbl_names, 'offset')}; % offset version preferred
load(spk_case, 'All_SpikesPerMove_Tbl');


%% OPTIONAL: Case-specific cleaning (remove duplicates if needed)

% for 3_23_2023 case - remove duplicates / only plot for primary electrode
if strcmp(char(CaseDate), '03_23_2023')
    % All_SpikesPerMove_Tbl = All_SpikesPerMove_Tbl(158:end,1:13); % Comment or adjust as needed
    All_SpikesPerMove_Tbl = All_SpikesPerMove_Tbl(168:end,1:13);
elseif strcmp(char(CaseDate), '04_05_2023')
    All_SpikesPerMove_Tbl = [All_SpikesPerMove_Tbl(1:68,1:11); All_SpikesPerMove_Tbl(133:197,1:11); All_SpikesPerMove_Tbl(266:329,1:11)];
end


%% ===== Functions =====

%% Run MaxSpkDuration_Raster_PSTH

% outputs max spike segment duration and a raster + PSTH for all trials
[Max_SpikeDuration_samples, spikesMatrix] = MaxSpkDuration_Raster_PSTH(CaseDate, All_SpikesPerMove_Tbl, AO_spike_fs, Case_FRKin_Dir);

Max_SpkDur_seconds = Max_SpikeDuration_samples/AO_spike_fs; % seconds
Max_SpkDus_ms = Max_SpkDur_seconds * 1000; % milliseconds

%% Zeta Test Functions

% https://github.com/JorritMontijn/zetatest?tab=readme-ov-file

% Point to folder that contains zetatest.m and the dependencies/ subfolder
zetaRepo = 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\zetatest';
addpath(genpath(zetaRepo));   % genpath adds dependencies/
rehash

% % quick test
% mv = 'HAND OC'; dz = 't';
% move_tbl = All_SpikesPerMove_Tbl(strcmp(All_SpikesPerMove_Tbl.MoveType,mv) & contains(All_SpikesPerMove_Tbl.move_trial_ID,dz),:);
% 
% [spkT, evTimes, useMaxDur] = makeZetaInputs_fromAOStartStopTimes( ...
%     move_tbl, AO_spike_fs, ...
%     'PreWindow_s', 0.050, ...
%     'PostWindow_s', 0.000);
% 
% [pZ,sZ] = zetatest(spkT, evTimes, useMaxDur, 500, 0);
% fprintf('Sanity OK: %s-%s  p=%.3g  Z=%.2f\n', mv, dz, pZ, sZ.dblZETA);


%% Run Zeta test for each MoveType × STN depth using true per-trial durations

% https://elifesciences.org/articles/71969
% https://github.com/JorritMontijn/zetatest/blob/main/zetatest.m
% zetatest.m: Calculates the Zenith of Event-based Time-locked Anomalies (ZETA) for spike times of a single neuron. Outputs a p-value.
%   % binning-free method for determining whether a neuron shows any time-locked modulation of spiking activity
%   % detects whether a single neuron is responsive to a stimulus/event in a statistically robust way and avoids binning and parameter selection

% run wrapper + helper function for zetatest.m
% loop through all possible electrodes
[ZETA_Summary, all_sZETA, all_sRate, all_sLatencies] = runZETA_byDepthMove_actualDurations( ...
    All_SpikesPerMove_Tbl, AO_spike_fs, ...
    'UseMaxDur_s',  1.60, ...     % 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8
    'Resamples',    5000, ...
    'PlotFlag',     0, ...
    'RestrictRange',[-inf inf], ...
    'DirectQuantile', false, ...
    'JitterSize',    2, ...
    'Stitch',        true, ...
    'PreWindow_s',   0.050, ...   % 50 ms lead-in
    'PostWindow_s',  0.000);


%% Save Outputs to "Zeta Testing" folder

Zeta_outDir = fullfile(Case_FRKin_Dir, 'Zeta Testing');
if ~exist(Zeta_outDir,'dir'); mkdir(Zeta_outDir); end

ZetaSummary_csv = fullfile(Zeta_outDir, sprintf('%s_ZETA_Summary.csv', CaseDate));
ZetaAll_mat = fullfile(Zeta_outDir, sprintf('%s_ZETA_AllOutputs.mat', CaseDate));

writetable(ZETA_Summary, ZetaSummary_csv);
save(ZetaAll_mat, 'ZETA_Summary', 'all_sZETA', 'all_sRate', 'all_sLatencies', '-v7.3');

fprintf('ZETA test outputs saved:\n  %s\n  %s\n', ZetaSummary_csv, ZetaAll_mat);


%% Zeta test outputs of interest (p-vals and z-vals):

%	- dblZetaP; p-value based on Zenith of Event-based Time-locked Anomalies
%   % Low ZETA-test p-values indicate that the neuron's firing pattern is 
%   % statistically unlikely to be observed if the neuron is not modulated 
%   % by the event of interest.

%	- sZETA; structure with fields:
	%	- dblZETA; responsiveness z-score (i.e., >2 is significant), ZETA (ζ) 
    %       % use p with the standard normal's quantile function Φ−1 to obtain a corrected ZETA, ζ, that is interpretable as a z-score
    %       % z-score: ζ = Φ−1(1−p2) 
	%	- dblD; temporal deviation value underlying ZETA
    %       % d, a time-invariant mean-normalized version of δ, a neuron's deviation from a temporally non-modulated spiking rate at a time point
    %       % Zenith of Event-based Time-locked Anomalies (ZETA) = most extreme value, that is the maximum of the absolute values: max(|d|)
    %       % Statistical significance of ZETA ... extreme value theory ... distribution of maximum values is known as a Gumbel distribution (Gumbel, 1941)
	%	- dblP; p-value corresponding to ZETA 
    %       % cumulative Gumbel distribution at sample maximum, ζr: 
    %       % p-value: p = 1−F(ζr;m,β)  ...    % same as dblZetaP
	%	- dblZetaT; time corresponding to ZETA
    % ...

%   - vecLatencies; different latency estimates, number determined by intLatencyPeaks.
    % ...
%	- sRate; structure with fields:
    % ...
%   - sLatencies; structure containing latencies (copy of sZETA.vecLatencies)
    % ...


%% Compute instantaneous firing rate (IFR) + PSTH per MoveType × STN depth using getIFR 
% uses the temporal deviations upon which ZETA is based 
% d: temporal deviation vector 
% δ: scaled to lie between –1 and +1, value depends on the difference between the fractional position of a spike (from 0 to 1) and the linear interval from x,y=[0,0] to [τ,1]. 

% https://github.com/JorritMontijn/zetatest/blob/main/getIFR.m
% getIFR.m: Calculates the instantaneous firing rate (IFR) without running the ZETA-test. Use this as you would a PSTH function.

IFR_outDir = fullfile(Zeta_outDir, 'IFR_PSTH');  

[IFR_PSTH_Summary, all_IFR] = runIFR_PSTH_byDepthMove(All_SpikesPerMove_Tbl, AO_spike_fs, ...
    'UseMaxDur_s', 1.60, ...
    'PadITI_s',    0.005, ...
    'SpikeField',  'C1', ...
    'StartField',  'TTL_spk_idx_Start', ...
    'StopField',   'TTL_spk_idx_End', ...
    'PreWindow_s', 0.050, ...
    'PostWindow_s',0.000, ...
    'BinSize_ms',  10, ...
    'IFR_SmoothSd',2, ...
    'IFR_MinScale',[], ...
    'IFR_Base',    1.5, ...
    'DoPlot',      true, ...
    'SaveDir',     IFR_outDir, ...
    'CaseDate',    CaseDate);

% change electrode when necessary when calling function 
% 'SpikeField', 'C1'
% or loop through all possible electrodes


% Save summary table + all structs
if ~exist(IFR_outDir,'dir'); mkdir(IFR_outDir); end
writetable(IFR_PSTH_Summary, fullfile(IFR_outDir, sprintf('%s_IFR_PSTH_Summary.csv', CaseDate)));
save(fullfile(IFR_outDir, sprintf('%s_IFR_PSTH_All.mat', CaseDate)), 'IFR_PSTH_Summary', 'all_IFR', '-v7.3');
fprintf('IFR/PSTH saved to %s\n', IFR_outDir);

cd(IFR_outDir)


%% Plot ZETA z-score scatter plots
                    
MasterZETA = aggregate_ZETA_and_plot(FR_Kin_Dir, ...
    'IO_DataDir', IO_DataDir, ...
    'ZetaCsvNamePattern','*_ZETA_Summary.csv', ...
    'SigZ', 2, 'SigP', 0.05, ...
    'YMax', 5); 


%% PCA of ZETA z-scores per Subject / Spike Field (C1, C2, ..., Cn) for each MoveType per STN Depth

% replicate London et al., 2021 paper's Fig 2
% PCA on spike response profiles (ZetaZ scores)
% run through MasterZeta

% Category: each MoveType per STN depth (E.g., Hand OC x dorsal STN)

% x-dim (rows) = time (in PSTH bins), must be conserved across subjects
% y-dim (cols) = subject neuron/spike field/unit

% distinct PCA per category


%% zetatstest.m: Calculates the time-series version of ZETA, for data such as calcium imaging or EEG recordings.

% run on all_IFR per MoveType x STN depth

%% MUA

%% Heat plot 


%%
