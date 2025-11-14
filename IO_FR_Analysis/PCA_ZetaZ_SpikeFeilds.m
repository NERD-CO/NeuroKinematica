%% PCA of ZETA z-scores per MoveType × STN Depth - Subject Spike Fields (Clusters: C1, C2, ..., Cn)
% clear; clc;

% replicate London et al., 2021 paper's Fig 2 
% PCA on spike response profiles (ZetaZ scores)


%% Dimensions

% x-dim (rows) = time (in PSTH bins, 10ms), must be conserved across subjects
% y-dim (cols) = subject neuron/spike field/unit


% X: fixed time axis for PCA (10-ms bins) comes from runIFR_PSTH_byDepthMove 
%    via the PSTH centers it computes (centers, saved both into 
%    IFR_PSTH_Summary.PSTH_TimeCenters_s and each all_IFR{k}.centers). 
%    Use those as the PCA x-axis.

%% Navigate to MasterZeta file location

FR_Kin_Dir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\Ephys_Kinematics\FR_Kinematic_Analyses';
Aggr_ZETA_dir = fullfile(FR_Kin_Dir, 'Aggregate Zeta Plots');
cd(Aggr_ZETA_dir);

% Load master ZETA file (produced by aggregate_ZETA_and_plot)
MasterZETA = readtable('MasterZeta_AllSubjects.csv');
MasterZETA.Properties.VariableNames;


%% Extract metadata
% Extract ZetaZ values from MasterZeta file for each MoveType × STN depth 
% across subject spike clusters

ZetaZ     = MasterZETA.ZetaZ;       % [timeBins × units]
MoveTypes = MasterZETA.MoveType;    % 1 × units cell array
Depths    = MasterZETA.Depth;       % 1 × units cell array

uniqMoves = unique(MoveTypes, 'stable');
uniqDepths = unique(Depths, 'stable');


% run distinct PCA per category
% Category: each MoveType per STN depth (E.g., Hand OC x dorsal STN)

% Plot PC1 per movement category in tiled layout form 
% subplot row: STN depth, color code based on movement category 
% (like in the '...AllCategories_ByDepth_Tiles') plot created by the aggregate_ZETA_and_plot function)


%%




