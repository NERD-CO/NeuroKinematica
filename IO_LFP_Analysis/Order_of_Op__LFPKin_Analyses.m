%% Order_of_Op__LFPKin_Analyses

clear; close all; clc;

%% Environment / Directory Set-up

% specify directory where case-specific data files are located
curPCname = getenv('COMPUTERNAME');

switch curPCname
    case 'DESKTOP-I5CPDO7'  % PC_1
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'DSKTP-JTLAB-EMR'  % Lab Desktop
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_LFP_Analysis'
    case 'NSG-M-H8J3X34'    % PC_2
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        addpath 'C:\GitHub\NeuroKinematica\IO_LFP_Analysis'
end

cd(IO_DataDir)
Subject_AO = readtable('Subject_AO.xlsx');


%% Config - Define datastream sampling rates (Alpha Omega and Video sampling fs)

TTL_fs = 44000;         % Hz, Alpha Omega TTL clock
AO_spike_fs = 44000;    % Hz, Alpha Omega spike sampling rate
AO_LFP_fs = 1375;       % Hz, Alpha Omega LFP sampling rate
DLC_fs = 100;           % fps, Video/DLC frame rate


%% Config - Define offset duration and useOffset_LFP function

pre_offset_ms = 50; % milliseconds
offset_seconds = pre_offset_ms / 1000; % seconds

% Calculate number of TTL samples
offset_TTLs = round(TTL_fs * offset_seconds); % ensure value is integer

% Calculate number TTL samples in AO_LFP sample domain by downsampling TTL_fs
offset_LFPs = round((offset_TTLs/TTL_fs)*AO_LFP_fs); % ensure value is integer

useOffset = true;
% If useOffset == true or pre_offset_ms>0, useOffset_LFP function returns the #
% of samples to pre-pad in LFP domain based on a set pre-trial offset time
% (and a meta struct for reference).
% If useOffset == false or pre_offset_ms<=0, useOffset_LFP function returns 0.


%% Config - Ephys Case Input

% isolate a specific CaseDate / studyID (StudyNum in Subject_AO csv)
CaseDate = '03_23_2023';

% '03_09_2023'; % studyID = 1, ptID 1

% '03_23_2023'; % studyID = 2, ptID 2    * % Use for INS 2026
% '04_05_2023'; % studyID = 3, ptID 2    * % Use for INS 2026

% '04_13_2023_bilateral'; ptID 3
% studyID = 4(L), 5(R),

% '05_11_2023'; % studyID = 6, ptID 4
% '05_18_2023_a'; % studyID = 7, ptID 4

% '05_18_2023_b_bilateral';
% LSTN: studyID = 8, ptID = 5    % Use for INS 2026
% RSTN: studyID = 9, ptID = 5

% '05_31_2023';  % studyID = 10, ptID 6

% '06_08_2023_bilateral'; ptID = 7
% LSTN: studyID = 11,
% RSTN: studyID = 12(R),

% '07_13_2023_bilateral';
% studyID = 15(L), 16(R), ptID = 9


%% Config - Movement Case Input

% Specify case ID
Move_CaseID = 'IO_03_23_2023_LSTN';

% 'IO_03_09_2023_RSTN'; % studyID = 1, ptID 1 (processed, incomplete case)

% 'IO_03_23_2023_LSTN'; % studyID = 2, ptID 2 (processed, complete case) *
% 'IO_04_05_2023_RSTN'; % studyID = 3, ptID 2 (processed, complete case) *

% 'IO_04_13_2023_LSTN'; % studyID = 4, ptID 3 (processed, complete case)
% 'IO_04_13_2023_RSTN'; % studyID = 5, ptID 3

% 'IO_05_11_2023_LSTN'; % studyID = 6, ptID 4 (processed, incomplete case)
% 'IO_05_18_2023_a_RSTN'; % studyID = 7, ptID 4

% 'IO_05_18_2023_b_LSTN'; % studyID = 8, ptID 5 (processed, complete case) *
% 'IO_05_18_2023_b_RSTN'; % studyID = 9, ptID 5

% 'IO_05_31_2023_LSTN'; % studyID = 10, ptID 6

% 'IO_06_08_2023_LSTN'; % studyID = 11, ptID = 7 (processed, complete case)
% 'IO_06_08_2023_RSTN'; % studyID = 12, ptID = 7 (processed, incomplete case)

% 'IO_07_13_2023_LSTN'; % studyID = 15, ptID = 9
% 'IO_07_13_2023_RSTN'; % studyID = 16, ptID = 9

ephys_offset = 1;


%% Config - Input and Output Data Dirs

Case_DataDir = [IO_DataDir, filesep, CaseDate];
MoveDataDir = [IO_DataDir, filesep, 'Processed DLC'];

% Case-specific Input dirs
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology'];       % directory for processed ephys data and spike clusters
Move_CaseDir = [MoveDataDir, filesep, Move_CaseID];                         % directory for processed DLC data and Movement Indices

% Case-specific Output dir
ephysTbl_Dir = [Case_DataDir, filesep, 'DLC_Ephys'];                        % directory where all ephys per move-rep tables are located

% Data folder paths
KinematicsDir = fullfile(IO_DataDir, 'Kinematic Analyses');
Case_KinDir = fullfile(KinematicsDir, Move_CaseID);
Ephys_Kin_Dir = fullfile(IO_DataDir, 'Ephys_Kinematics');
LFP_Kin_Dir = fullfile(Ephys_Kin_Dir, 'LFP_Kinematic_Analyses');
Case_LFP_Kin = fullfile(LFP_Kin_Dir, CaseDate); % case-specific results from run_FR_KinematicCorr saved here


%% Config - Handle bilateral cases and hemisphere selection

% create metadata sheet to pull from (for future users)

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
else
    CaseDate_hem = ''; % No hemisphere for unilateral cases
end

% Append hemisphere folder if needed
if ~isempty(CaseDate_hem)
    ProcDataDir = fullfile(ProcDataDir, CaseDate_hem);
    fprintf('[INFO] Hemisphere-specific input ephys directory set: %s\n', ProcDataDir);
    ephysTbl_Dir = fullfile(ephysTbl_Dir, CaseDate_hem);
    fprintf('[INFO] Hemisphere-specific output directory set: %s\n', ephysTbl_Dir);

else
    fprintf('[INFO] Using base ephys input directory: %s\n', ProcDataDir);
    fprintf('[INFO] Using base output directory: %s\n', ephysTbl_Dir);
end


%% Load all LFPs per Movement Tbl and Offset_meta struct

cd(ephysTbl_Dir)

if useOffset && pre_offset_ms>0
    load(sprintf('filtProc_All_LFPsPerMove_offset%ims.mat', pre_offset_ms), 'All_LFPsPerMove_Tbl_filt');
else
    load('filtProc_All_LFPsPerMove_N0offset.mat', 'All_LFPsPerMove_Tbl_filt');
end


%% Extract processed LFP columns

% Extract LFP column names
allVars_LFPtbl     = string(All_LFPsPerMove_Tbl_filt.Properties.VariableNames);
lfpCols_raw = allVars_LFPtbl(startsWith(allVars_LFPtbl,"LFP_E") & endsWith(allVars_LFPtbl,"_filt"));     % 1375 Hz, post spectrumInterpolation
lfpCols_proc= allVars_LFPtbl(startsWith(allVars_LFPtbl,"LFP_E") & endsWith(allVars_LFPtbl,"_proc500"));  % 500 Hz, fully preprocessed

if isempty(lfpCols_raw) || isempty(lfpCols_proc)
    warning('Expected LFP_E*_filt and/or LFP_E*_proc500 not found.');
end













%% Welch PSDs (save freq & power per trial & channel) for FOOOF/specparam

% raw (post spectrumInterpolation) vs processed (HP+LP+DS to 500 Hz)
% Robust to short trials by capping window/overlap/NFFT per trial.

% Define which row to visualize
rowIdx = 1; % testRow

% Pick a raw column and its matching processed column by name
colRaw_1375 = lfpCols_raw(1);                 % e.g., "LFP_E1_filt"
colProc_500 = replace(colRaw_1375,"_filt","_proc500");

% Make sure the processed partner exists
if ~ismember(colProc_500, lfpCols_proc)
    error('Processed partner column %s not found for %s.', colProc_500, colRaw_1375);
end

lfp_x_raw  = All_LFPsPerMove_Tbl_filt.(char(colRaw_1375)){rowIdx};   % 1375 Hz
lfp_x_proc = All_LFPsPerMove_Tbl_filt.(char(colProc_500)){rowIdx};   % 500 Hz

% Parameters
Fs_raw     = AO_LFP_fs;   % 1375
Fs_proc    = 500;         % 500
fMax       = 250;
fs_target  = 1.5;
win_s      = 4;
win_min    = 128;
nfft_min   = 128;
overlapFrac= 0.5;

% Compute Welch safely for each
[P_raw,F_raw]   = computeWelchSafe(lfp_x_raw,  Fs_raw,  fMax, fs_target, win_s, win_min, nfft_min, overlapFrac);
[P_proc,F_proc] = computeWelchSafe(lfp_x_proc, Fs_proc, fMax, fs_target, win_s, win_min, nfft_min, overlapFrac);

% Plot
P_raw_dB  = pow2db(P_raw);
P_proc_dB = pow2db(P_proc);

figure('Name','Welch PSD — Raw vs Processed (0–100 Hz)');
plot(F_raw, P_raw_dB, 'LineWidth', 1); hold on;
plot(F_proc, P_proc_dB, 'LineWidth', 1);
xlim([0 100]); grid on;
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
title(sprintf('Welch PSD (Row %d) — %s vs %s', rowIdx, char(colRaw_1375), char(colProc_500)), 'Interpreter','none');
legend('Raw @ 1375 Hz','Processed @ 500 Hz','Location','best');

figure('Name','Welch PSD — Raw vs Processed (0–250 Hz)');
plot(F_raw, P_raw_dB, 'LineWidth', 1); hold on;
plot(F_proc, P_proc_dB, 'LineWidth', 1);
xlim([0 250]); grid on;
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
title(sprintf('Welch PSD (Row %d) — %s vs %s', rowIdx, char(colRaw_1375), char(colProc_500)), 'Interpreter','none');
legend('Raw @ 1375 Hz','Processed @ 500 Hz','Location','best');


%% Band isolation & Hilbert transformation (append column metrics)

% theta band
% Low beta bandpass (13–20 Hz)
% High beta bandpass (21-35 Hz)
% gamma band

% Uses _proc500 columns as input (Fs = 500). Creates:
%   *_<band>_bp       : band-passed signal
%   *_<band>_amp      : instantaneous amplitude (abs(hilbert))
%   *_<band>_pow      : instantaneous power (amp.^2)
%   *_<band>_zamp     : z-scored amplitude (per-trial, optional but handy)

Fs_proc = 500;

%  bands of interest 
bands = struct( ...
    'theta',   [4 7], ...
    'L_beta',   [13 20], ...
    'H_beta',   [21 35], ...
    'gamma',   [36 90], ...
    'hfo',     [130 200] ...
);

%  find processed LFP columns 
if isempty(lfpCols_proc)
    warning('No LFP_E*_proc500 columns found. Run preprocessing first.');
else
    fprintf('[INFO] Band-extraction from %d processed LFP columns.\n', numel(lfpCols_proc));
end

% Design bandpass filter (4th order Butterworth IIR, zero-phase later)
BPfilt = @(lo,hi) designfilt('bandpassiir', 'FilterOrder', 4, ...
                           'HalfPowerFrequency1', lo, ...
                           'HalfPowerFrequency2', hi, ...
                           'SampleRate', Fs_proc);

% Pre-build filters to avoid re-alloc
bp = struct();
fieldNames = fieldnames(bands);
for k = 1:numel(fieldNames)
    curBand = bands.(fieldNames{k});
    % clamp to Nyquist
    curBand(2) = min(curBand(2), Fs_proc/2 - 1);
    if curBand(1) >= curBand(2)-1
        error('Band %s limits invalid after clamping.', fieldNames{k});
    end
    bp.(fieldNames{k}) = BPfilt(curBand(1), curBand(2));
end

% create output columns (cells)
nRows = height(All_LFPsPerMove_Tbl_filt);
for col_i = 1:numel(lfpCols_proc)
    base = erase(lfpCols_proc(col_i), "_proc500");  % e.g., "LFP_E1"
    for k = 1:numel(fieldNames)
        band = fieldNames{k};
        All_LFPsPerMove_Tbl_filt.(base + "_" + band + "_bp")   = cell(nRows,1);
        All_LFPsPerMove_Tbl_filt.(base + "_" + band + "_amp")  = cell(nRows,1);
        All_LFPsPerMove_Tbl_filt.(base + "_" + band + "_pow")  = cell(nRows,1);
        All_LFPsPerMove_Tbl_filt.(base + "_" + band + "_zamp") = cell(nRows,1);
    end
end


% Hilbert transformation @ spec. freq. (inst. power)
% process row-by-row

for col_i = 1:numel(lfpCols_proc)
    inCol = lfpCols_proc(col_i);
    base  = erase(inCol, "_proc500");

    for curBand = 1:nRows
        LFPvec_x = All_LFPsPerMove_Tbl_filt.(inCol){curBand};
        if isempty(LFPvec_x) || ~isvector(LFPvec_x) || ~isnumeric(LFPvec_x), continue; end
        LFPvec_x = double(LFPvec_x(:));
        if ~any(isfinite(LFPvec_x)), continue; end

        for k = 1:numel(fieldNames)
            band = fieldNames{k};
            try
                xbp = filtfilt(bp.(band), LFPvec_x);               % bandpassed
            catch ME
                warning('Row %d, %s (%s): bandpass failed: %s', ...
                    curBand, inCol, band, ME.message);
                continue
            end
            xa   = abs(hilbert(xbp));                       % inst. amplitude
            xp   = xa.^2;                                   % inst. power
            xz   = (xa - mean(xa,'omitnan')) ./ std(xa,0,'omitnan'); % z-amp

            All_LFPsPerMove_Tbl_filt.(base + "_" + band + "_bp"){curBand}   = xbp;
            All_LFPsPerMove_Tbl_filt.(base + "_" + band + "_amp"){curBand}  = xa;
            All_LFPsPerMove_Tbl_filt.(base + "_" + band + "_pow"){curBand}  = xp;
            All_LFPsPerMove_Tbl_filt.(base + "_" + band + "_zamp"){curBand} = xz;
        end
    end
end

fprintf('[DONE] Added bandpassed, amplitude, power, and z-amplitude columns.\n');


%% Save updated All_LFPsPerMove_Tbl_filt with computed characteristics per band

% Resave All_LFPsPerMove_Tbl_filt
cd(ephysTbl_Dir)

if useOffset && pre_offset_ms > 0
    outName = sprintf('bp_All_LFPsPerMove_offset%ims.mat', pre_offset_ms);
else
    outName = 'bp_All_LFPsPerMove_N0offset.mat';
end

save(outName, "All_LFPsPerMove_Tbl_filt");
fprintf('[SAVED] Bandpassed LFPsPerMove table written to %s\n', fullfile(ephysTbl_Dir, outName));


%% Run run_MovementFeatureAnalysis_IO_v2 or load kinTbl & kinSummaryTbl

cd(KinematicsDir)

% fprintf('[INFO] Loading movement data from: %s\n', MoveDataDir);
% [kinTbl, kinSummaryTbl] = run_MovementFeatureAnalysis_IO_v2(IO_DataDir, MoveDataDir, Move_CaseID);


%% PSDs on different context

% Extract LFP segments per 
% Rest, H O/C, A Pro/Sup, E Flex/Exten






%% ===== Export per-trial spectra for FOOOF/specparam =====
Fs = 500; fMax = 100; fs_target=1.5; win_s=4; win_min=128; nfft_min=128; overlapFrac=0.5;

out = struct();
out.meta.case   = CaseDate;
out.meta.moveID = Move_CaseID;

chanList = cellstr(procCols);  % list of channels
out.F = [];                    % frequency base (common, if possible)
out.P_by_trial = cell(numel(chanList), 1);

for col_i = 1:numel(chanList)
    col = chanList{col_i};
    Fref = []; Pcell = cell(nRows,1);
    for curBand = 1:nRows
        LFPvec_x = All_LFPsPerMove_Tbl_filt.(col){curBand};
        if isempty(LFPvec_x) || ~isvector(LFPvec_x) || ~isnumeric(LFPvec_x), continue; end
        [P,F] = computeWelchSafe(LFPvec_x, Fs, fMax, fs_target, win_s, win_min, nfft_min, overlapFrac);
        if isempty(P), continue; end
        if isempty(Fref)
            Fref = F;
        elseif numel(F) ~= numel(Fref) || max(abs(F-Fref)) > 1e-9
            P = interp1(F, P, Fref, 'linear', 'extrap');
        end
        Pcell{curBand} = P;
    end
    out.P_by_trial{col_i} = Pcell;
    if isempty(out.F), out.F = Fref; end
end

save(fullfile(ephysTbl_Dir, sprintf('Welch_forFOOOF_%s_%s.mat', CaseDate, Move_CaseID)), '-struct', 'out');
fprintf('[SAVED] Welch spectra for FOOOF/specparam.\n');

%% FOOOFs on different contexts

% * The aperiodic exponent of subthalamic field potentials reflects excitation/inhibition balance in Parkinsonism: https://pubmed.ncbi.nlm.nih.gov/36810199/; https://elifesciences.org/articles/82467#data
% * Parameterizing neural power spectra into periodic and aperiodic components: https://www.nature.com/articles/s41593-020-00744-x
% * Electrophysiological Frequency Band Ratio Measures Conflate Periodic and Aperiodic Neural Activity: https://www.eneuro.org/content/7/6/ENEURO.0192-20.2020
% * Electrophysiological Frequency Band Ratio Measures Conflate Periodic and Aperiodic Neural Activity: https://onlinelibrary.wiley.com/doi/abs/10.1111/ejn.15361
% * Time-resolved parameterization of aperiodic and periodic brain activity: Time-resolved parameterization of aperiodic and periodic brain activity | eLife
% * Subthalamic nucleus encoding steers adaptive therapies for gait in Parkinson's disease: https://www.medrxiv.org/content/10.1101/2025.08.20.25333478v1



%% Hilbert transformation spec. freq. (inst. power)

% theta band
% low-beta band
% high-beta band
% gamma band
% HFO

% at specific movement timepoints

%% Spike-field coherence

% inst. beta power at specific time/context (move index)

%% LFP bursting analysis / continuous wavelet transformation


%% Helper functions

% computeWelchSafe
function [P,F] = computeWelchSafe(x, Fs, fMax, df_target, win_s, win_min, nfft_min, overlapFrac)

% Robust Welch PSD for vectors of any length.
% Falls back to a periodogram when the segment would exceed x.

% Inputs:
%   x           : signal vector (row/col ok)
%   Fs          : sampling rate (Hz)
%   fMax        : upper frequency to keep (Hz)
%   df_target   : desired bin spacing ~Fs/nfft (Hz)
%   win_s       : target window length (seconds), e.g., 4
%   win_min     : minimum window length in samples, e.g., 128
%   nfft_min    : minimum NFFT, e.g., 128
%   overlapFrac : fraction overlap, e.g., 0.5
%
% Outputs:
%   P, F : one-sided Welch PSD and frequency vector (limited to <= fMax)

    % shape, finite, and DC removal
    x = double(x(:));
    x = x(isfinite(x));
    if isempty(x)
        P = []; F = [];
        return
    end
    x  = x - mean(x,'omitnan');
    Lx = numel(x);

    % choose window length that does not exceed the signal length
    win_len = min(Lx, max(win_min, round(win_s * Fs)));

    % if the signal is still too short to support a Welch segment, fallback
    if Lx < win_min
        nfft = max(nfft_min, 2^nextpow2(Lx));
        [P,F] = periodogram(x, [], nfft, Fs, 'onesided');
        keep = F <= fMax;
        P = P(keep); F = F(keep);
        return
    end

    % window, overlap, NFFT
    w      = hann(win_len, 'periodic');
    nover  = floor(overlapFrac * win_len);
    nover  = min(nover, win_len - 1);  % must be < win_len
    nfft   = max(nfft_min, 2^nextpow2(max(win_len, round(Fs/df_target))));
    nfft   = max(nfft, win_len);       % MATLAB allows nfft >= win_len

    % Welch PSD
    [P,F] = pwelch(x, w, nover, nfft, Fs, 'onesided');

    % limit to fMax
    keep = F <= fMax;
    P = P(keep); F = F(keep);
end


