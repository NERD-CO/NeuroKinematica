%% Order_of_Op__LFPKin_Processing

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

%% Config - Define offset duration for useOffset_LFP function

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

%% Config - Define epoch duration for UniformEpochs_LFP function

epochDur_ms = 500; % milliseconds
Epoch_dur_seconds = epochDur_ms / 1000; % seconds

% Calculate number of TTL samples
epochDur_TTLs = round(TTL_fs * Epoch_dur_seconds); % ensure value is integer

% Calculate number TTL samples in AO_LFP sample domain by downsampling TTL_fs
epochDur_LFPs = round((epochDur_TTLs/TTL_fs)*AO_LFP_fs); % ensure value is integer

UniformEpochs = true;
% If UniformEpochs == true or epochDur_ms > 0, useOffset_LFP function returns the #
% of samples to pre-pad in LFP domain based on a set pre-trial offset time
% (and a meta struct for reference).
% If UniformEpochs == false or epochDur_ms <= 0, useOffset_LFP function returns 0.


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


%% Config - Input and Output Data Dirs

Case_DataDir = [IO_DataDir, filesep, CaseDate];
MoveDataDir = [IO_DataDir, filesep, 'Processed DLC'];

% Case-specific Input dirs
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology'];       % directory for processed ephys data and spike clusters
Move_CaseDir = [MoveDataDir, filesep, Move_CaseID];                         % directory for processed DLC data and Movement Indices

% Case-specific Output dir
ephysTbl_Dir = [Case_DataDir, filesep, 'DLC_Ephys'];                        % directory where all ephys per move-rep tables are located


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


%% Run align_LFPsPerMove_TTL function

All_LFPsPerMove_Tbl = align_LFPsPerMove_TTL(Subject_AO, ProcDataDir, Move_CaseDir, ephysTbl_Dir, TTL_fs, AO_LFP_fs, pre_offset_ms, useOffset);

%% Run align_LFPsPerMove_uniformEpochDur function
% Add 500 ms from start of each MoveRep

All_LFPsPerMove_Tbl_uniformEpochs = align_LFPsPerMove_uniformEpochDur(Subject_AO, ProcDataDir, Move_CaseDir, ephysTbl_Dir, TTL_fs, AO_LFP_fs, pre_offset_ms, useOffset, UniformEpochs, epochDur_ms);


%% Load all LFPs per Movement Tbl and Offset_meta struct

cd(ephysTbl_Dir)
if useOffset && pre_offset_ms>0
    load(sprintf('All_LFPsPerMove_offset%ims.mat', pre_offset_ms), 'All_LFPsPerMove_Tbl');
else
    load('All_LFPsPerMove_N0offset.mat', 'All_LFPsPerMove_Tbl');
end

load('Offset_meta_struct.mat', 'meta_Offset')


%% Config - Inputs for spectrumInterpolation (Miguel's function :))

% This function interpolates around the frequency of interest (Fl) and
% replaces its and some neighbors using a constant value.

% function inputs
% data = All_LFPsPerMove_Tbl.LFPs; % data: column vector of the data that needs to be filtered
Fs = AO_LFP_fs; % Fs: Sampling Frequency (in Hz) of the data, % loop on harmonics of 60 Hz
Fl = 60; % Fl: Line frequency (in Hz), center of our interpolation
% loop on harmonics of 60 Hz (notch & comb)
neighborsToSample = 4; % Hz, 4 or 5
neighborsToReplace = 2; % Hz, 1 or 2

% neighborsToSample: This parameter is in Hz, and tells this function how
% large of a window (in Hz) to use when picking the constant value to
% replace frequency content with.

% neighborsToReplace: This parameter is also in Hz, and tells this function
% which neighbors need to be replaced with the constant value that was
% determined by neighborsToSample. Generally, neighborsToReplace <
% neighborsToSample in order to get a better spectral estimate.


% pick which columns to filter (LFP channels) ---
LFPsPerMoveTbl_vars = string(All_LFPsPerMove_Tbl.Properties.VariableNames);
lfpCols = LFPsPerMoveTbl_vars(startsWith(LFPsPerMoveTbl_vars,"LFP_E"));     % e.g., LFP_E1, LFP_E2, ...

% create new columns with "_filt" (set false to overwrite originals)
makeNewCols = true;


%% Run wrapper_applySpectrumInterpolation_LFPs function

% run wrapper function for spectrumInterpolation and output updated All_LFPsPerMove_Tbl
All_LFPsPerMove_Tbl_filt = wrapper_applySpectrumInterpolation_LFPs( All_LFPsPerMove_Tbl, ...
    lfpCols, AO_LFP_fs, Fl, neighborsToSample, neighborsToReplace, ...
    makeNewCols);

% save All_LFPsPerMove_Tbl_filt
cd(ephysTbl_Dir)

if useOffset && pre_offset_ms > 0
    outName = sprintf('filt_All_LFPsPerMove_offset%ims.mat', pre_offset_ms);
else
    outName = 'filt_All_LFPsPerMove_N0offset.mat';
end

save(outName, "All_LFPsPerMove_Tbl_filt");

%% filters

% 1) high-pass at 0.5 or 1 Hz
% low freq. noise + drift removal
% use designfilt + filtfilt

% 1) low-pass at 250 Hz
% Adjust: LP cutoff to ~200–230 Hz;
% anti-alias + keep LFP band
% use designfilt + filtfilt

% 3) downsample / resample 1375 to 500 Hz

% 4) Beta bandpass (13–30 Hz) after resampling

%% Question (in initial notes)

% 1-3 kHz (?)

%% Filtering (after spectrumInterpolation)

% Input:  All_LFPsPerMove_Tbl_filt (with LFP_E*_filt columns, Fs = 1375 Hz)
% Output: Preprocessed vectors (high-pass, low-pass, downsampled to 500 Hz)

Fs_in  = AO_LFP_fs;   % original LFP sampling rate = 1375 Hz
Fs_out = 500;         % downsampled rate

% Design filters (once)
hpFilt = designfilt('highpassiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency', 1, ...      % cutoff at ~1 Hz
    'SampleRate', Fs_in);

lpFilt = designfilt('lowpassiir', ...
    'FilterOrder', 8, ...
    'HalfPowerFrequency', 250, ...    % cutoff at 250 Hz (or ~200–230 Hz)
    'SampleRate', Fs_in);

% Pick which columns to filter
lfpCols_filt = string(All_LFPsPerMove_Tbl_filt.Properties.VariableNames);
lfpCols_filt = lfpCols_filt(startsWith(lfpCols_filt,"LFP_E") & endsWith(lfpCols_filt,"_filt"));

% Loop through each LFP column and row
for col_i = 1:numel(lfpCols_filt)
    colName = lfpCols_filt(col_i);
    outName = replace(colName,"_filt","_proc500");  % output column name

    % Preallocate
    All_LFPsPerMove_Tbl_filt.(outName) = cell(height(All_LFPsPerMove_Tbl_filt),1);

    for row_i = 1:height(All_LFPsPerMove_Tbl_filt)
        LFP_vec = All_LFPsPerMove_Tbl_filt.(colName){row_i};
        if isempty(LFP_vec) || ~isnumeric(LFP_vec) || ~isvector(LFP_vec)
            All_LFPsPerMove_Tbl_filt.(outName){row_i} = LFP_vec;
            continue
        end

        % --- Filtering steps
        LFP_vec = double(LFP_vec(:));                  % ensure column, double
        LFP_vec_hp = filtfilt(hpFilt, LFP_vec);        % high-pass drift removal
        LFP_vec_lp = filtfilt(lpFilt, LFP_vec_hp);     % low-pass anti-alias
        LFP_vec_ds = resample(LFP_vec_lp, 4, 11);      % downsample 1375 → 500 Hz

        % Store result
        All_LFPsPerMove_Tbl_filt.(outName){row_i} = LFP_vec_ds;
    end
end


%% Visualization (Quality check / QA)

testRow = 6;   % change index to preview other trials
rawVec   = All_LFPsPerMove_Tbl_filt.(lfpCols_filt(1)){testRow};
filtVec  = All_LFPsPerMove_Tbl_filt.(replace(lfpCols_filt(1),"_filt","_proc500")){testRow};

figure;
subplot(2,1,1);
plot(rawVec); title('Raw LFP (post spectrumInterpolation, 1375 Hz)');
subplot(2,1,2);
plot(filtVec); title('Preprocessed LFP (HP+LP+downsampled to 500 Hz)');








%% Cut these visualizations:
%% pspectrum PSD comparison: raw (post spectrumInterpolation) vs preprocessed (HP+LP+DS to 500 Hz)

% pick a row and an LFP column to visualize
rowIdx = testRow;  % adjust to inspect other trials
Fs_res = 2;  % frequency resolution (target bin size)

% grab the first filtered LFP column name, e.g., "LFP_E1_filt"
if isempty(lfpCols_filt)
    warning('No LFP_E*_filt columns found. Skipping PSD plots.');
else
    colFilt_1375 = lfpCols_filt(1);
    colProc_500   = replace(colFilt_1375,"_filt","_proc500");

    % Use char() for widest MATLAB compatibility with dynamic fieldnames
    lfp_x_raw = All_LFPsPerMove_Tbl_filt.(char(colFilt_1375)){rowIdx};   % Fs_raw = 1375
    lfp_x_proc = All_LFPsPerMove_Tbl_filt.(char(colProc_500)){rowIdx};   % Fs_proc = 500

    if isempty(lfp_x_raw) || isempty(lfp_x_proc) || ~isvector(lfp_x_raw) || ~isvector(lfp_x_proc)
        warning('Selected row/columns contain empty or non-vector data. Choose another row.');
    else
        lfp_x_raw = double(lfp_x_raw(:));
        lfp_x_proc = double(lfp_x_proc(:));

        minDur = numel(lfp_x_raw)/AO_LFP_fs;  % seconds of data
        % if duration < 2 s, relax target Fs_res
        if minDur < 2
            Fs_res = 5;
        end

        % Require finite samples
        if ~any(isfinite(lfp_x_raw)) || ~any(isfinite(lfp_x_proc))
            warning('Selected signals contain only non-finite values. Skipping PSD.');
        else
            Fs_raw = AO_LFP_fs;          % 1375
            Fs_proc = 500;                % downsampled
            fMax   = 250;                % max for band of interest
            epsF   = 1e-6;               % margin from Nyquist

            % Clip FrequencyLimits to Nyquist for each signal
            lim_raw = [0, min(fMax, Fs_raw/2 - epsF)];
            lim_pre = [0, min(fMax, Fs_proc/2 - epsF)];

            % PSDs
            % Try with fixed FrequencyResolution; if it fails, retry auto
            try
                [P_raw,F_raw] = pspectrum(lfp_x_raw, Fs_raw, ...
                    'FrequencyLimits', lim_raw, ...
                    'FrequencyResolution', Fs_res);
            catch
                warning('pspectrum(raw): falling back to automatic resolution.');
                [P_raw,F_raw] = pspectrum(lfp_x_raw, Fs_raw, 'FrequencyLimits', lim_raw);
            end

            try
                [P_proc,F_proc] = pspectrum(lfp_x_proc, Fs_proc, ...
                    'FrequencyLimits', lim_pre, ...
                    'FrequencyResolution', Fs_res);
            catch
                warning('pspectrum(preproc): falling back to automatic resolution.');
                [P_proc,F_proc] = pspectrum(lfp_x_proc, Fs_proc, 'FrequencyLimits', lim_pre);
            end

            % Convert power to dB
            P_raw_dB = pow2db(P_raw);
            P_proc_dB = pow2db(P_proc);

            % Plot overlay up to 100 Hz (detail around line noise & beta)
            figure('Name','PSD Raw vs Processed (0–100 Hz)');
            plot(F_raw, P_raw_dB, 'LineWidth', 1); hold on;
            plot(F_proc, P_proc_dB, 'LineWidth', 1);
            xlim([0 100]); grid on;
            xlabel('Frequency (Hz)'); ylabel('Power (dB)');
            title(sprintf('PSD (Row %d) — %s vs %s', rowIdx, char(colFilt_1375), char(colProc_500)), 'Interpreter','none');
            legend('Raw @ 1375 Hz','Processed @ 500 Hz','Location','best');

            % Plot overlay up to 250 Hz (full range of interest)
            figure('Name','PSD Raw vs Processed (0–250 Hz)');
            plot(F_raw, P_raw_dB, 'LineWidth', 1); hold on;
            plot(F_proc, P_proc_dB, 'LineWidth', 1);
            xlim([0 250]); grid on;
            xlabel('Frequency (Hz)'); ylabel('Power (dB)');
            title(sprintf('PSD (Row %d) — %s vs %s', rowIdx, char(colFilt_1375), char(colProc_500)), 'Interpreter','none');
            legend('Raw @ 1375 Hz','Processed @ 500 Hz','Location','best');

            % Optional: time-domain traces
            t_raw = (0:numel(lfp_x_raw)-1)/Fs_raw;
            t_pre = (0:numel(lfp_x_proc)-1)/Fs_proc;

            figure('Name','Time Series Raw vs Processed');
            subplot(2,1,1);
            plot(t_raw, lfp_x_raw); grid on;
            xlabel('Time (s)'); ylabel('Amplitude');
            title(sprintf('Raw (1375 Hz) — %s (Row %d)', char(colFilt_1375), rowIdx), 'Interpreter','none');
            subplot(2,1,2);
            plot(t_pre, lfp_x_proc); grid on;
            xlabel('Time (s)'); ylabel('Amplitude');
            title(sprintf('Processed (500 Hz) — %s (Row %d)', char(colProc_500), rowIdx), 'Interpreter','none');
        end
    end
end


%% Welch PSD comparison: raw (post spectrumInterpolation) vs preprocessed (HP+LP+DS to 500 Hz)

% RAW @ 1375 Hz; for ~1-2 Hz bin spacing, nfft = 1024 gives ~1.34 Hz bins
% PROC @ 500 Hz; for ~1-2 Hz bin spacing, nfft = 256 → ~1.95 Hz bins

% choose which row/columns to visualize
rowIdx = testRow;   % change to inspect another movement rep segment (row)

if isempty(lfpCols_filt)
    warning('No LFP_E*_filt columns found. Skipping Welch PSD.');
else
    colFilt_1375 = lfpCols_filt(1);                              % e.g., "LFP_E1_filt"
    colProc_500   = replace(colFilt_1375,"_filt","_proc500");    % e.g., "LFP_E1_proc500"

    % pull vectors (use char() for dynamic fieldname compatibility)
    lfp_x_raw = All_LFPsPerMove_Tbl_filt.(char(colFilt_1375)){rowIdx};   % Fs_raw = 1375
    lfp_x_proc = All_LFPsPerMove_Tbl_filt.(char(colProc_500)){rowIdx};   % Fs_pre = 500

    if isempty(lfp_x_raw) || isempty(lfp_x_proc) || ~isvector(lfp_x_raw) || ~isvector(lfp_x_proc)
        warning('Selected row has empty/non-vector data. Choose another row.');
    else
        lfp_x_raw = double(lfp_x_raw(:));
        lfp_x_proc = double(lfp_x_proc(:));

        % remove DC (optional but recommended)
        lfp_x_raw = lfp_x_raw - mean(lfp_x_raw,'omitnan');
        lfp_x_proc = lfp_x_proc - mean(lfp_x_proc,'omitnan');

        if ~any(isfinite(lfp_x_raw)) || ~any(isfinite(lfp_x_proc))
            warning('Signals contain non-finite values only. Skipping Welch PSD.');
        else
            % Welch parameters (auto-tune for ~1–2 Hz/bin)
            Fs_raw = AO_LFP_fs;   % 1375 Hz
            Fs_proc = 500;        % downsampled

            fs_target = 1.5;      % desired bin spacing (Hz)
            fMax      = 250;      % plot limit (Hz)

            % helper for nfft given Fs and target df
            computeNFFT = @(Fs,df) 2^nextpow2(max(128, round(Fs/df)));

            nfft_raw = computeNFFT(Fs_raw, fs_target);   % ~1.3–1.6 Hz/bin
            nfft_proc = computeNFFT(Fs_proc, fs_target);   % ~1.0–2.0 Hz/bin

            % window lengths: ~4 s (or smaller if signal short), >= nfft
            win_raw = min(numel(lfp_x_raw), max(nfft_raw, round(4*Fs_raw)));
            win_pre = min(numel(lfp_x_proc), max(nfft_proc, round(4*Fs_proc)));

            % ensure minimum practical sizes
            win_raw = max(256, win_raw);
            win_pre = max(128, win_pre);

            % if still longer than the signal, clip to signal length
            win_raw = min(win_raw, numel(lfp_x_raw));
            win_pre = min(win_pre, numel(lfp_x_proc));

            w_raw  = hann(win_raw,'periodic');
            w_pre  = hann(win_pre,'periodic');
            olap_raw = floor(0.5*win_raw);     % 50% overlap
            olap_pre = floor(0.5*win_pre);

            % Welch PSDs (one-sided)
            [P_raw, F_raw] = pwelch(lfp_x_raw, w_raw, olap_raw, nfft_raw, Fs_raw, 'onesided');
            [P_proc, F_proc] = pwelch(lfp_x_proc, w_pre, olap_pre, nfft_proc, Fs_proc, 'onesided');

            % keep up to fMax
            keep_raw = F_raw <= fMax;
            keep_proc = F_proc <= fMax;

            F_raw = F_raw(keep_raw);
            F_proc = F_proc(keep_proc);
            P_raw_dB = pow2db(P_raw(keep_raw));
            P_proc_dB = pow2db(P_proc(keep_proc));

            % Plot 0–100 Hz
            figure('Name','Welch PSD — Raw vs Processed (0–100 Hz)');
            plot(F_raw, P_raw_dB, 'LineWidth', 1); hold on;
            plot(F_proc, P_proc_dB, 'LineWidth', 1);
            xlim([0 100]); grid on;
            xlabel('Frequency (Hz)'); ylabel('Power (dB)');
            title(sprintf('Welch PSD (Row %d) — %s vs %s', ...
                rowIdx, char(colFilt_1375), char(colProc_500)), 'Interpreter','none');
            legend('Raw @ 1375 Hz','Processed @ 500 Hz','Location','best');

            % Plot 0–250 Hz (full range of interest)
            figure('Name','Welch PSD — Raw vs Processed (0–250 Hz)');
            plot(F_raw, P_raw_dB, 'LineWidth', 1); hold on;
            plot(F_proc, P_proc_dB, 'LineWidth', 1);
            xlim([0 250]); grid on;
            xlabel('Frequency (Hz)'); ylabel('Power (dB)');
            title(sprintf('Welch PSD (Row %d) — %s vs %s', ...
                rowIdx, char(colFilt_1375), char(colProc_500)), 'Interpreter','none');
            legend('Raw @ 1375 Hz','Processed @ 500 Hz','Location','best');
        end
    end
end


%% Save fully filtered and processed LFP outputs (HP+LP+DS to 500 Hz)

% Resave All_LFPsPerMove_Tbl_filt
cd(ephysTbl_Dir)

if useOffset && pre_offset_ms > 0
    outName = sprintf('filtProc_All_LFPsPerMove_offset%ims.mat', pre_offset_ms);
else
    outName = 'filtProc_All_LFPsPerMove_N0offset.mat';
end

save(outName, "All_LFPsPerMove_Tbl_filt");
fprintf('[SAVED] Processed table written to %s\n', fullfile(ephysTbl_Dir,outName));
