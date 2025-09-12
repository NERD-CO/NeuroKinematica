%% Order_of_Op__LFPKin_Analyses

% clear; clc;
close all;

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
fs_downsamp = 500;      % Hz, Downsampled LFP to 500 Hz
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

%% Config - Define epoch duration for UniformEpochs_LFP function

epochDur_ms = 1000; % milliseconds
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

% Data folder paths
KinematicsDir = fullfile(IO_DataDir, 'Kinematic Analyses');
Case_KinDir = fullfile(KinematicsDir, Move_CaseID);
Ephys_Kin_Dir = fullfile(IO_DataDir, 'Ephys_Kinematics');
LFP_Kin_Dir = fullfile(Ephys_Kin_Dir, 'LFP_Kinematic_Analyses');
% Case_LFP_Kin = fullfile(LFP_Kin_Dir, CaseDate); % case-specific results from run_FR_KinematicCorr saved here


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


%% Load LFPsPerMove Tbl and offset meta structs

cd(ephysTbl_Dir)

if useOffset && pre_offset_ms > 0 && UniformEpochs && epochDur_ms > 0
    load(sprintf('filtProc_All_LFPsPerMove_pre%ims_post%ims.mat', pre_offset_ms, epochDur_ms), 'All_LFPsPerMove_Tbl_filt');
else
    load('filtProc_All_LFPsPerMove_N0offset.mat', 'All_LFPsPerMove_Tbl_filt');
end

% meta structs (ms, sec, ttl_samp, lfp_samp):
load('Offset_meta_struct.mat', 'meta_Offset')       % pre-trial offset from start of rep
load('UniformEpoch_meta_struct', 'meta_epochDur')   % post-trial duration from start rep


%% ------- Start here (unless workspace is clear) after Order_of_Op__LFPKin_Processing --------

% config / inputs for processing and analyses scripts should be identical

% after finalizing analysis functions here, make into seperate functions &
% update scripts (with single config / input step)


%% Extract processed LFP columns

% Extract LFP column names
allVars_LFPtbl     = string(All_LFPsPerMove_Tbl_filt.Properties.VariableNames);
lfpCols_raw = allVars_LFPtbl(startsWith(allVars_LFPtbl,"LFP_E") & endsWith(allVars_LFPtbl,"_filt"));     % 1375 Hz, post spectrumInterpolation
lfpCols_proc= allVars_LFPtbl(startsWith(allVars_LFPtbl,"LFP_E") & endsWith(allVars_LFPtbl,"_proc500"));  % 500 Hz, fully preprocessed

if isempty(lfpCols_raw) || isempty(lfpCols_proc)
    warning('Expected LFP_E*_filt and/or LFP_E*_proc500 not found.');
end


%% Test Plot - random LFP_E chan and row

% Pick a processed LFP column by E_channel name
lfp_colProc_500 = lfpCols_proc(1);                % e.g., "LFP_E1_proc500"

% Input - define single trial (row) to vizualize
test_rowIDX = 27;    % adjust
proc_LFP_vec = All_LFPsPerMove_Tbl_filt.(char(lfp_colProc_500)){test_rowIDX};

% Ensure frequency and time params are loaded
fs_downsamp = 500;              % downsampled rate = 500 Hz
t_step_proc = 1/fs_downsamp;    % 0.0020 sec
ts_LFP_proc = 0:t_step_proc:(length(proc_LFP_vec)-1)/fs_downsamp; % time vec

% Plot test trial
figure;
plot(ts_LFP_proc, proc_LFP_vec)
xlabel('Time (s)');
ylabel('LFP amplitude (µV)');
title('Processed LFP trial segment');
grid on;

% compute power spectral density (PSD) of LFP segment(s) - pspectrum function
[Pxx,Fxx] = pspectrum(proc_LFP_vec, fs_downsamp,'FrequencyLimits',[0 50],'FrequencyResolution', 3);

figure;
pspectrum(proc_LFP_vec, fs_downsamp, "power",'FrequencyLimits',[0 50],'FrequencyResolution', 3) % 2 or 3

% normalize using the common logarithm via decibel conversion - pow2db function
Pxx_db = pow2db(Pxx);

% Plot PSD for visualization
figure;
%plot(freq, 10*log10(power));
plot(Fxx, Pxx_db);
xlim([0 50])
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density of Processed LFP');

%% Test Plot per STN depth & trial number

% Input - define STN depth & move_trial # of trial to test:
STN_depth = 'dorsal';   % adjust
move_trial_num = '3';   % adjust
rep_idx = 5;

switch STN_depth
    case 'dorsal'
        LFP_test_depth = 't'; % top/dorsal STN
        trialID = [LFP_test_depth, move_trial_num];
        LFP_t_r1 = find(contains(All_LFPsPerMove_Tbl_filt.move_trial_ID, trialID), 1, 'first'); % Extract the first row index where this is true
        LFP_t_rn = find(contains(All_LFPsPerMove_Tbl_filt.move_trial_ID, trialID), 1, 'last'); % Extract the last row index where this is true
        LFP_t = All_LFPsPerMove_Tbl_filt.(char(lfp_colProc_500))(LFP_t_r1:LFP_t_rn); % Extract LFP segments for trial t3
        LFP_test_row = LFP_t{rep_idx};
    case 'central'
        LFP_test_depth = 'c'; % central STN
        trialID = [LFP_test_depth, move_trial_num];
        LFP_c_r1 = find(contains(All_LFPsPerMove_Tbl_filt.move_trial_ID, trialID), 1, 'first'); % Extract the first row index where this is true
        LFP_c_rn = find(contains(All_LFPsPerMove_Tbl_filt.move_trial_ID, trialID), 1, 'last'); % Extract the last row index where this is true
        LFP_c = All_LFPsPerMove_Tbl_filt.(char(lfp_colProc_500))(LFP_c_r1:LFP_c_rn); % Extract LFP segments for trial c2
        LFP_test_row = LFP_c{rep_idx};
    case 'ventral'
        LFP_test_depth = 'b'; % bottom/ventral STN
        trialID = [LFP_test_depth, move_trial_num];
        LFP_b_r1 = find(contains(All_LFPsPerMove_Tbl_filt.move_trial_ID, trialID), 1, 'first'); % Extract the first row index where this is true
        LFP_b_rn = find(contains(All_LFPsPerMove_Tbl_filt.move_trial_ID, trialID), 1, 'last'); % Extract the last row index where this is true
        LFP_b = All_LFPsPerMove_Tbl_filt.(char(lfp_colProc_500))(LFP_b_r1:LFP_b_rn); % Extract LFP segments for trial b2
        LFP_test_row = LFP_b{rep_idx};
end

% Display test trial ID prefix;
disp(trialID);

% Plot test trial
figure;
plot(ts_LFP_proc, LFP_test_row)
xlabel('Time (s)');
ylabel('LFP amplitude (µV)');
title('Processed LFP trial segment', ...
    [STN_depth, ' STN, trial #', move_trial_num, ' rep #', num2str(rep_idx)]);
grid on;

% compute power spectral density (PSD) of LFP segment(s) - pspectrum function
[Pxx,Fxx] = pspectrum(LFP_test_row, fs_downsamp,'FrequencyLimits',[0 50],'FrequencyResolution', 3);

figure;
pspectrum(LFP_test_row,fs_downsamp,"power",'FrequencyLimits',[0 50],'FrequencyResolution', 3) % 2 or 3

% normalize using the common logarithm via decibel conversion - pow2db function
Pxx_db = pow2db(Pxx);

% Plot PSD for visualization
figure;
%plot(freq, 10*log10(power));
plot(Fxx, Pxx_db);
xlim([0 50])
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density of Processed LFP', ...
    [STN_depth, ' STN, trial #', move_trial_num, ' rep #', num2str(rep_idx)]);


%% Run plot_LFP_byDepthTrial function for a test trial rep at each depth

% adjust function to create subplot:

dorsalSTN_rowIDX = plot_LFP_byDepthTrial(All_LFPsPerMove_Tbl_filt, 'dorsal', ...
                   move_trial_num, rep_idx, fs_downsamp);
centralSTN_rowIDX = plot_LFP_byDepthTrial(All_LFPsPerMove_Tbl_filt, 'central', ... 
                    move_trial_num, rep_idx, fs_downsamp);
ventralSTN_rowIDX = plot_LFP_byDepthTrial(All_LFPsPerMove_Tbl_filt, 'ventral', ... 
                    move_trial_num, rep_idx, fs_downsamp);


%% Beta Bandpass Filtering + log normalizing

lfp_data_test = proc_LFP_vec; 

% extract beta band using a 4th order IIR bandpass filter between 13-30 Hz 
% - designfilt function
beta_bandpass_filt = designfilt('bandpassiir','FilterOrder', 4, ...
    'HalfPowerFrequency1',13,'HalfPowerFrequency2',30, ...
    'SampleRate', fs_downsamp);

% pass band-filtered LFP signal through a zero-phase digital filter 
% - filtfilt function (to minimize phase distortion and transients)
beta_zerophase_filt = filtfilt(beta_bandpass_filt, lfp_data_test);


% compute power spectral density (PSD) of LFP segment(s) - pspectrum function
[beta_Pxx,beta_Fxx] = pspectrum(beta_zerophase_filt, fs_downsamp,'FrequencyLimits',[0 50],'FrequencyResolution',3); % 2 or 3

% normalize using the common logarithm via decibel conversion - pow2db function
betaPxxP = pow2db(beta_Pxx);

% plot normalized PSD of beta band
figure;
%plot(freq, 10*log10(power));
plot(beta_Fxx, beta_Pxx);
xlim([0 40])  % Limit x-axis to 50 Hz for better visualization
ylim([0 60])
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('PSD of LFP beta band (13-30 Hz)');


%% Hilbert transorm of beta-filtered signal

% compute instantaneous phase and frequency of the LFP beta band via Hilbert transform - hilbert function
hilbert_transformed = hilbert(beta_zerophase_filt);
inst_phase = angle(hilbert_transformed);
inst_freq = diff(unwrap(inst_phase))/(2*pi*t_step_proc);

LFP_beta_amp     = abs(hilbert(beta_zerophase_filt));    % inst. amplitude
LFP_beta_power   = LFP_beta_amp.^2;                      % inst. power
LFP_beta_zamp   = (LFP_beta_amp - mean(LFP_beta_amp,'omitnan')) ./ std(LFP_beta_amp,0,'omitnan'); % z-amp


%% Band isolation & Hilbert transformation (append column metrics)

% theta band
% Low beta bandpass (13–20 Hz)
% High beta bandpass (21-35 Hz)
% gamma band
% HFO

% Uses _proc500 columns as input (Fs = 500). Creates:
%   *_<band>_bp       : band-passed signal
%   *_<band>_amp      : instantaneous amplitude (abs(hilbert))
%   *_<band>_pow      : instantaneous power (amp.^2)
%   *_<band>_zamp     : z-scored amplitude (per-trial, optional but handy)


%  bands of interest
bands = struct( ...
    'theta',   [4 7], ...
    'L_beta',  [13 20], ...
    'H_beta',  [21 35], ...
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
    'SampleRate', fs_downsamp);

% Pre-build filters to avoid re-alloc
bp_filters = struct();
fieldNames = fieldnames(bands);
for k = 1:numel(fieldNames)
    curBand = bands.(fieldNames{k});
    % clamp to Nyquist
    curBand(2) = min(curBand(2), fs_downsamp/2 - 1);
    if curBand(1) >= curBand(2)-1
        error('Band %s limits invalid after clamping.', fieldNames{k});
    end
    bp_filters.(fieldNames{k}) = BPfilt(curBand(1), curBand(2));
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

    for row_i = 1:nRows
        LFPvec = All_LFPsPerMove_Tbl_filt.(inCol){row_i};
        if isempty(LFPvec) || ~isvector(LFPvec) || ~isnumeric(LFPvec), continue; end
        LFPvec = double(LFPvec(:));
        if ~any(isfinite(LFPvec)), continue; end

        for k = 1:numel(fieldNames)
            band = fieldNames{k};
            try
                LFP_vec_bp = filtfilt(bp_filters.(band), LFPvec);           % bandpassed
            catch ME
                warning('Row %d, %s (%s): bandpass failed: %s', ...
                    row_i, inCol, band, ME.message);
                continue
            end
            LFP_vec_amp   = abs(hilbert(LFP_vec_bp));                       % inst. amplitude
            LFP_vec_p   = LFP_vec_amp.^2;                                   % inst. power
            LFP_vec_zamp   = (LFP_vec_amp - mean(LFP_vec_amp,'omitnan')) ./ std(LFP_vec_amp,0,'omitnan'); % z-amp

            All_LFPsPerMove_Tbl_filt.(base + "_" + band + "_bp"){row_i}   = LFP_vec_bp;
            All_LFPsPerMove_Tbl_filt.(base + "_" + band + "_amp"){row_i}  = LFP_vec_amp;
            All_LFPsPerMove_Tbl_filt.(base + "_" + band + "_pow"){row_i}  = LFP_vec_p;
            All_LFPsPerMove_Tbl_filt.(base + "_" + band + "_zamp"){row_i} = LFP_vec_zamp;
        end
    end
end

fprintf('[DONE] Added bandpassed, amplitude, power, and z-amplitude columns.\n');


%% Save updated All_LFPsPerMove_Tbl_filt with computed characteristics per band

% Resave All_LFPsPerMove_Tbl_filt
cd(ephysTbl_Dir)

if useOffset && pre_offset_ms > 0 && UniformEpochs && epochDur_ms > 0
    bp_outName = sprintf('bp_All_LFPsPerMove_pre%ims_post%ims.mat', pre_offset_ms, epochDur_ms);
else
    bp_outName = 'bp_All_LFPsPerMove_N0offset.mat';
end

save(bp_outName, "All_LFPsPerMove_Tbl_filt");
fprintf('[SAVED] Bandpassed LFPsPerMove table written to %s\n', fullfile(ephysTbl_Dir, bp_outName));


%% Auto-split STN depths

dorsalSTN_tbl  = All_LFPsPerMove_Tbl_filt(contains(All_LFPsPerMove_Tbl_filt.move_trial_ID,'t'),:);
centralSTN_tbl = All_LFPsPerMove_Tbl_filt(contains(All_LFPsPerMove_Tbl_filt.move_trial_ID,'c'),:);
ventralSTN_tbl = All_LFPsPerMove_Tbl_filt(contains(All_LFPsPerMove_Tbl_filt.move_trial_ID,'b'),:);

LFPs_perSTNdepth_Tables = struct('dorsal', dorsalSTN_tbl, 'central', centralSTN_tbl, 'ventral', ventralSTN_tbl);


%% Plot PSDs for different movement contexts per STN depth

moveType_ids = unique(All_LFPsPerMove_Tbl_filt.MoveType);
% Rest, Hand OC, Arm Pro/Sup, Arm Exten/Flex

% Preallocate field names for moveTbls


% Loop though moveTypes and fill moveTbls per moveType per depth
for moveT_i = 1:height(moveType_ids)
    moveType_str = moveType_ids{moveT_i};
    dorsal_moveTbl{moveT_i} = dorsalSTN_tbl(contains(dorsalSTN_tbl.MoveType, moveType_str), :);
    central_moveTbl{moveT_i} = centralSTN_tbl(contains(centralSTN_tbl.MoveType, moveType_str), :);
    ventral_moveTble{moveT_i} = ventralSTN_tbl(contains(ventralSTN_tbl.MoveType, moveType_str), :);

    % Extract LFP segments per move type per depth


end




%% Compute LFP power for each moveType per depth

% Hilbert transformation spec. freq. (inst. power)
% theta band
% low-beta band
% high-beta band
% gamma band
% HFO

% inst. power at specific time/context (move index)

%% Run run_MovementFeatureAnalysis_IO_v2 or load kinTbl & kinSummaryTbl

% cd(KinematicsDir)

% fprintf('[INFO] Loading movement data from: %s\n', MoveDataDir);
% [kinTbl, kinSummaryTbl] = run_MovementFeatureAnalysis_IO_v2(IO_DataDir, MoveDataDir, Move_CaseID);


%% FOOOF/spec param on different contexts

% * The aperiodic exponent of subthalamic field potentials reflects excitation/inhibition balance in Parkinsonism: https://pubmed.ncbi.nlm.nih.gov/36810199/; https://elifesciences.org/articles/82467#data
% * Parameterizing neural power spectra into periodic and aperiodic components: https://www.nature.com/articles/s41593-020-00744-x
% * Electrophysiological Frequency Band Ratio Measures Conflate Periodic and Aperiodic Neural Activity: https://www.eneuro.org/content/7/6/ENEURO.0192-20.2020
% * Electrophysiological Frequency Band Ratio Measures Conflate Periodic and Aperiodic Neural Activity: https://onlinelibrary.wiley.com/doi/abs/10.1111/ejn.15361
% * Time-resolved parameterization of aperiodic and periodic brain activity: Time-resolved parameterization of aperiodic and periodic brain activity | eLife
% * Subthalamic nucleus encoding steers adaptive therapies for gait in Parkinson's disease: https://www.medrxiv.org/content/10.1101/2025.08.20.25333478v1


%% Spike-field coherence

% at specific movement timepoints

%% LFP bursting analysis / time-frequency domain


%% Continuous Wavelet Transformation
% look into math behind cwt fun

% bursting analysis, based on Torrecillos et al., 2018 (test for LFPs in 1 trial / STN depth)
% time-frq transform via complex Morlet wavelets (f0/σf = 7, f0 = 1-45 Hz, steps of 0.25 Hz).

% define 'lfp_data' as a LFP data matrix with dimensions [samples x trials]
lfp_data_test = proc_LFP_vec; 
lfp_data_test_beta = beta_zerophase_filt;
time = (1:length(lfp_data_test)) / fs_downsamp;

% time-domain bin size determined by 'VoicesPerOctave' (higher VPO = finer freq. res; lower time res?)
% continuous wavelet transform
[cwt_power, cwt_freq] = cwt(lfp_data_test, fs_downsamp ,'FrequencyLimits',[1 125],'VoicesPerOctave',14); % try 14, mults of 7
figure;
cwt(lfp_data_test,fs_downsamp,'FrequencyLimits',[1 125],'VoicesPerOctave',14); % try 14, mults of 7
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([5 50])
title('LFP in T-F domain during motor trial rep', 'Continuous Wavelet Transform');

% Normalize power for each frequency band
mean_power = mean(cwt_power, 2);
std_power = std(cwt_power, 0, 2);
normalized_power = (cwt_power - mean_power) ./ std_power;


%% Beta Peak Selection

% beta_range = find(f0 >= 13 & f0 <= 30);
beta_range = find(cwt_freq >= 13 & cwt_freq <= 30);
[~, peak_idx] = max(mean(normalized_power(beta_range, :), 2));
beta_peak_frequency = cwt_freq(beta_range(peak_idx));

% Compute beta power time courses for each trial using a 6-Hz-wide frequency band centered on the beta peak frequency
beta_power_time_courses = mean(normalized_power(beta_range(peak_idx) + (-3:3), :), 1);

% Plot the beta power time course
figure('Renderer', 'opengl');
plot(time, beta_power_time_courses);
xlabel('Time (s)');
ylabel('Normalized Beta Power');
title('LFP Beta Power (untrimmed)')


%% Trim beta power lfp and time data based on known/computed offset

% 1) find index of artifact, chop signal before that artifact,  
% run following burst detection code 

% Define a threshold as a percentage of the maximum beta power
power_threshold = 0.1 * max(beta_power_time_courses);

% Find the index where the power first exceeds the threshold from the end
drop_off_index = find(beta_power_time_courses(end:-1:1) > power_threshold, 1, 'first');
if isempty(drop_off_index)
    drop_off_index = length(beta_power_time_courses);  % Use the full length if no drop-off is found
else
    drop_off_index = length(beta_power_time_courses) - drop_off_index + 1;  % Convert to index from the start
end

% Trim the beta power time courses and time data based on known/computed offset
beta_power_time_courses = beta_power_time_courses(1:drop_off_index);
time_trimmed = time(1:drop_off_index);

%% continuous wavelet transform on trimmed LFP

lfp_data_L_trimmed = lfp_data_test(1:drop_off_index);
figure;
cwt(lfp_data_L_trimmed,fs_downsamp,'FrequencyLimits',[1 125],'VoicesPerOctave',14); % try 14, mults of 7'
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([5 50])
title('LFP (trimmed) during motor trial rep', 'Continuous Wavelet Transform');


%% Beta Burst Detection and Plotting for trimmed LFP

% Compute mean beta power
mean_beta_power = mean(beta_power_time_courses); % complex double
mean_beta_power_real = mean(real(beta_power_time_courses)); % real part
mean_beta_power_imag = mean(imag(beta_power_time_courses)); % imaginary part

mean_beta_power_magnitude = mean(abs(beta_power_time_courses)); % mean of magnitudes

% Set threshold at the 75th percentile of the mean beta power
threshold = prctile(mean_beta_power_magnitude, 75);

% Define beta bursts as time points exceeding the threshold for more than 3 oscillatory cycles (3 oscillatory cycles = ~100ms)
min_duration = round(0.1 * fs_downsamp); % Minimum duration of a burst in samples (100 ms)
% min_duration = round(0.066 * fs);  % Reduce minimum duration to 2 osc. cycles(66 ms)

above_threshold = beta_power_time_courses > threshold;
above_threshold = above_threshold(:);  % Ensure it's a column vector
bursts = zeros(size(above_threshold));  % Initialize bursts array
[start_indices, end_indices] = findConsecutiveIndices_fun(above_threshold);  % Find start and end indices of consecutive regions above threshold
for k = 1:length(start_indices)
    if (end_indices(k) - start_indices(k) + 1) >= min_duration
        bursts(start_indices(k):end_indices(k)) = 1;  % Mark burst regions
    end
end

% Plot the beta power time course
figure('Renderer', 'opengl');
plot(time_trimmed, beta_power_time_courses);
xlabel('Time (s)');
ylabel('Normalized Beta Power');
title('LFP Beta Power')


%% Plot the beta power time course and beta burst overlays using xregion (test with L_streamOfInt_s1)

% https://www.mathworks.com/help/matlab/ref/xregion.html

figure('Renderer', 'opengl');
p = plot(time_trimmed, beta_power_time_courses); % store handle to the beta power plot
hold on;

% Overlay the detected bursts on the same plot using xregion
for i = 1:length(start_indices)
    burst_start_time = time_trimmed(start_indices(i));
    burst_end_time = time_trimmed(end_indices(i));
    xregion(burst_start_time, burst_end_time, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8]);
end

% Disable automatic legend updates to prevent xregion from adding entries
legend('auto update','off');

% Create a dummy patch handle 'h' for the detected bursts to add to the legend
h = patch(NaN, NaN, [0.8, 0.8, 0.8], 'FaceAlpha', 0.3);

% Add the beta power plot handle 'p' to the legend
legend([p, h], 'Beta Power', 'Detected Bursts');

% Plot
xlabel('Time (s)');
ylabel('Normalized Beta Power');
title('Beta Bursting Dynamics');
hold off;

