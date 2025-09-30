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


%% Config - Define pre-trial offset duration for useOffset_LFP function

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

%% Run useOffset_LFP helper function

[offset_LFP_samples, meta_Offset] = useOffset_LFP(TTL_fs, AO_LFP_fs, pre_offset_ms, useOffset);

% If useOffset == true or pre_offset_ms>0, useOffset_LFP function returns the #
% of samples to pre-pad in LFP domain based on a set pre-trial offset time
% (and a meta struct for reference).
% If useOffset == false or pre_offset_ms<=0, useOffset_LFP function returns 0.

%% Config - Define epoch duration for UniformEpochs_LFP function

% Define # of milliseconds to include from start of a move rep
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


%% Run UniformEpochs_LFP function

[epochDur_LFP_samples, meta_epochDur] = UniformEpochs_LFP(TTL_fs, AO_LFP_fs, epochDur_ms, UniformEpochs);


%% Config - Case-specific Inputs

CaseDate = '07_06_2023_bilateral'; % Adjust as needed

% '03_23_2023'; % NER 2025
% '04_05_2023'; % NER 2025
% '05_18_2023_b_bilateral'; % NER 2025
% '05_31_2023';
% '06_08_2023_bilateral'; % NER 2025
% '08_23_2023'; % NANS 2026
% '07_06_2023_bilateral';

MoveDir_CaseID = 'IO_07_06_2023_LSTN'; % Adjust as needed

% 'IO_03_23_2023_LSTN'; % NER 2025
% 'IO_04_05_2023_RSTN'; % NER 2025
% 'IO_05_18_2023_b_LSTN'; % NER 2025
% 'IO_05_31_2023_LSTN';
% 'IO_06_08_2023_LSTN'; % NER 2025
% 'IO_08_23_2023_RSTN'; % NANS 2026
% 'IO_07_06_2023_LSTN';


%% Config - Input and Output Data Dirs

% Case-specific Input dirs
Ephys_CaseDir = [IO_DataDir, filesep, CaseDate];
ProcDataDir = [Ephys_CaseDir, filesep, 'Processed Electrophysiology'];      % directory for processed ephys data and spike clusters

MoveDataDir = [IO_DataDir, filesep, 'Processed DLC'];
Move_CaseDir = [MoveDataDir, filesep, MoveDir_CaseID];                      % directory for processed DLC data and Movement Indices

% Case-specific Output dir
ephysTbl_Dir = [Ephys_CaseDir, filesep, 'DLC_Ephys'];                       % directory where all ephys per move-rep tables are located
if ~exist(ephysTbl_Dir,'dir')
    mkdir(ephysTbl_Dir);
end

%% Define new output directory:

% Ephys_Kinematics
LFP_Kin_Dir = [IO_DataDir, filesep, 'Ephys_Kinematics', filesep, 'LFP_Kinematic_Analyses'];
Case_LFPKin_Dir = [LFP_Kin_Dir, filesep, CaseDate];
if ~exist(Case_LFPKin_Dir,'dir')
    mkdir(Case_LFPKin_Dir);
end

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
    Case_LFPKin_Dir = fullfile(Case_LFPKin_Dir, CaseDate_hem);
    if ~exist(Case_LFPKin_Dir,'dir')
        mkdir(Case_LFPKin_Dir);
    end
    fprintf('[INFO] Hemisphere-specific output directory set: %s\n', Case_LFPKin_Dir);

else
    fprintf('[INFO] Using base ephys input directory: %s\n', ProcDataDir);
    fprintf('[INFO] Using base output directory: %s\n', Case_LFPKin_Dir);
end



%% ------ Functions ------

%% Run align_LFPsPerMove_TTL function

% All_LFPsPerMove_Tbl = align_LFPsPerMove_TTL(Subject_AO, ProcDataDir, Move_CaseDir, ephysTbl_Dir, TTL_fs, AO_LFP_fs, pre_offset_ms, useOffset);

%% Run align_LFPsPerMove_uniformEpochDur function

% Add 500-1000 ms from start of each MoveRep

% All_LFPsPerMove_Tbl_uniformEpochs = align_LFPsPerMove_uniformEpochDur(Subject_AO, ProcDataDir, Move_CaseDir, ephysTbl_Dir, TTL_fs, AO_LFP_fs, pre_offset_ms, useOffset, UniformEpochs, epochDur_ms);
All_LFPsPerMove_Tbl_uniformEpochs = align_LFPsPerMove_uniformEpochDur(Subject_AO, ProcDataDir, Move_CaseDir, Case_LFPKin_Dir, TTL_fs, AO_LFP_fs, pre_offset_ms, useOffset, UniformEpochs, epochDur_ms);

% meta structs (ms, sec, ttl_samp, lfp_samp):
% meta_Offset    :  pre-trial offset from start of rep
% meta_epochDur  :  post-trial duration from start rep


%% Inputs for spectrumInterpolation (Miguel's function :))

% This function interpolates around the frequency of interest (Fl) and
% replaces its and some neighbors using a constant value.

% function inputs
% data = All_LFPsPerMove_Tbl_uniformEpochs.LFP; % data: column vector of the data that needs to be filtered
Fs = AO_LFP_fs; % Fs: Sampling Frequency (in Hz) of the data, % loop on harmonics of 60 Hz
Fl = 60; % Fl: Line frequency (in Hz), center of our interpolation, loop on harmonics of 60 Hz (notch & comb)
neighborsToSample = 4; % Hz, 4 or 5
neighborsToReplace = 2; % Hz, 1 or 2

% neighborsToSample: This parameter is in Hz, and tells this function how
% large of a window (in Hz) to use when picking the constant value to
% replace frequency content with.

% neighborsToReplace: This parameter is also in Hz, and tells this function
% which neighbors need to be replaced with the constant value that was
% determined by neighborsToSample. Generally, neighborsToReplace <
% neighborsToSample in order to get a better spectral estimate.


%% Run wrapper_applySpectrumInterpolation_LFPs function

% pick which columns to filter (LFP channels) ---
LFPsPerMoveTbl_vars = string(All_LFPsPerMove_Tbl_uniformEpochs.Properties.VariableNames);
lfpCols = LFPsPerMoveTbl_vars(startsWith(LFPsPerMoveTbl_vars,"LFP_E"));     % e.g., LFP_E1, LFP_E2, ...

% create new columns with "_filt" after running each through spectrumInterpolation
makeNewCols = true;

% run wrapper function for spectrumInterpolation and output updated All_LFPsPerMove_Tbl
All_LFPsPerMove_Tbl_filt = wrapper_applySpectrumInterpolation_LFPs(All_LFPsPerMove_Tbl_uniformEpochs, ...
    lfpCols, AO_LFP_fs, Fl, neighborsToSample, neighborsToReplace, ...
    makeNewCols);


%% save All_LFPsPerMove_Tbl with LFP_E column vectors filtered by spectrumInterpolation

cd(Case_LFPKin_Dir)

if useOffset && pre_offset_ms > 0 && UniformEpochs && epochDur_ms > 0
    specI_filt_outName = sprintf('specI_All_LFPsPerMove_pre%ims_post%ims.mat', pre_offset_ms, epochDur_ms);
else
    specI_filt_outName = 'specI_All_LFPsPerMove_N0offset.mat';
end

save(specI_filt_outName, "All_LFPsPerMove_Tbl_filt");


%% filters (notes)

% 1) high-pass at 0.5 or 1 Hz
% low freq. noise + drift removal
% use designfilt + filtfilt

% 1) low-pass at 250 Hz
% anti-alias + keep LFP band
% use designfilt + filtfilt

% 3) downsample / resample 1375 to 500 Hz

% Question (in initial notes)
% 1-3 kHz (?)


%% Filtering (after spectrumInterpolation)

% Input:  All_LFPsPerMove_Tbl_filt (with LFP_E*_filt columns, Fs = 1375 Hz)
% Output: Preprocessed vectors (high-pass, low-pass, downsampled to 500 Hz)

fs_AO  = AO_LFP_fs;   % original LFP sampling rate = 1375 Hz
fs_downsamp = 500;    % downsampled rate = 500 Hz

fs_AO_Nyquist = fs_AO./2;       % 687.5 Hz
fs_DS_Nyquist = fs_downsamp./2; % 250 Hz

% Resampling Ratio
resamp_Ratio = fs_downsamp./fs_AO; % (new_fs/orig_fs) = 500/1375  =  0.3636

% Get the rational approximation
[up_factor, down_factor] = rat(resamp_Ratio); % 4/11
disp(['Resampling Ratio, Fractional Representation: ', num2str(up_factor), '/', num2str(down_factor)]);


% Design filters (once)
hpFilt = designfilt('highpassiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency', 1, ...      % highpass at ~1 Hz
    'SampleRate', fs_AO);

lpFilt = designfilt('lowpassiir', ...
    'FilterOrder', 8, ...
    'HalfPowerFrequency', 250, ...    % lowpass at 250 Hz (or ~200–230 Hz)
    'SampleRate', fs_AO);

% Pick which columns to filter (post-spectrumInterpolation: '_filt')
lfpCols_filt = string(All_LFPsPerMove_Tbl_filt.Properties.VariableNames);
lfpCols_filt = lfpCols_filt(startsWith(lfpCols_filt,"LFP_E") & endsWith(lfpCols_filt,"_filt"));

% Loop through each LFP column and row
for col_i = 1:numel(lfpCols_filt)
    colName = lfpCols_filt(col_i);
    proc_col_outName = replace(colName,"_filt","_proc500");  % output column name

    % Preallocate
    All_LFPsPerMove_Tbl_filt.(proc_col_outName) = cell(height(All_LFPsPerMove_Tbl_filt),1);

    for row_i = 1:height(All_LFPsPerMove_Tbl_filt)
        LFP_vec = All_LFPsPerMove_Tbl_filt.(colName){row_i};
        if isempty(LFP_vec) || ~isnumeric(LFP_vec) || ~isvector(LFP_vec)
            All_LFPsPerMove_Tbl_filt.(proc_col_outName){row_i} = LFP_vec;
            continue
        end

        % Filtering steps
        LFP_vec = double(LFP_vec(:));                                       % ensure column, double
        LFP_vec_hp = filtfilt(hpFilt, LFP_vec);                             % high-pass, lowF drift removal
        LFP_vec_lp = filtfilt(lpFilt, LFP_vec_hp);                          % low-pass anti-alias
        LFP_vec_downsamp = resample(LFP_vec_lp, up_factor, down_factor);    % downsample 1375 → 500 Hz (resamp ratio = 4/11)

        % Store result
        All_LFPsPerMove_Tbl_filt.(proc_col_outName){row_i} = LFP_vec_downsamp;
    end
end


%% Visualization (Quality check)

testRow = 6;   % change index to preview other trials
raw_LFP_vec   = All_LFPsPerMove_Tbl_filt.(lfpCols_filt(1)){testRow};
proc_LFP_vec  = All_LFPsPerMove_Tbl_filt.(replace(lfpCols_filt(1),"_filt","_proc500")){testRow};

% fs_AO = AO_LFP_fs; % 1375 Hz
t_step_AO = 1/fs_AO;
ts_LFP_AO = 0:t_step_AO:(length(raw_LFP_vec)-1)/fs_AO;

% fs_downsamp = 500; % 500 Hz
t_step_proc = 1/fs_downsamp;
ts_LFP_proc = 0:t_step_proc:(length(proc_LFP_vec)-1)/fs_downsamp;

% visualize
figure;
subplot(2,1,1);
plot(ts_LFP_AO, raw_LFP_vec); title('LFP (post spectrumInterpolation, 1375 Hz)');
xlabel('Time (s)');
ylabel('Preprocessed LFP'); % (uV or unscaled)
subplot(2,1,2);
plot(ts_LFP_proc, proc_LFP_vec); title('Filtered & Processed LFP (post hpfilt + lpfilt + downsampling, 500 Hz)');
xlabel('Time (s)');
ylabel('Processed LFP'); % (uV or unscaled)


%% Questions

% 1) Why/how are the y-axis magnitudes this different?
% is this ok/expected or is there something wrong in my filtering methods?

% 2) Are the LFP units in microVolts (uV) and being plotted as such,
% or are the y-axis units different/unscaled?

% 3) Do I need to adjust the lowpass cutoff from 250 Hz to ~200–230 Hz
% to be below the Nyquist of the downsampled fs (500 Hz)?


%% Plot Power Spectral Density of Raw LFP (post spectrumInterpolation)

% compute power spectral density (PSD) of LFP segment(s) - pspectrum function
[rawPxx,rawFxx] = pspectrum(raw_LFP_vec,fs_AO,'FrequencyLimits',[0 90],'FrequencyResolution', 3);

figure;
pspectrum(raw_LFP_vec,fs_AO,"power",'FrequencyLimits',[0 90],'FrequencyResolution', 3) % 2 or 3

% normalize using the common logarithm via decibel conversion - pow2db function
rawPxx_db = pow2db(rawPxx);

% plot PSD for visualization
figure;
%plot(freq, 10*log10(power));
plot(rawFxx, rawPxx_db);
xlim([0 80])
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('PSD of Preprocessed LFP (post spectrumInterpolation)');


%% Plot Power Spectral Density of Processed LFP (post hpfilt + lpfilt + downsampling, 500 Hz)

% compute power spectral density (PSD) of LFP segment(s) - pspectrum function
[procPxx,procFxx] = pspectrum(proc_LFP_vec,fs_downsamp,'FrequencyLimits',[0 100],'FrequencyResolution', 3);

figure;
pspectrum(proc_LFP_vec,fs_downsamp,"power",'FrequencyLimits',[0 100],'FrequencyResolution', 3) % 2 or 3

% normalize using the common logarithm via decibel conversion - pow2db function
procPxx_db = pow2db(procPxx);

% plot PSD for visualization
figure;
%plot(freq, 10*log10(power));
plot(procFxx, procPxx_db);
xlim([0 80])
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('PSD of Processed LFP (post hp/lp filtering and downsampling)');


%% Save fully filtered and processed LFP outputs (HP+LP+DS to 500 Hz)

% Resave All_LFPsPerMove_Tbl_filt
cd(Case_LFPKin_Dir)

if useOffset && pre_offset_ms > 0 && UniformEpochs && epochDur_ms > 0
    filtProc_outName = sprintf('filtProc_All_LFPsPerMove_pre%ims_post%ims.mat', pre_offset_ms, epochDur_ms);
else
    filtProc_outName = 'filtProc_All_LFPsPerMove_N0offset.mat';
end

save(filtProc_outName, "All_LFPsPerMove_Tbl_filt");
fprintf('[SAVED] Processed table written to %s\n', fullfile(Case_LFPKin_Dir, filtProc_outName));


%% Go to Order_of_Op__LFPKin_Analyses script