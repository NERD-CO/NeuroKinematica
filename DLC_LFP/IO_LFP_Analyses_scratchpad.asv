%% IO LFP analyses

% specify directory where case-specific data files are located
curPCname = getenv('COMPUTERNAME');

switch curPCname
    case 'DESKTOP-I5CPDO7'  % PC_1
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'DSKTP-JTLAB-EMR'  % Lab Desktop
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'NSG-M-H8J3X34'    % PC_2
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
end

cd(IO_DataDir)
Subject_AO = readtable('Subject_AO.xlsx');

%% Hardcode case-specific ephys data directories

% isolate a specific CaseDate / studyID (StudyNum in Subject_AO csv)
CaseDate = '03_09_2023'; % studyID = 1
% CaseDate = '03_23_2023'; % studyID = 2

case_ID = CaseDate;

% define case-specific data directory locations
Case_DataDir = [IO_DataDir, filesep, CaseDate];
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology'];       % directory where processed ephys .mat data should be saved
ephysTbl_Dir = [Case_DataDir, filesep, 'DLC_MER'];                          % directory where all ephys per move-rep tables are located


%% Define case-specific directories for kinematic data and movement indices

% define kinematic data directory
MoveDataDir = [IO_DataDir, filesep, 'Kinematic Analyses'];

% specify case ID
Move_CaseID = 'IO_03_09_2023_RSTN'; % studyID = 1
% Move_CaseID = 'IO_03_23_2023_LSTN'; % studyID = 2

% isolate case-specific kinematic data directory
Move_CaseDir = [MoveDataDir, filesep, Move_CaseID];

% data subfolders:
Move_CaseMats = [Move_CaseDir, filesep, 'mat folder'];
Move_CaseVideos = [Move_CaseDir, filesep, 'video folder'];

% for kinematic analyses
cd(Move_CaseMats)
moveMat = dir('*.mat');
moveMat_names = {moveMat.name};

% Generate list of Motor Index CSVs
cd(Move_CaseVideos)
moveCSV = dir('*.csv');
moveCSV_names = {moveCSV.name};

% Generate list of Motor Index CSVs (filters for CSVs that contain 'Move' string)
moveCSV_list = moveCSV_names(contains(moveCSV_names,'Move'));

%% define datastream sampling rates (Alpha Omega and Video sampling fs)

TTL_fs = 44000; % Hz
AO_spike_fs = 44000; % Hz
AO_LFP_fs = 1375; % Hz
DLC_fs = 100; % fps


%% define offset duration
offset_ms = 50; % milliseconds
offset_seconds = offset_ms / 1000; % seconds

% Calculate number of TTL samples
offset_TTLs = round(TTL_fs * offset_seconds); % ensure value is integer

% Calculate number TTL samples in AO_LFP sample domain by downsampling TTL_fs
offset_TTLs_LFP = round((offset_TTLs/TTL_fs)*AO_LFP_fs); % ensure value is integer

% for future function input: useOffset; when = 1, use offset; when = 0, don't

%% Load all processed LFPs %%% fix this block later

cd(ProcDataDir)

% list of filename
LFPmatfiles = dir('*.mat');
LFPmatnames = {LFPmatfiles.name};

for LFP_mat_name = 1:length(LFPmatnames)

    cd(ProcDataDir)
    load(LFPmatnames{LFP_mat_name},'ProcEphys')

    fileparts = split(LFPmatnames{LFP_mat_name},'_');
    ProcName = fileparts{2};

    ProcFile = LFPmatnames{contains(LFPmatnames, ProcName)};
    load(ProcFile, 'ProcEphys')

    % account for mult. electrode channels
    electrod_names = fieldnames(ProcEphys.LFP); % get num
    for e_names = 1:length(electrod_names)

        LFP_raw = ProcEphys.LFP.(electrod_names{e_names}).rawData; % dynamically index within a struct

        % Find row of ao_MAT_file that corresponds with trial
        SubjectAO_row = Subject_AO(contains(Subject_AO.ao_MAT_file,ProcName),:);

        switch SubjectAO_row.stn_loc{1}(1)
            case 'd'
                depthName = 't'
            case 'v'
                depthName = 'b'
            otherwise
                depthName = 'c'
        end

        motor_trial_ID = [depthName, num2str(SubjectAO_row.trialNum),'_', SubjectAO_row.depth{1}];

        % Generate list of Motor Index CSVs
        cd(Move_CaseVideos)
        moveCSV = dir('*.csv');
        moveCSV_names = {moveCSV.name};

        % Generate list of Motor Index CSVs (filters for CSVs that contain 'Move' string)
        moveCSV = moveCSV_names(contains(moveCSV_names,'Move'));

        if ~any(contains(moveCSV, motor_trial_ID))  % checking if logical is false
            continue
        end

        moveTbl_name = moveCSV{contains(moveCSV, motor_trial_ID)};
        moveTbl = readtable(moveTbl_name);

        electrode_LFP_name = [ProcName,'_', electrod_names{e_names}]; %temp var
        LFP_raw_depthName = [ProcName,'_', motor_trial_ID]; % temp var

    end

end


%% Load all LFPs per Movement Tbl %%% fix this block later

addpath 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\AO_Processing'

cd(ephysTbl_Dir)

% list of filename
Tblmatfiles = dir('*.mat');
Tblmatnames = {Tblmatfiles.name};

use_offset = 1;
for Tblmat_name = 1:length(Tblmatnames)
    fileparts = split(Tblmatnames{Tblmat_name},'_');
    if use_offset == 1
        Tbl_name = fileparts{2:3};
    else
        Tbl_name = fileparts{2};
    end
end

LFPsPerMove_Tbl = Tblmatnames{contains(Tblmatnames, 'All_LFPsPerMove_offset')};
load(LFPsPerMove_Tbl, 'All_LFPsPerMove_Tbl')


%% LFP frequency domain analysis 

% Define LFP segment(s) of interest:
depth_trial_IDs = unique(All_LFPsPerMove_Tbl.move_trial_ID);
for depth_trial_i = 1:length(depth_trial_IDs)
    depth_trial = All_LFPsPerMove_Tbl.move_trial_ID{depth_trial_i};

    %%% fill in later

end

% depth_trial = 'b1_d0p4'
LFP_depth_trial_b1_array = table2array(All_LFPsPerMove_Tbl(1:18,8));
LFP_depth_trial_b1 = horzcat(LFP_depth_trial_b1_array{:});
TTL_idx_b1 = All_LFPsPerMove_Tbl{1:18,5};

% depth_trial = 'b3_d0p4'
LFP_depth_trial_b3_array = table2array(All_LFPsPerMove_Tbl(37:54,8));
LFP_depth_trial_b3 = horzcat(LFP_depth_trial_b3_array{:});
TTL_idx_b3 = All_LFPsPerMove_Tbl{37:54,5};

% depth_trial = 'c2_d2p0'
LFP_depth_trial_c3_array = table2array(All_LFPsPerMove_Tbl(73:end,8));
LFP_depth_trial_c3 = horzcat(LFP_depth_trial_c3_array{:});
TTL_idx_c3 = All_LFPsPerMove_Tbl{73:end,5};

% define test trial
LFP_depth_test1 = LFP_raw(:, TTL_idx_c3(1):TTL_idx_c3(end)); % LFP_depth_trial_c3


%% Target signal (LFP) sampling rate (samples per sec)
fs = AO_LFP_fs; % 1375 Hz
t_step = 1/fs;
ts_LFP = 0:t_step:(length(LFP_depth_test1)-1)/fs;

% compute power spectral density (PSD) of LFP segment(s) - pspectrum function
[rPxx,rFxx] = pspectrum(LFP_depth_test1,fs,'FrequencyLimits',[0 400],'FrequencyResolution',3); 
figure;
pspectrum(LFP_depth_test1,fs,"spectrogram",'FrequencyLimits',[0 400],'FrequencyResolution',3) % 2 or 3

% normalize using the common logarithm via decibel conversion - pow2db function
rPxxP = pow2db(rPxx);

% plot PSD for visualization
figure;
%plot(freq, 10*log10(power));
plot(rFxx, rPxxP);
xlim([0 400]) 
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density of Unfiltered LFP');


%% LFP filtering
%%% fix this %%%

% remove 60 Hz line noise and harmonics via 2nd-order infinite impulse response (IIR) notch and comb filters

% Notch Filter centered at 60Hz
f_notch = 60;  % Notch frequency
Q = 35;   % Quality factor, q = ω0/bw where ω0 is the frequency to remove from the signal
[b_notch, a_notch] = iirnotch(f_notch/(fs/2), f_notch/(fs/2)/Q); % 2nd order infinite impulse response (IIR) notch filter
LFP_dt1_notchfilt = filtfilt(b_notch, a_notch, LFP_depth_test1); % Apply notch filter to remove 60Hz line noise

% Comb Filter centered at 60Hz and its harmonics (120Hz, 180Hz, 240Hz, 300Hz, 360Hz, etc.) 
bw = f_notch / Q;  % Bandwidth, bw = (fo/(fs/2))/q;
n_harmonics = floor((fs/2) / f_notch);  % Number of harmonics within Nyquist frequency
[b_comb, a_comb] = iircomb(floor(fs/f_notch), bw/(fs/2), 'notch');  % Design comb filter
LFP_dt1_combfilt = filtfilt(b_comb, a_comb, LFP_dt1_notchfilt); % Apply comb filter to remove harmonics


% attenuate low frequency components in LFP via 4th order Butterworth or IIR high-pass filter

% High-pass filter 
cutoff_freq = 5; % Hz, attenuate frequency components below cutoff frequency 
hp_filter = designfilt('highpassiir', 'FilterOrder', 4, ... 
    'HalfPowerFrequency', cutoff_freq, 'SampleRate', fs); % Design 4th order Butterworth or IIR hp-filter
LFP_dt1_hpfilt = filtfilt(hp_filter, LFP_dt1_combfilt); % Apply high-pass filter to notch/comb-filtered LFP signal

% Plot the filtered LFP signal
figure;
plot(ts_LFP, LFP_dt1_hpfilt);
xlabel('Time (s)');
ylabel('Amplitude');
title('High-Pass Filtered LFP Signal');


% compute power spectral density (PSD) of LFP segment(s) - pspectrum function
[fPxx,fFxx] = pspectrum(LFP_dt1_hpfilt,fs,'FrequencyLimits',[0 100],'FrequencyResolution',3); 
figure;
pspectrum(LFP_dt1_hpfilt,fs,"spectrogram",'FrequencyLimits',[0 100],'FrequencyResolution',3) % 2 or 3

% normalize using the common logarithm via decibel conversion - pow2db function
fPxxP = pow2db(fPxx);

% plot PSD for visualization
figure;
%plot(freq, 10*log10(power));
plot(fFxx, fPxxP);
xlim([0 100]) 
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density of Filtered LFP');


%% Beta Bandpass Filtering + log normalizing

% extract beta band using a 4th order IIR bandpass filter between 13-30 Hz - designfilt function
beta_bandpass_filt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',13,'HalfPowerFrequency2',30, ...
    'SampleRate',fs);

% pass beta band-filtered LFP signal through a zero-phase digital filter - filtfilt function (to minimize phase distortion and transients)
beta_zerophase_filt = filtfilt(beta_bandpass_filt, LFP_dt1_hpfilt);
% beta_zerophase_filt2 = [beta_zerophase_filt, 0];

% compute power spectral density (PSD) of LFP segment(s) - pspectrum function
[betaPxx,betaFxx] = pspectrum(beta_zerophase_filt,fs,'FrequencyLimits',[0 50],'FrequencyResolution',3); % 2 or 3
figure;
pspectrum(beta_zerophase_filt,fs,"spectrogram","FrequencyLimits",[0 50], "FrequencyResolution", 3)
figure
pspectrum(beta_zerophase_filt,fs,"spectrogram","FrequencyLimits",[0 50], "TimeResolution", 0.50)

% normalize using the common logarithm via decibel conversion - pow2db function
betaPxxP = pow2db(betaPxx);

% plot normalized PSD of beta band
figure;
%plot(freq, 10*log10(power));
plot(betaFxx, betaPxx);
xlim([0 50])  % Limit x-axis to 50 Hz for better visualization
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density of LFP beta band (13-30 Hz)');


%% Hilbert transorm of beta-filtered signal
% %%% fix this

% % compute instantaneous phase and frequency of the LFP beta band via Hilbert transform - hilbert function
% hilbert_transformed = hilbert(beta_zerophase_filt);
% inst_phase = angle(hilbert_transformed);
% inst_freq = diff(unwrap(inst_phase))/(2*pi*t_step);
% 
% figure;
% %plot(ts_LFP,beta_zerophase_filt);
% plot(ts_LFP(2:end), inst_freq);
% %xlim([0 round(max(ts_LFP))])
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% title('Instantaneous Frequency of Filtered LFP Beta Band (13-30 Hz)');


%% compute PSD of LFP, normalized via decibel conversion, and plot for visualization 

% simplify
% compute power spectral density (PSD) of LFP segment(s) - pspectrum function
[lfpPxx,lfpFxx] = pspectrum(LFP_dt1_hpfilt,fs,'FrequencyLimits',[0 50],'FrequencyResolution',3); 

figure;
pspectrum(LFP_dt1_hpfilt,fs,"spectrogram","FrequencyLimits",[0 50], "FrequencyResolution", 3)
figure
pspectrum(LFP_dt1_hpfilt,fs,"spectrogram","FrequencyLimits",[0 50], "TimeResolution", 0.50)

% normalize using the common logarithm via decibel conversion - pow2db function
lfpPxx = pow2db(lfpPxx);

% plot PSD for visualization
figure;
%plot(freq, 10*log10(power));
plot(lfpFxx, lfpPxx);
xlim([0 50])  % Limit x-axis to 50 Hz for better visualization
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density of Filtered LFP');


%% LFP time-frequency domain analysis
% bursting analysis, based on Torrecillos et al., 2018 (test for LFPs in 1 trial / STN depth)
% time-frq transform via complex Morlet wavelets (f0/σf = 7, f0 = 1-45 Hz, steps of 0.25 Hz).

% define 'lfp_data' as a LFP data matrix with dimensions [samples x trials]
lfp_data_test = LFP_dt1_hpfilt;

% Frequency-Time Decomposition using complex Morlet wavelets
f0 = 1:0.25:45; % center frequency range from 1-45Hz in steps of 0.25Hz
time = (1:length(lfp_data_test)) / fs;
cwt_power = zeros(length(f0), length(lfp_data_test));

for i = 1:length(f0)
    % create a sinusoidal wave modulated by a Gaussian window
    sinusoidal_wave = exp(2 * pi * 1i * f0(i) .* time); % complex cosine component of wavelet with an imaginary sine component
    Guassian_window = exp(-time.^2 / (2 * (f0(i) / 7)^2)); % tapers the sinusoid; (f0(i) / 7)^2 controls width of Gaussian
    wavelet = sinusoidal_wave .* Guassian_window; % applies wavelet across f0 range to create filters convolved with the LFP signal
    cwt_power(i, :) = abs(conv(lfp_data_test, wavelet, 'same')).^2; % square to compute power
end

% time-domain bin size determined by 'VoicesPerOctave' (higher VPO = finer freq. res; lower time res?)
% continuous wavelet transform
figure;
cwt(lfp_data_test,fs,'FrequencyLimits',[1 125],'VoicesPerOctave',14); % try 14, mults of 7
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Continuous Wavelet Transform');

% Normalize power for each frequency band
mean_power = mean(cwt_power, 2);
std_power = std(cwt_power, 0, 2);
normalized_power = (cwt_power - mean_power) ./ std_power;

% Beta Peak Selection
beta_range = find(f0 >= 13 & f0 <= 30);
[~, peak_idx] = max(mean(normalized_power(beta_range, :), 2));
beta_peak_frequency = f0(beta_range(peak_idx));

% Compute beta power time courses for each trial using a 6-Hz-wide frequency band centered on the beta peak frequency
beta_power_time_courses = mean(normalized_power(beta_range(peak_idx) + (-3:3), :), 1);

% Plot the beta power time course
figure('Renderer', 'opengl');
plot(time, beta_power_time_courses);
xlabel('Time (s)');
ylabel('Normalized Beta Power');
title('LFP Beta Power (untrimmed)')

%% Trim beta power lfp and time data based on known/computed offset

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
cwt(lfp_data_L_trimmed,fs,'FrequencyLimits',[1 125],'VoicesPerOctave',14); % try 14, mults of 7'
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Continuous Wavelet Transform');

% Plot the scalogram
[cwtCoeffs, frequencies] = cwt(lfp_data_L_trimmed, fs);
figure;
surface(time_trimmed, frequencies, abs(cwtCoeffs).^2);
shading interp;
axis tight;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('CWT Scalogram');
% colorbar;
% ylabel(colorbar, 'Power (µV^2/Hz)');


%% Beta Burst Detection and Plotting for trimmed LFP

% Compute mean beta power
mean_beta_power = mean(beta_power_time_courses);

% Set threshold at the 75th percentile of the mean beta power
threshold = prctile(mean_beta_power, 75);

% Define beta bursts as time points exceeding the threshold for more than 3 oscillatory cycles (3 oscillatory cycles = ~100ms)
min_duration = round(0.1 * fs); % Minimum duration of a burst in samples (100 ms)
% min_duration = round(0.066 * fs);  % Reduce minimum duration to 2 osc. cycles(66 ms)

above_threshold = beta_power_time_courses > threshold;
above_threshold = above_threshold(:);  % Ensure it's a column vector
bursts = zeros(size(above_threshold));  % Initialize bursts array
[start_indices, end_indices] = findConsecutiveIndices(above_threshold);  % Find start and end indices of consecutive regions above threshold
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


% Plot the beta power time course and beta burst overlays using xregion (test with L_streamOfInt_s1)
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

%% Compute and Plot PSD on trimmed LFP

pspectrum(lfp_data_L_trimmed,fs,'FrequencyLimits',[0 50],'FrequencyResolution',3) % 2 or 3
[bPxx,bFxx] = pspectrum(lfp_data_L_trimmed,250,'FrequencyLimits',[0 50],'FrequencyResolution',3); %
bPxxP = pow2db(bPxx);

% Plot the power spectral density
figure;
%plot(freq, 10*log10(power));
plot(bFxx, bPxxP);
xlim([0 50])  % Limit x-axis to 50 Hz for better visualization
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density of LFP');










%% Functions

function [beta_power_time_courses, beta_peak_frequency, bursts, time] = burstingAnalysis(lfp_data, fs)
% time-frq transform via complex Morlet wavelets (f0/σf = 7, f0 = 1-45 Hz, steps of 0.25 Hz)

% Define frequency range for Morlet wavelets
f0 = 1:0.25:45; % center frequency range from 1:45Hz in steps of 0.25Hz

% Frequency-Time Decomposition using complex Morlet wavelets
time = (1:length(lfp_data)) / fs;
cwt_power = zeros(length(f0), length(lfp_data));

for i = 1:length(f0)
    % create a sinusoidal wave modulated by a Gaussian window with a constant ratio of f0/σf = 7; σf = standard deviation of Gaussian in the freq domain
    sinusoidal_wave = exp(2 * pi * 1i * f0(i) .* time); % complex cosine component of wavelet with an imaginary sine component
    Guassian_window = exp(-time.^2 / (2 * (f0(i) / 7)^2)); % tapers the sinusoid; (f0(i) / 7)^2 controls width of Gaussian
    wavelet = sinusoidal_wave .* Guassian_window; % applies wavelet across f0 range to create filters convolved with the LFP signal
    cwt_power(i, :) = abs(conv(lfp_data, wavelet, 'same')).^2; % square to compute power
end
conv()
% time-domain bin size determined by 'VoicesPerOctave' (higher VPO = finer freq. res; lower time res?)
% continuous wavelet transform
figure;
cwt(lfp_data,fs,'FrequencyLimits',[1 125],'VoicesPerOctave',14); % try 14, mults of 7
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Continuous Wavelet Transform (untrimmed)');

% Normalize power for each frequency band
mean_power = mean(cwt_power, 2); % 2nd dim is the time dimension
std_power = std(cwt_power, 0, 2);
normalized_power = (cwt_power - mean_power) ./ std_power;

% Beta Peak Selection
beta_range = find(f0 >= 13 & f0 <= 30);
[~, peak_idx] = max(mean(normalized_power(beta_range, :), 2));
beta_peak_frequency = f0(beta_range(peak_idx));

% Burst Detection
% Compute beta power time courses for each trial using a 6-Hz-wide frequency band centered on the beta peak frequency
beta_power_time_courses = mean(normalized_power(beta_range(peak_idx) + (-3:3), :), 1);

% Compute mean beta power
mean_beta_power = mean(beta_power_time_courses);

% Set threshold at the 75th percentile of the mean beta power
threshold = prctile(mean_beta_power, 75);

% Define beta bursts as time points exceeding the threshold for more than 3 oscillatory cycles (3 oscillatory cycles = ~100ms)
min_duration = round(0.1 * fs); % Minimum duration of a burst in samples (100 ms)

above_threshold = beta_power_time_courses > threshold;
above_threshold = above_threshold(:);  % Ensure it's a column vector
bursts = zeros(size(above_threshold));  % Initialize bursts array
[start_indices, end_indices] = findConsecutiveIndices(above_threshold);  % Find start and end indices of consecutive regions above threshold
for k = 1:length(start_indices)
    if (end_indices(k) - start_indices(k) + 1) >= min_duration
        bursts(start_indices(k):end_indices(k)) = 1;  % Mark burst regions
    end
end

% Adjust burst indices based on the trimmed length of the time array
valid_indices = (start_indices <= length(time)) & (end_indices <= length(time));
start_indices = start_indices(valid_indices);
end_indices = end_indices(valid_indices);

end


function [start_indices, end_indices] = findConsecutiveIndices(logical_array)
% Find start and end indices of consecutive true regions in a logical array
d = diff([0; logical_array(:); 0]);
start_indices = find(d > 0);
end_indices = find(d < 0) - 1;
end


function dlcDataInterpolated = interpolateDLCData(dlcData, tsOriginal, tsTarget)
% Function to interpolate DLC data from original timestamps to target LFP timestamps
variableNames = dlcData.Properties.VariableNames;
dlcDataInterpolated = table();
for varIdx = 1:length(variableNames)
    dlcDataInterpolated.(variableNames{varIdx}) = interp1(tsOriginal, dlcData.(variableNames{varIdx}), tsTarget, 'linear', 'extrap');
end
end


function plotKinematicsAndLFP(ts, dlcData, lfpData, streamIndex, sessionID)
% Function to plot kinematics and LFP data, and compute PSD
figure;
tiledlayout(2, 1, "TileSpacing", "tight");
nexttile;
plot(ts, smoothdata(dlcData.fTip1_x, 'gaussian', 70));
xlabel('Time (s)');
ylabel('fTip1, X deflection');
nexttile;
plot(ts, lfpData);
xlabel('Time (s)');
ylabel('LFP (uV)');
sgtitle(sprintf('Kinematics and LFP for Stream %d, Session %s', streamIndex, sessionID));

% Power Spectral Density Analysis
figure;
[p, f] = pspectrum(lfpData, AO_LFP_fs, 'FrequencyLimits', [0 50], 'FrequencyResolution', 3);
plot(f, pow2db(p));
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title(sprintf('PSD of LFP Stream %d, Session %s', streamIndex, sessionID));
end


function bursts = betaBurstDetection(cfs, freqs, time, beta_range, fs)
% Compute mean and std of power in the beta range
beta_power = mean(abs(cfs(beta_range,:)).^2, 1);
mean_beta_power = mean(beta_power);
std_beta_power = std(beta_power);

% Detect bursts as periods where power is significantly above mean
burst_threshold = mean_beta_power + 2 * std_beta_power; % Example threshold
bursts = beta_power > burst_threshold;
end


















