% Clinical LFP analyses - LSTN

%% Combine Percept LFP with DLC Video

% Data ingredients:
% 1) LFP data - JSON Session Reports (multiple rows (stream recordings) per report [metadata informs ID of row])
% 2) Movement data (.mat files) and Movement Indices (.csvs)

clear; close all; clc;

%% Directory set-up - Navigate b/t machines
pcName = getenv('COMPUTERNAME');

switch pcName

    case 'DSKTP-JTLAB-EMR'   %%% ER Desktop

        mainDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\DLC_LFP';

    case 'NSG-M-FQBPFK3'     %%% ER PC 1

        mainDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\DLC_LFP';

    case 'NSG-M-H8J3X34'     %%% ER PC 2

        mainDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\DLC_LFP';
end

%% Analyze data isolated by casedate and hemisphere

% Define casedate_hem
casedate_hem = '09_12_2023_LSTN_v2';
% casedate_hem = '09_12_2023_RSTN_v2';

mainDir2 = [mainDir , filesep , casedate_hem];

cd(mainDir2)

%% Navigate to LFP time domain data in JSON Session Reports

% Load JSON Session Reports
% JSON_name1 = 'Report_Json_Session_Report_20230912T115956.json';
% JSON_name2 = 'Report_Json_Session_Report_20230912T115939.json';

% Create array of JSON Session Report file names
JSON_filenames = {JSON_name1, JSON_name2};

% Initialize structure to store the first row of outTAB for each JSON file
session_StartTimes = struct();

% Define perceive_ecg function params and add function paths
fs = 250;
plotit = 0; % 0 = don't plot, 1 = plot
addpath 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\perceive-master'

% loop through each JSON file
for json_i = 1:length(JSON_filenames)

    % load current JSON Session Report
    currentJSON_name = JSON_filenames{json_i};
    currentJSON = jsondecode(fileread(currentJSON_name));

    % calculate the date/time for the current JSON
    sessDate_field_1 = currentJSON.SessionDate;

    % Convert the string to a datetime object
    dateTimeObj_1 = datetime(sessDate_field_1, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

    % Set the time zone to UTC (Coordinated Universal Time) since the input is in UTC
    dateTimeObj_1.TimeZone = 'UTC';

    % Convert to Mountain Time
    dateTimeObj_Mountain_1 = datetime(dateTimeObj_1, 'TimeZone', 'America/Denver');

    % Extract the time component in AM/PM format
    timeComponent_AMPM = datetime(dateTimeObj_Mountain_1,'Format','hh:mm:ss a z');
    timeComponent_DATE = datetime(dateTimeObj_Mountain_1,'Format','dd-MMM-yyyy');

    % convert timedomain to table
    BSTD_1 = currentJSON.BrainSenseTimeDomain; % struct
    BSTD_1_table = struct2table(BSTD_1);

    % plot raw and ecg-filtered time domain data for each row of current JSON file
    for BSTD_i = 1:size(BSTD_1_table, 1) % Loop through each row in current JSON

        % filter out ECG for each row of time domain data in current JSON
        tempData_1 = transpose(BSTD_1_table.TimeDomainData{BSTD_i}); % Transpose raw data for current row
        ecg = perceive_ecg(tempData_1, fs, plotit);

    end

    % navigate to time domain data
    [outTAB] = getBSLFPtimes(BSTD_1_table.FirstPacketDateTime);

    % Determine session start time - display first row of outTAB
    disp(outTAB(1,:));

    % Store the first row of outTAB in the session_StartTimes structure
    session_StartTimes.(sprintf('File%d', json_i)) = outTAB(1,:);

end

% JSON_name1 (...5956.json), FullNAT {'12-Sep-2023 10:17:12'},  Off Med
% JSON_name2 (...5939.json), FullNAT {'12-Sep-2023 11:31:25'},  On Med


%% Process OFF Med JSON Session Report 1st

% load JSON Session Report
js_1_name = 'Report_Json_Session_Report_20230912T115956.json' % Off Med
js_1 = jsondecode(fileread(js_1_name));

js_2_name = 'Report_Json_Session_Report_20230912T115939.json'; % On Med
js_2 = jsondecode(fileread(js_2_name));

% calculate the date/time
sessDate_field_1 = js_1.SessionDate;

% Convert the string to a datetime object
dateTimeObj_1 = datetime(sessDate_field_1, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

% Set the time zone to UTC (Coordinated Universal Time) since the input is in UTC
dateTimeObj_1.TimeZone = 'UTC';

% Convert to Mountain Time
dateTimeObj_Mountain_1 = datetime(dateTimeObj_1, 'TimeZone', 'America/Denver');

% convert BrainSense timedomain data to table
BSTD_1 = js_1.BrainSenseTimeDomain; % struct
BSTD_1_table = struct2table(BSTD_1); % table

% navigate to time domain data
BSLFP_times_1 = getBSLFPtimes(BSTD_1_table.FirstPacketDateTime);
BSLFP_times_1_table = cell2table(BSLFP_times_1.FullNAT);
disp(BSLFP_times_1(1,3)); % ensure correct file by start-time info


%% test - filter out ecg for first LFP stream

plot(BSTD_1_table.TimeDomainData{1}); % plot raw time domain data for row 1

streamOfInt_transposed = transpose(BSTD_1_table.TimeDomainData{1}); % transpose raw data for row 1
ecg = perceive_ecg(streamOfInt_transposed,fs,plotit); % run perceive_ecg function

hold on
plot(ecg.cleandata) % plot ecg-filtered time domain data for row 1
hold off


%% Load Left & Right STN BrainSense LFP streaming sessions in OFF Med JSON Session Report

% BrainSense timedomain data table
streaming_TAB_1 = BSTD_1_table;

% Trim by STN of interest - LSTN, RSTN
stream_LEFT_1 = streaming_TAB_1(contains(streaming_TAB_1.Channel,'LEFT'),:); % L STN, R body
stream_RIGHT_1 = streaming_TAB_1(contains(streaming_TAB_1.Channel,'RIGHT'),:); % R STN, L body

% Determine duration (in seconds) of each stream
stream_LEFT_1_times = getDAYtimes(stream_LEFT_1.FirstPacketDateTime, stream_LEFT_1.TimeDomainData);
LEFT_sessDurations_seconds_1 = stream_LEFT_1_times.Duration;

stream_RIGHT_1_times = getDAYtimes(stream_RIGHT_1.FirstPacketDateTime, stream_RIGHT_1.TimeDomainData);
RIGHT_sessDurations_seconds_1 = stream_RIGHT_1_times.Duration;


%% Separate LFP streams and corresponding DLC sessions of interest by condition
% 9/12/2023 BrainSense LFP streams (12) in JSON Session Report 1 - Off Med

% baseline
% js1_Row	sessDur(s)	DLC_sessID	DLC_sessDur(s) 	Hemisphere	Notes: case videos (10) for JSON Session Report 1 - Off Med
% 1	        32	        NA	            NA	        L	        % • L_baseline, (Off Med, Off Stim @ 0 mA)
% 2	        32	        NA		        NA          R	        % • R_baseline, (Off Med, Off Stim @ 0 mA)

% Off Stim
% js1_Row	sessDur(s)	DLC_sessID	DLC_sessDur(s) 	Hemisphere	Notes: case videos (10) for JSON Session Report 1 - Off Med
% 3	        47	        session001	    56	        L	        % • session001, set1 (Off Med, Off Stim @ 0 mA)
% 4	        41	        session002		40          R	        % • session002, set1 (Off Med, Off Stim @ 0 mA)
% 5	        38	        session003	    48	        L	        % • session003, set2 (Off Med, Off Stim @ 0 mA)
% 6	        34	        session004		44          R	        % • session004, set2 (Off Med, Off Stim @ 0 mA)

% Ramp
% js1_Row	sessDur(s)	DLC_sessID	DLC_sessDur(s) 	Hemisphere	Notes: case videos (10) for JSON Session Report 1 - Off Med
% 7	        322	        session005	    336	        L	        % • session005 (Off Med, Stim Ramping)
% 8	        251	        session006		262         R	        % • session006 (Off Med, Stim Ramping)

% On Stim
% js1_Row	sessDur(s)	DLC_sessID	DLC_sessDur(s) 	Hemisphere	Notes: case videos (10) for JSON Session Report 1 - Off Med
% 9	        36	        session007	    45	        L	        % • session007, set1 (Off Med, On Stim @ max mA)
% 10	    29	        session008		37          R	        % • session008, set1 (Off Med, On Stim @ max mA)
% 11	    28	        session009	    35	        L	        % • session009, set2 (Off Med, On Stim @ max mA)
% 12	    26	        session010		36          R       	% • session010, set2 (Off Med, On Stim @ max mA)


%% Corresponding movement data

DLC_csv_dir = [mainDir2 , filesep , 'csv folder']; % contains dlc label timeseries data as csv files
DLC_mat_dir = [mainDir2 , filesep , 'mat folder']; % contains dlc label timeseries data as mat files
DLC_video_dir = [mainDir2 , filesep , 'video folder', filesep , 'Converted for dual GUI']; % contains labeled videos + MovementIndex csv


%% Isolate dlc outputs of interest

% Generate list of dlc-video-labeled MAT files
cd(DLC_mat_dir)
mainMat = dir('*.mat');
mainMAT2 = {mainMat.name};

% Generate list of Motor Index CSVs
cd(DLC_video_dir)
moveCSV = dir('*.csv');
moveCSV2 = {moveCSV.name};
% Generate list of Motor Index CSVs (filters for CSVs that contain 'Move' string)
moveCSV = moveCSV2(contains(moveCSV2,'Move'));

% EUC indicies
cd(mainDir2)
eucINDICIES = readtable("EUC_Indicies_1.xlsx");

% Find all rows in eucINDICIES where the 'moveID' is 'Tablet'
tablet_indices = find(strcmp(eucINDICIES.moveID, 'Tablet'));

% Define row indices in eucINDICIES for each 'moveID'
HandMov_indices = find(strcmp(eucINDICIES.moveID, 'HAND OC'));
PronSup_indices = find(strcmp(eucINDICIES.moveID, 'HAND PS'));
FlexExtend_indices = find(strcmp(eucINDICIES.moveID, 'ARM EF'));


%% Isolate specific rows of interest in JSON Session Report 1 - Off Med

% Left STN, Right body
% Specific Rows of interest
L_rowfromTab_bsln = 1;
L_streamOfInt_bsln = stream_LEFT_1.TimeDomainData{L_rowfromTab_bsln}; % L_baseline (Off Med, Off Stim @ 0 mA)

L_rowfromTab_s1 = 3;
L_streamOfInt_s1 = stream_LEFT_1.TimeDomainData{L_rowfromTab_s1}; % L_set1 (Off Med, Off Stim @ 0 mA)

L_rowfromTab_s3 = 5;
L_streamOfInt_s3 = stream_LEFT_1.TimeDomainData{L_rowfromTab_s3}; % L_set2 (Off Med, Off Stim @ 0 mA)

L_rowfromTab_s5 = 7;
L_streamOfInt_s5 = stream_LEFT_1.TimeDomainData{L_rowfromTab_s5}; % L_ramp (Off Med, Stim Ramping)

L_rowfromTab_s7 = 9;
L_streamOfInt_s7 = stream_LEFT_1.TimeDomainData{L_rowfromTab_s7}; % L_set1 (Off Med, On Stim @ max mA)

L_rowfromTab_s9 = 11;
L_streamOfInt_s9 = stream_LEFT_1.TimeDomainData{L_rowfromTab_s9}; % L_set2 (Off Med, On Stim @ max mA)


% Right STN, Left body
% Specific Rows of interest
R_rowfromTab_bsln = 2;
R_streamOfInt_bsln = stream_RIGHT_1.TimeDomainData{R_rowfromTab_bsln}; % R_baseline, (Off Med, Off Stim @ 0 mA)

R_rowfromTab_s2 = 4;
R_streamOfInt_s2 = stream_RIGHT_1.TimeDomainData{R_rowfromTab_s2}; % R_set1 (Off Med, Off Stim @ 0 mA)

R_rowfromTab_s4 = 6;
R_streamOfInt_s4 = stream_RIGHT_1.TimeDomainData{R_rowfromTab_s4}; % R_set2 (Off Med, Off Stim @ 0 mA)

R_rowfromTab_s6 = 8;
R_streamOfInt_s6 = stream_LEFT_1.TimeDomainData{R_rowfromTab_s6}; % R_ramp (Off Med, Stim Ramping)

R_rowfromTab_s8 = 10;
R_streamOfInt_s8 = stream_LEFT_1.TimeDomainData{R_rowfromTab_s8}; % R_set1 (Off Med, On Stim @ max mA)

R_rowfromTab_s10 = 12;
R_streamOfInt_s10 = stream_LEFT_1.TimeDomainData{R_rowfromTab_s10}; % R_set2 (Off Med, On Stim @ max mA)


%% Sanity Check Plotting

% L_streamsofInt_OffMed = {L_streamOfInt_bsln, L_streamOfInt_s1, L_streamOfInt_s3, L_streamOfInt_s5, L_streamOfInt_s7, L_streamOfInt_s9}; % w/ baseline
L_streamsofInt_OffMed = {L_streamOfInt_s1, L_streamOfInt_s3, L_streamOfInt_s7, L_streamOfInt_s9}; % w/o baseline or ramp

% non-ecg-filtered
% maybe hp filter

% Titles for each plot in L_streamsofInt_OffMed
plotTitles = {
    'LSTN set1 (Off Med, Off Stim @ 0 mA)',
    'LSTN set2 (Off Med, Off Stim @ 0 mA)',
    'LSTN set1 (Off Med, On Stim @ max mA)',
    'LSTN set2 (Off Med, On Stim @ max mA)'};

figure
for L_i = 1:length(L_streamsofInt_OffMed)
    ts_LFP = 0:1/250:(height(L_streamsofInt_OffMed{L_i})-1)/250;

    subplot(4,1,L_i)
    plot(ts_LFP, L_streamsofInt_OffMed{L_i})
    title(plotTitles{L_i})
    xlabel('Time (s)')
    ylabel('LFP Amplitude (uV)')
    grid on
end


% R_streamsofInt_OffMed = {R_streamOfInt_bsln, R_streamOfInt_s2, R_streamOfInt_s4, R_streamOfInt_s6, R_streamOfInt_s8, R_streamOfInt_s10}; % w/ baseline
R_streamsofInt_OffMed = {R_streamOfInt_s2, R_streamOfInt_s4, R_streamOfInt_s8, R_streamOfInt_s10}; % w/o baseline or ramp

% Titles for each plot in L_streamsofInt_OffMed
plotTitles = {
    'RSTN set1 (Off Med, Off Stim @ 0 mA)',
    'RSTN set2 (Off Med, Off Stim @ 0 mA)',
    'RSTN set1 (Off Med, On Stim @ max mA)',
    'RSTN set2 (Off Med, On Stim @ max mA)'};

figure
for R_i = 1:length(R_streamsofInt_OffMed)
    ts_LFP = 0:1/250:(height(R_streamsofInt_OffMed{R_i})-1)/250;

    subplot(4,1,R_i)
    plot(ts_LFP, R_streamsofInt_OffMed{R_i})
    title(plotTitles{R_i})
    xlabel('Time (s)')
    ylabel('LFP Amplitude (uV)')
    grid on
end


%% Filter, Hilbert Transform (test w/ L_streamOfInt_s1 (Off Med, Off Stim @ 0 mA))

% L_streamOfInt_s1 (LFP stream 3) corresponds with DLC session001 (movement video 1), (Off Med, Off Stim @ 0 mA)
% Define LFP stream of interest:
L_streamOfInt_OffOff = L_streamOfInt_s1; % raw LFP

fs = 250;
t = 1/250;

% Target signal (LFP) sampling rate at 250 Hz (samples per sec)
ts_LFP = 0:1/250:(height(L_streamOfInt_OffOff)-1)/250;

ecgClean = perceive_ecg(transpose(L_streamOfInt_OffOff),250,0);

cleanLFP_OffOff = ecgClean.cleandata;

% Calculate instantaneous phase and freq. using Hilbert transform
hilbert = hilbert(cleanLFP_OffOff);
inst_phase = angle(hilbert);
inst_freq = diff(unwrap(inst_phase))/(2*pi*t);

% Filter instantaneous frequency for 13-30 Hz
beta_bandpass_filt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',13,'HalfPowerFrequency2',30, ...
    'SampleRate',fs);
inst_freq_filtered = filtfilt(beta_bandpass_filt, inst_freq);
inst_freq_filtered2 = [inst_freq_filtered , 0];

figure;
plot(ts_LFP,inst_freq_filtered2);
xlim([0 round(max(ts_LFP))])
ylabel('frequency (Hz)')
xlabel('time (seconds)')
title('Instantaneous Frequency of Filtered LFP Beta Band (13-30 Hz)');

% Apply the bandpass filter
betaLFP = filtfilt(beta_bandpass_filt, L_streamOfInt_OffOff);

% % Compute the power of the beta band signal
% betaLFP_power = abs(hilbert(betaLFP)).^2;
% 
% % Calculate the mean and standard deviation
% mean_beta_power = mean(betaLFP_power);
% std_beta_power = std(betaLFP_power);

%% LFP filtering
%%% fix this %%%

% remove 60 Hz line noise and harmonics via 2nd-order infinite impulse response (IIR) notch and comb filters
% Notch Filter centered at 60Hz
f_notch = 60;  % Notch frequency
Q = 35;   % Quality factor, q = ω0/bw where ω0 is the frequency to remove from the signal
[b_notch, a_notch] = iirnotch(f_notch/(fs/2), f_notch/(fs/2)/Q); % 2nd order infinite impulse response (IIR) notch filter
LFP_notchfilt = filtfilt(b_notch, a_notch, cleanLFP_OffOff); % Apply notch filter to remove 60Hz line noise

% Comb Filter centered at 60Hz and its harmonics (120Hz, 180Hz, 240Hz, 300Hz, 360Hz, etc.) 
bw = f_notch / Q;  % Bandwidth, bw = (fo/(fs/2))/q;
n_harmonics = floor((fs/2) / f_notch);  % Number of harmonics within Nyquist frequency
[b_comb, a_comb] = iircomb(floor(fs/f_notch), bw/(fs/2), 'notch');  % Design comb filter
LFP_combfilt = filtfilt(b_comb, a_comb, cleanLFP_OffOff); % Apply comb filter to remove harmonics


% attenuate low frequency components in LFP via 4th order Butterworth or IIR high-pass filter
% High-pass filter 
cutoff_freq = 2; % Hz, attenuate frequency components below cutoff frequency 
hp_filter = designfilt('highpassiir', 'FilterOrder', 4, ... 
    'HalfPowerFrequency', cutoff_freq, 'SampleRate', fs); % Design 4th order Butterworth or IIR hp-filter
LFP_hpfilt = filtfilt(hp_filter, LFP_combfilt); % Apply high-pass filter to notch/comb-filtered LFP signal

% Plot the filtered LFP signal
figure;
plot(ts_LFP, LFP_hpfilt);
xlabel('Time (s)');
ylabel('Amplitude');
title('High-Pass Filtered LFP Signal');


% compute power spectral density (PSD) of LFP segment(s) - pspectrum function
[fPxx,fFxx] = pspectrum(LFP_hpfilt,fs,'FrequencyLimits',[0 100],'FrequencyResolution',3); 

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
beta_zerophase_filt = filtfilt(beta_bandpass_filt, LFP_hpfilt);
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


%% Compute and Plot Power Spectral Density (PSD), (test w/ L_streamOfInt_s1 (Off Med, Off Stim @ 0 mA))
% Use pspectrum to compute the power spectral density
% [power, freq] = pspectrum(cleanLFP_s1, fs);

pspectrum(cleanLFP_OffOff,fs,'FrequencyLimits',[0 100],'FrequencyResolution',3) % 2 or 3
[bPxx,bFxx] = pspectrum(cleanLFP_OffOff,250,'FrequencyLimits',[0 100],'FrequencyResolution',3); %
bPxxP = pow2db(bPxx);

% Plot the power spectral density
figure;
%plot(freq, 10*log10(power));
plot(bFxx, bPxxP);
xlim([0 100])  % Limit x-axis to 50 Hz for better visualization
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density of LFP (untrimmed)');


%% bursting analysis, based on Torrecillos et al., 2018 (test w/ L_streamOfInt_s1 (Off Med, Off Stim @ 0 mA))
% time-frq transform via complex Morlet wavelets (f0/σf = 7, f0 = 1-45 Hz, steps of 0.25 Hz).

% define 'lfp_data' as a LFP data matrix with dimensions [samples x trials]
lfp_data_L = LFP_hpfilt;

fs = 250; % sampling rate
f0 = 1:0.25:45; % center frequency range from 1-45Hz in steps of 0.25Hz

% Frequency-Time Decomposition using complex Morlet wavelets
time = (1:length(lfp_data_L)) / fs;
cwt_power = zeros(length(f0), length(lfp_data_L));

for i = 1:length(f0)
    % create a sinusoidal wave modulated by a Gaussian window
    sinusoidal_wave = exp(2 * pi * 1i * f0(i) .* time); % complex cosine component of wavelet with an imaginary sine component
    Guassian_window = exp(-time.^2 / (2 * (f0(i) / 7)^2)); % tapers the sinusoid; (f0(i) / 7)^2 controls width of Gaussian
    wavelet = sinusoidal_wave .* Guassian_window; % applies wavelet across f0 range to create filters convolved with the LFP signal
    cwt_power(i, :) = abs(conv(lfp_data_L, wavelet, 'same')).^2; % square to compute power
end

% time-domain bin size determined by 'VoicesPerOctave' (higher VPO = finer freq. res; lower time res?)
% continuous wavelet transform
figure;
cwt(lfp_data_L,fs,'FrequencyLimits',[1 125],'VoicesPerOctave',14); % try 14, mults of 7
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
time = time(1:drop_off_index);

%% continuous wavelet transform on trimmed LFP

lfp_data_L_trimmed = lfp_data_L(1:drop_off_index);
figure;
cwt(lfp_data_L_trimmed,fs,'FrequencyLimits',[1 125],'VoicesPerOctave',14); % try 14, mults of 7'
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Continuous Wavelet Transform');

% Plot the scalogram
[cwtCoeffs, frequencies] = cwt(lfp_data_L_trimmed, fs);
figure;
surface(time, frequencies, abs(cwtCoeffs).^2);
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
plot(time, beta_power_time_courses);
xlabel('Time (s)');
ylabel('Normalized Beta Power');
title('LFP Beta Power')


% Plot the beta power time course and beta burst overlays using xregion (test with L_streamOfInt_s1)
% https://www.mathworks.com/help/matlab/ref/xregion.html

figure('Renderer', 'opengl');
p = plot(time, beta_power_time_courses); % store handle to the beta power plot
hold on;

% Overlay the detected bursts on the same plot using xregion
for i = 1:length(start_indices)
    burst_start_time = time(start_indices(i));
    burst_end_time = time(end_indices(i));
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


%% Loop through each LFP stream in Left STN - Apply Filtering, PSD, CWT, and Beta Bursting code to all LFP streams
for i = 1:length(L_streamsofInt_OffMed)

    lfp_data_L = L_streamsofInt_OffMed{i};

    % Target signal (LFP) sampling rate at 250 Hz (samples per sec)
    ts_LFP = 0:1/250:(height(lfp_data_L)-1)/250;

    % Filter ECG
    ecgClean_L = perceive_ecg(transpose(lfp_data_L),250,0);
    cleanLFP_L = ecgClean_L.cleandata;

    % Beta Bursting Analysis
    [beta_power_time_courses, beta_peak_frequency, bursts, time] = burstingAnalysis(cleanLFP_L, fs);


    % Trimming: Define a threshold as a percentage of the maximum beta power
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
    time = time(1:drop_off_index);

    % Adjust burst indices based on the trimmed length of the time array
    valid_indices = (start_indices <= length(time)) & (end_indices <= length(time));
    start_indices = start_indices(valid_indices);
    end_indices = end_indices(valid_indices);

    % Trim the LFP data
    lfp_data_L_trimmed = cleanLFP_L(1:drop_off_index);

    % Compute power spectrum of LFP stream sampled at fs
    pspectrum(lfp_data_L_trimmed, 250, 'FrequencyLimits', [0 50],'FrequencyResolution',3); % 2 or 3
    [bPxx,bFxx] = pspectrum(lfp_data_L_trimmed,250,'FrequencyLimits',[0 50],'FrequencyResolution',3); %
    bPxxP = pow2db(bPxx); % power to decibel conversion

    % Plot the power spectral density (PSD)
    figure;
    %plot(freq, 10*log10(power));
    plot(bFxx, bPxxP);
    xlim([0 50])  % Limit x-axis to 50 Hz for better visualization
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title(sprintf('Power Spectral Density, LSTN LFP Stream %d', i));

    % Continuous wavelet transform on trimmed lfp
    figure;
    cwt(lfp_data_L_trimmed, fs, 'FrequencyLimits', [1 125], 'VoicesPerOctave', 14);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(sprintf('Continuous Wavelet Transform, LSTN LFP Stream %d', i));

    % Plot the beta power time course
    figure('Renderer', 'opengl');
    plot(time, beta_power_time_courses);
    xlabel('Time (s)');
    ylabel('Normalized Beta Power');
    title(sprintf('LFP Beta Power, LSTN LFP Stream %d', i))

    % Plot the beta power time course and beta burst overlays using xregion
    figure('Renderer', 'opengl');
    p = plot(time, beta_power_time_courses); % store handle to the beta power plot
    hold on;

    % Overlay the detected bursts on the same plot using xregion
    for j = 1:length(start_indices)
        if j <= length(time) % Check if the index is within the valid range of the time array
            burst_start_time = time(start_indices(j));
            burst_end_time = time(min(end_indices(j), length(time))); % Ensure the end index does not exceed the length of the time array
            xregion(burst_start_time, burst_end_time, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8]);
        end
    end

    legend('auto update', 'off'); % Disable automatic legend updates
    h = patch(NaN, NaN, [0.8, 0.8, 0.8], 'FaceAlpha', 0.3); % Create a dummy patch handle for the legend
    legend([p, h], 'Beta Power', 'Detected Bursts'); % Add the beta power plot handle to the legend

    xlabel('Time (s)');
    ylabel('Normalized Beta Power');
    title(sprintf('Beta Bursting Dynamics, LSTN LFP Stream %d', i));
    hold off;
end





%%  % Isolate LFP and Video of iterest (LFP stream 3, DLC session 1) - LSTN Off Med, Off Stim
% Kinematics and LFP - LSTN Off Med, Off Stim

% Off Stim
% js1_Row	sessDur(s)	DLC_sessID	DLC_sessDur(s) 	Hemisphere	Notes: case videos (10) for JSON Session Report 1 - Off Med
% 3	        47	        session001	    56	        L	        % • session001, set1 (Off Med, Off Stim @ 0 mA)
% 4	        41	        session002		40          R	        % • session002, set1 (Off Med, Off Stim @ 0 mA)

% L_streamOfInt_s1 (LFP stream 3) corresponds with DLC session001 (movement video 1), (Off Med, Off Stim @ 0 mA)

% Define LFP stream of interest:
L_streamOfInt_OffOff = L_streamOfInt_s1;

% Extract DLC video of interest
dlc_session_id_OffOff = eucINDICIES.videoID{1};

% Load mat file corresponding to DLC video session
cd(DLC_mat_dir)
load('dlcDAT_20230912_session001_idea08_resnet50_rightCam-0000DLC.mat','outDATA')
dlc_labDAT_OffOff = outDATA;

% Extract dlc Tablet start and stop frame index
tablet_start_index_OffOff = eucINDICIES.StartInd(tablet_indices(1));
tablet_stop_index_OffOff = eucINDICIES.StopInd(tablet_indices(1));

% Trim mat file corresponding to DLC video session (based on offset b/t Video start/stop to Tablet start/stop index)
dlc_trimmed_OffOff = dlc_labDAT_OffOff(tablet_start_index_OffOff:tablet_stop_index_OffOff,:); % 2821x40 table

%% Synchronize corresponding LFP stream (trimmed by Tablet start index) s.t. LFP_start index syncs with Tablet_start index

% Sampling Rate Conversions and Calculations
totalNumSamples_Vid_OffOff = height(dlc_trimmed_OffOff);
totalNumSecs_OffOff = totalNumSamples_Vid_OffOff/60; % 60 fps
totalNumSamples_LFP_OffOff = floor(totalNumSecs_OffOff*250); % 250 samples per second

% Original signal (video) sampling rate at 60 Hz (fps)
ts_DLC_OffOff = 0:1/60:(height(dlc_trimmed_OffOff)-1)/60;

% Target signal (LFP) sampling rate at 250 Hz (samples per sec)
ts_LFP_OffOff = 0:1/250:(height(L_streamOfInt_OffOff)-1)/250;

% Interpolate - upsample kinematic data to match LFP sampling rate
allColNs = dlc_trimmed_OffOff.Properties.VariableNames;
dlc_trimmed_OffOff_int = table;
for coli = 1:width(dlc_trimmed_OffOff)

    % Tmp col
    tmpCol = allColNs{coli};

    % % Interpolation
    % x_250 = interp1(ts_DLC, transpose(dlc_lab2use.(tmpCol)), ts_LFP, 'spline');
    x_250 = interp1(ts_DLC_OffOff, dlc_trimmed_OffOff.(tmpCol), ts_LFP_OffOff);

    dlc_trimmed_OffOff_int.(tmpCol) = transpose(x_250);
end

trimFrame_int_OffOff = linspace(tablet_start_index_OffOff,tablet_stop_index_OffOff, length(ts_LFP_OffOff));


% Kinematics and LFP plot (LFP stream 3, DLC session 1) - Off Med, Off Stim
figure;
tiledlayout(2,1,"TileSpacing","tight")

xTime = ts_LFP_OffOff;

% kinematics
fTip1_X_smooth = smoothdata(dlc_trimmed_OffOff_int.fTip1_x,'gaussian',70);
nexttile
plot(xTime,fTip1_X_smooth);
xlim([0 round(max(ts_LFP_OffOff))])
xlabel('Time (s)')
ylabel('fTip1, X deflection (mm)')
title('Kinematics: Finger Tip X Position');

% LFP
ecg = perceive_ecg(transpose(L_streamOfInt_OffOff),250,0);
nexttile
plot(xTime,ecg.cleandata);
xlim([0 round(max(ts_LFP_OffOff))])
xlabel('Time (s)')
ylabel('LFP Amplitude (µV)');
title('Time Synced STN LFP Data'); 

% continuous wavelet transform
figure;
cwt(ecg.cleandata,fs,'FrequencyLimits',[1 125],'VoicesPerOctave',21); % try 14, mults of 7
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Continuous Wavelet Transform');


%% Hand Open/Close segment (LFP stream 3, DLC session 1), LSTN Off Med, Off Stim - HOC

% Extract start and stop frame indices for 'Hand OC' in DLC session data
HOC_start_index_OffOff = eucINDICIES.StartInd(HandMov_indices(1));
HOC_stop_index_OffOff = eucINDICIES.StopInd(HandMov_indices(1));

% Extract the kinematic data within the 'Hand OC' start and stop frame indices 
dlc_trimmed_HOC_OffOff = dlc_labDAT_OffOff(HOC_start_index_OffOff:HOC_stop_index_OffOff,:);

% Smooth changes in fTip1_x position data (euclidian distances) within frame indices
fTip1_X_smooth_HOC_OffOff = smoothdata(dlc_trimmed_HOC_OffOff.fTip1_x, 'gaussian', 70);

% Normalize kinematic data
fTip1_X_normalized = fTip1_X_smooth_HOC_OffOff - min(fTip1_X_smooth_HOC_OffOff) + (0.01 * fTip1_X_smooth_HOC_OffOff);

% Time vector for kinematic data
xTime_Kinematics = linspace(0, length(fTip1_X_smooth_HOC_OffOff)/60, length(fTip1_X_smooth_HOC_OffOff));

% Plot normalized kinematic data over time
figure;
plot(xTime_Kinematics, fTip1_X_normalized);
xlim([0, round(max(xTime_Kinematics)-1)])
ylim([0, round(max(fTip1_X_normalized)+1)]) 
xlabel('Time (s)')
ylabel('Normalized fTip1, X deflection (mm)')
title('Kinematics: Finger Tip X Position');

%% Create synchronized Kinematics and LFP plot for HOC task (LFP stream 3, DLC session 1) - HOC
% 'trimFrame_int_OffOff' has the indices for the LFP data after being interpolated

% Find the closest matching indices in the interpolated kinematic data
[~, startInd_Kin] = min(abs(trimFrame_int_OffOff - HOC_start_index_OffOff));
[~, stopInd_Kin] = min(abs(trimFrame_int_OffOff - HOC_stop_index_OffOff));

% Trim the interpolated kinematic data to match the 'Hand OC' segment
dlc_trimmed_HOC_OffOff_int = dlc_trimmed_OffOff_int(startInd_Kin:stopInd_Kin, :);

% Find the corresponding indices in the LFP time vector
[~, startIND_LFP] = min(abs(ts_LFP_OffOff - ts_DLC_OffOff(HOC_start_index_OffOff)));
[~, endIND_LFP] = min(abs(ts_LFP_OffOff - ts_DLC_OffOff(HOC_stop_index_OffOff)));

% Extract the LFP data for Hand OC
LFP_HOC_OffOff = L_streamOfInt_OffOff(startIND_LFP:endIND_LFP);

% Calculate the corresponding time vector for the trimmed kinematic data
xTime_Kinematics_HOC = linspace(0, (stopInd_Kin - startInd_Kin + 1) / fs, stopInd_Kin - startInd_Kin + 1);

% Normalize the trimmed kinematic data so that the min is 0 and max is 50
fTip1_X_smooth_HOC = smoothdata(dlc_trimmed_HOC_OffOff_int.fTip1_x, 'gaussian', 70);
fTip1_X_normalized_HOC = (fTip1_X_smooth_HOC - min(fTip1_X_smooth_HOC)) * 25 / (max(fTip1_X_smooth_HOC) - min(fTip1_X_smooth_HOC));

% HOC Kinematics and LFP plot (LFP stream 3, DLC session 1) - Off Med, Off Stim
figure;
tiledlayout(2,1,"TileSpacing","tight")

% Plot normalized kinematic data 
figure;
nexttile
plot(xTime_Kinematics_HOC, fTip1_X_normalized_HOC);
xlim([0, round(max(xTime_Kinematics_HOC)-1)])
ylim([0, round(max(fTip1_X_normalized_HOC)+1)]) 
xlabel('Time (s)');
ylabel('fTip1, X deflection (mm)');
title('Normalized Kinematics: Hand Movements');

% Plot LFP
xTime_LFP = linspace(0, length(LFP_HOC_OffOff)/250, length(LFP_HOC_OffOff));
ecg = perceive_ecg(LFP_HOC_OffOff', 250, 0);
nexttile
plot(xTime_LFP, ecg.cleandata);
xlim([0, round(max(xTime_LFP)-1)])
% ylim([round(min(ecg.cleandata)-1), round(max(ecg.cleandata)+1)])
ylim([-45, 27])
ylabel('LFP Amplitude (µV)');
title('Time Synced STN LFP Data');


%% Compute and Plot Power Spectral Density (PSD), (test w/ L_streamOfInt_s1 (Off Med, Off Stim @ 0 mA)) - HOC

% Filter ecg
ts_LFP_HOC_OffOff = 0:1/250:(height(LFP_HOC_OffOff)-1)/250;
ecgClean_HOC_OffOff = perceive_ecg(transpose(LFP_HOC_OffOff),250,0);
cleanLFP_HOC_OffOff = ecgClean_HOC_OffOff.cleandata;

pspectrum(cleanLFP_HOC_OffOff,250,'FrequencyLimits',[0 50],'FrequencyResolution',3) % 2 or 3
[bPxx_HOC_OffOff, bFxx_HOC_OffOff] = pspectrum(cleanLFP_HOC_OffOff,250,'FrequencyLimits',[0 50],'FrequencyResolution',3); %
bPxxP_HOC_OffOff = pow2db(bPxx_HOC_OffOff);

% Plot the power spectral density
figure;
%plot(freq, 10*log10(power));
plot(bFxx_HOC_OffOff, bPxxP_HOC_OffOff);
xlim([0 50])  % Limit x-axis to 50 Hz for better visualization
ylim([-10 25])  % Limit y-axis to 25 db for better visualization
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density');


%% bursting analysis, based on Torrecillos et al., 2018 (test w/ L_streamOfInt_s1 (Off Med, Off Stim @ 0 mA)) - HOC
% time-frq transform via complex Morlet wavelets (f0/σf = 7, f0 = 1-45 Hz, steps of 0.25 Hz).

% define 'lfp_data' as a LFP data matrix with dimensions [samples x trials]
lfp_data_HOC_OffOff = cleanLFP_HOC_OffOff;

fs = 250; % sampling rate
f0 = 1:0.25:45; % center frequency range from 1-45Hz in steps of 0.25Hz

% Frequency-Time Decomposition using complex Morlet wavelets
time = (1:length(lfp_data_HOC_OffOff)) / fs;
cwt_power = zeros(length(f0), length(lfp_data_HOC_OffOff));

for i = 1:length(f0)
    % create a sinusoidal wave modulated by a Gaussian window
    sinusoidal_wave = exp(2 * pi * 1i * f0(i) .* time); % complex cosine component of wavelet with an imaginary sine component
    Guassian_window = exp(-time.^2 / (2 * (f0(i) / 7)^2)); % tapers the sinusoid; (f0(i) / 7)^2 controls width of Gaussian
    wavelet = sinusoidal_wave .* Guassian_window; % applies wavelet across f0 range to create filters convolved with the LFP signal
    cwt_power(i, :) = abs(conv(lfp_data_HOC_OffOff, wavelet, 'same')).^2; % square to compute power
end

% time-domain bin size determined by 'VoicesPerOctave' (higher VPO = finer freq. res; lower time res?)
% continuous wavelet transform
figure;
cwt(lfp_data_HOC_OffOff,fs,'FrequencyLimits',[1 125],'VoicesPerOctave',14); % try 14, mults of 7
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
beta_power_time_courses_HOC_OffOff = mean(normalized_power(beta_range(peak_idx) + (-3:3), :), 1);

% Plot the beta power time course
figure('Renderer', 'opengl');
plot(time, beta_power_time_courses_HOC_OffOff);
xlim([0, round(max(time)-1)])
xlabel('Time (s)');
ylabel('Normalized Beta Power (µVp)');
title('LFP Beta Power')


%% Beta Burst Detection and Plotting for trimmed LFP - LSTN Off Med, Off Stim - HOC

% Compute mean beta power and st. dev beta power
mean_beta_power_HOC_OffOff = mean(beta_power_time_courses_HOC_OffOff);
std_beta_power_HOC_OffOff = std(beta_power_time_courses_HOC_OffOff);

% Set threshold at the 75th percentile of the mean beta power
threshold = prctile(mean_beta_power_HOC_OffOff, 75);

% Define beta bursts as time points exceeding the threshold for more than 3 oscillatory cycles (3 oscillatory cycles = ~100ms)
min_duration = round(0.1 * fs); % Minimum duration of a burst in samples (100 ms)
% min_duration = round(0.066 * fs);  % Reduce minimum duration to 2 osc. cycles(66 ms)

above_threshold = beta_power_time_courses_HOC_OffOff > threshold;
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
plot(time, beta_power_time_courses_HOC_OffOff);
xlim([0, round(max(time)-1)])
xlabel('Time (s)');
ylabel('Normalized Beta Power');
title('LFP Beta Power')

% Plot the beta power time course and beta burst overlays using xregion (test with L_streamOfInt_s1)
% https://www.mathworks.com/help/matlab/ref/xregion.html

figure('Renderer', 'opengl');
p = plot(time, beta_power_time_courses_HOC_OffOff); % store handle to the beta power plot
hold on;

% Overlay the detected bursts on the same plot using xregion
for i = 1:length(start_indices)
    burst_start_time = time(start_indices(i));
    burst_end_time = time(end_indices(i));
    xregion(burst_start_time, burst_end_time, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8]);
end

% Disable automatic legend updates to prevent xregion from adding entries
legend('auto update','off');

% Create a dummy patch handle 'h' for the detected bursts to add to the legend
h = patch(NaN, NaN, [0.8, 0.8, 0.8], 'FaceAlpha', 0.3);

% Add the beta power plot handle 'p' to the legend
legend([p, h], 'Beta Power', 'Detected Bursts');

% Plot
xlim([0, round(max(time)-1)])
xlabel('Time (s)');
ylabel('Normalized Beta Power');
title('Beta Bursting Dynamics');
hold off;





%%  % Isolate LFP and Video of iterest (LFP stream 9, DLC session 7) - LSTN Off Med, On Stim
% Kinematics and LFP - LSTN Off Med, On Stim

% On Stim
% js1_Row	sessDur(s)	DLC_sessID	DLC_sessDur(s) 	Hemisphere	Notes: case videos (10) for JSON Session Report 1 - Off Med
% 9	        36	        session007	    45	        L	        % • session007, set1 (Off Med, On Stim @ max mA)
% 10	    29	        session008		37          R	        % • session008, set1 (Off Med, On Stim @ max mA)

% L_streamOfInt_s7 (LFP stream 9) corresponds with DLC session007, (Off Med, On Stim @ max mA)

% Define LFP stream of interest:
L_streamOfInt_OffOn = L_streamOfInt_s7;

% Extract DLC video of interest
dlc_session_id_OffOn = eucINDICIES.videoID{9};

% Extract dlc Tablet start frame index
tablet_start_index_OffOn = eucINDICIES.StartInd(tablet_indices(3));
tablet_stop_index_OffOn = eucINDICIES.StopInd(tablet_indices(3));

% Load mat file corresponding to DLC video session
cd(DLC_mat_dir)
load('dlcDAT_20230912_session007_idea08_resnet50_rightCam-0000DLC.mat','outDATA')
dlc_labDAT_OffOn = outDATA;

% Trim mat file corresponding to DLC video session (based on offset b/t Video start/stop to Tablet start/stop index)
dlc_trimmed_OffOn = dlc_labDAT_OffOn(tablet_start_index_OffOn:tablet_stop_index_OffOn,:);


%% Synchronize corresponding LFP stream (trimmed by Tablet start index) s.t. LFP_start index syncs with Tablet_start index

% Sampling Rate Conversions and Calculations
totalNumSamples_Vid_OffOn = height(dlc_trimmed_OffOn);
totalNumSecs_OffOn = totalNumSamples_Vid_OffOn/60; % 60 fps
totalNumSamples_LFP_OffOn = floor(totalNumSecs_OffOn*250); % 250 samples per second

% Original signal (video) sampling rate at 60 Hz (fps)
ts_DLC_OffOn = 0:1/60:(height(dlc_trimmed_OffOn)-1)/60;

% Target signal (LFP) sampling rate at 250 Hz (samples per sec)
ts_LFP_OffOn = 0:1/250:(height(L_streamOfInt_OffOn)-1)/250;

% Interpolate - upsample kinematic data to match LFP sampling rate
allColNs = dlc_trimmed_OffOn.Properties.VariableNames;
dlc_trimmed_OffOn_int = table;
for coli = 1:width(dlc_trimmed_OffOn)

    % Tmp col
    tmpCol = allColNs{coli};

    % % Interpolation
    % x_250 = interp1(ts_DLC, transpose(dlc_lab2use.(tmpCol)), ts_LFP, 'spline');
    x_250 = interp1(ts_DLC_OffOn, dlc_trimmed_OffOn.(tmpCol), ts_LFP_OffOn);

    dlc_trimmed_OffOn_int.(tmpCol) = transpose(x_250);
end

trimFrame_int_OffOn = linspace(tablet_start_index_OffOn,tablet_stop_index_OffOn, length(ts_LFP_OffOn));


%% Kinematics and LFP plot (LFP stream 9, DLC session 7) - Off Med, On Stim

figure;
tiledlayout(2,1,"TileSpacing","tight")

xTime = ts_LFP_OffOn;

% kinematics
fTip1_X_smooth = smoothdata(dlc_trimmed_OffOn_int.fTip1_x,'gaussian',70);
nexttile
plot(xTime,fTip1_X_smooth);
xlim([0 round(max(ts_LFP_OffOn))])
xlabel('Time (s)')
ylabel('fTip1, X deflection (mm)')
title('Kinematics: Finger Tip X Position');

% LFP
ecg = perceive_ecg(transpose(L_streamOfInt_OffOn),250,0);
nexttile
plot(xTime,ecg.cleandata);
xlim([0 round(max(ts_LFP_OffOn))])
xlabel('Time (s)')
ylabel('LFP Amplitude (µV)');
title('Time Synced STN LFP Data'); 

% continuous wavelet transform
figure;
cwt(ecg.cleandata,fs,'FrequencyLimits',[1 125],'VoicesPerOctave',21); % try 14, mults of 7
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Continuous Wavelet Transform');


%% HOC segments (LFP stream 9, DLC session 7) - Off Med, On Stim - HOC

% Extract start and stop frame indices for 'Hand OC' (LFP stream 9, DLC session 7)
HOC_start_index_OffOn = eucINDICIES.StartInd(HandMov_indices(3));
HOC_stop_index_OffOn = eucINDICIES.StopInd(HandMov_indices(3));

% Extract the kinematic data within the 'Hand OC' start and stop frame indices 
dlc_trimmed_HOC_OffOn = dlc_labDAT_OffOn(HOC_start_index_OffOn:HOC_stop_index_OffOn,:);

% Smooth changes in fTip1_x position data (euclidian distances) within frame indices
fTip1_X_smooth_HOC_OffOn = smoothdata(dlc_trimmed_HOC_OffOn.fTip1_x, 'gaussian', 70);

% Normalize kinematic data
fTip1_X_normalized = fTip1_X_smooth_HOC_OffOn - min(fTip1_X_smooth_HOC_OffOn) + (0.01 * fTip1_X_smooth_HOC_OffOn);

% Time vector for kinematic data
xTime_Kinematics = linspace(0, length(fTip1_X_smooth_HOC_OffOn)/60, length(fTip1_X_smooth_HOC_OffOn));

% Plot normalized kinematic data over time
figure;
plot(xTime_Kinematics, fTip1_X_normalized);
xlim([0, round(max(xTime_Kinematics)-1)])
ylim([0, round(max(fTip1_X_normalized)+1)]) 
xlabel('Time (s)')
ylabel('Normalized fTip1, X deflection (mm)')
title('Kinematics: Finger Tip X Position');


%% Create synchronized Kinematics and LFP plot for HOC task (LFP stream 9, DLC session 7) - Off Med, On Stim - HOC
% 'trimFrame_int_OffOn' has the indices for the LFP data after being interpolated

% Find the closest matching indices in the interpolated kinematic data
[~, startInd_Kin] = min(abs(trimFrame_int_OffOn - HOC_start_index_OffOn));
[~, stopInd_Kin] = min(abs(trimFrame_int_OffOn - HOC_stop_index_OffOn));

% Trim the interpolated kinematic data to match the 'Hand OC' segment
dlc_trimmed_HOC_OffOn_int = dlc_trimmed_OffOn_int(startInd_Kin:stopInd_Kin, :);

% Find the corresponding indices in the LFP time vector
[~, startIND_LFP] = min(abs(ts_LFP_OffOn - ts_DLC_OffOn(HOC_start_index_OffOn)));
[~, endIND_LFP] = min(abs(ts_LFP_OffOn - ts_DLC_OffOn(HOC_stop_index_OffOn)));

% Extract the LFP data for Hand OC
LFP_HOC_OffOn = L_streamOfInt_OffOn(startIND_LFP:endIND_LFP);

% Calculate the corresponding time vector for the trimmed kinematic data
xTime_Kinematics_HOC = linspace(0, (stopInd_Kin - startInd_Kin + 1) / fs, stopInd_Kin - startInd_Kin + 1);

% Normalize the trimmed kinematic data so that the min is 0 and max is 50
fTip1_X_smooth_HOC = smoothdata(dlc_trimmed_HOC_OffOn_int.fTip1_x, 'gaussian', 70);
fTip1_X_normalized_HOC = (fTip1_X_smooth_HOC - min(fTip1_X_smooth_HOC)) * 25 / (max(fTip1_X_smooth_HOC) - min(fTip1_X_smooth_HOC));

% Create synchronized Kinematics and LFP plot for HOC task (Off Med, On Stim @ max mA)
figure;
tiledlayout(2,1,"TileSpacing","tight")

% Plot normalized kinematic data 
nexttile
plot(xTime_Kinematics_HOC, fTip1_X_normalized_HOC);
xlim([0, round(max(xTime_Kinematics_HOC)-1)])
ylim([0, round(max(fTip1_X_normalized_HOC)+1)]) 
xlabel('Time (s)');
ylabel('fTip1, X deflection (mm)');
title('Normalized Kinematics: Hand Movements');

% Plot LFP
xTime_LFP = linspace(0, length(LFP_HOC_OffOn)/250, length(LFP_HOC_OffOn));
ecg = perceive_ecg(LFP_HOC_OffOn', 250, 0);
nexttile
plot(xTime_LFP, ecg.cleandata);
xlim([0, round(max(xTime_LFP)-1)])
% ylim([round(min(ecg.cleandata)-1), round(max(ecg.cleandata)+1)])
ylim([-45, 26])
ylabel('LFP Amplitude (µV)');
title('Time Synced STN LFP Data');


%% Compute and Plot Power Spectral Density (PSD), (LFP stream 9, DLC session 7) - Off Med, On Stim - HOC

% Filter ecg
ts_LFP_HOC_OffOn = 0:1/250:(height(LFP_HOC_OffOn)-1)/250;
ecgClean_HOC_OffOn = perceive_ecg(transpose(LFP_HOC_OffOn),250,0);
cleanLFP_HOC_OffOn = ecgClean_HOC_OffOn.cleandata;

pspectrum(cleanLFP_HOC_OffOn,250,'FrequencyLimits',[0 50],'FrequencyResolution',3) % 2 or 3
[bPxx_HOC_OffOn, bFxx_HOC_OffOn] = pspectrum(cleanLFP_HOC_OffOn,250,'FrequencyLimits',[0 50],'FrequencyResolution',3); %
bPxxP_HOC_OffOn = pow2db(bPxx_HOC_OffOn);

% Plot the power spectral density
figure;
%plot(freq, 10*log10(power));
plot(bFxx_HOC_OffOn, bPxxP_HOC_OffOn);
xlim([0 50])  % Limit x-axis to 50 Hz for better visualization
ylim([-10 25])  % Limit y-axis to 25 db for better visualization
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density');


%% bursting analysis, based on Torrecillos et al., 2018 (LFP stream 9, DLC session 7) - Off Med, On Stim - HOC
% time-frq transform via complex Morlet wavelets (f0/σf = 7, f0 = 1-45 Hz, steps of 0.25 Hz).

% define 'lfp_data' as a LFP data matrix with dimensions [samples x trials]
lfp_data_HOC_OffOn = cleanLFP_HOC_OffOn;

fs = 250; % sampling rate
f0 = 1:0.25:45; % center frequency range from 1-45Hz in steps of 0.25Hz

% Frequency-Time Decomposition using complex Morlet wavelets
time = (1:length(lfp_data_HOC_OffOn)) / fs;
cwt_power = zeros(length(f0), length(lfp_data_HOC_OffOn));

for i = 1:length(f0)
    % create a sinusoidal wave modulated by a Gaussian window
    sinusoidal_wave = exp(2 * pi * 1i * f0(i) .* time); % complex cosine component of wavelet with an imaginary sine component
    Guassian_window = exp(-time.^2 / (2 * (f0(i) / 7)^2)); % tapers the sinusoid; (f0(i) / 7)^2 controls width of Gaussian
    wavelet = sinusoidal_wave .* Guassian_window; % applies wavelet across f0 range to create filters convolved with the LFP signal
    cwt_power(i, :) = abs(conv(lfp_data_HOC_OffOn, wavelet, 'same')).^2; % square to compute power
end

% time-domain bin size determined by 'VoicesPerOctave' (higher VPO = finer freq. res; lower time res?)
% continuous wavelet transform
figure;
cwt(lfp_data_HOC_OffOn,fs,'FrequencyLimits',[1 125],'VoicesPerOctave',14); % try 14, mults of 7
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
beta_power_time_courses_HOC_OffOn = mean(normalized_power(beta_range(peak_idx) + (-3:3), :), 1);

% Plot the beta power time course
figure('Renderer', 'opengl');
plot(time, beta_power_time_courses_HOC_OffOn);
xlim([0, round(max(time)-1)])
xlabel('Time (s)');
ylabel('Normalized Beta Power (µVp)');
title('LFP Beta Power')

%% Beta Burst Detection and Plotting for trimmed LFP - LSTN Off Med, On Stim - HOC

% Compute mean beta power
mean_beta_power_HOC_OffOn = mean(beta_power_time_courses_HOC_OffOn);
std_beta_power_HOC_OffOn = std(beta_power_time_courses_HOC_OffOn);

% Set threshold at the 75th percentile of the mean beta power
threshold = prctile(mean_beta_power_HOC_OffOn, 75);

% Define beta bursts as time points exceeding the threshold for more than 3 oscillatory cycles (3 oscillatory cycles = ~100ms)
min_duration = round(0.1 * fs); % Minimum duration of a burst in samples (100 ms)
% min_duration = round(0.066 * fs);  % Reduce minimum duration to 2 osc. cycles(66 ms)

above_threshold = beta_power_time_courses_HOC_OffOn > threshold;
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
plot(time, beta_power_time_courses_HOC_OffOn);
xlim([0, round(max(time)-1)])
xlabel('Time (s)');
ylabel('Normalized Beta Power');
title('LFP Beta Power')


% Plot the beta power time course and beta burst overlays using xregion (test with L_streamOfInt_s1)
% https://www.mathworks.com/help/matlab/ref/xregion.html

figure('Renderer', 'opengl');
p = plot(time, beta_power_time_courses_HOC_OffOn); % store handle to the beta power plot
hold on;

% Overlay the detected bursts on the same plot using xregion
for i = 1:length(start_indices)
    burst_start_time = time(start_indices(i));
    burst_end_time = time(end_indices(i));
    xregion(burst_start_time, burst_end_time, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8]);
end

% Disable automatic legend updates to prevent xregion from adding entries
legend('auto update','off');

% Create a dummy patch handle 'h' for the detected bursts to add to the legend
h = patch(NaN, NaN, [0.8, 0.8, 0.8], 'FaceAlpha', 0.3);

% Add the beta power plot handle 'p' to the legend
legend([p, h], 'Beta Power', 'Detected Bursts');

% Plot
xlim([0, round(max(time)-1)])
xlabel('Time (s)');
ylabel('Normalized Beta Power');
title('Beta Bursting Dynamics');
hold off;









%% DLC session IDs (sorted by condition)

videoIDs = eucINDICIES.videoID;

% OFF Med
DLC_s1 = videoIDs{1:4}; % set 1 [OFF Med, OFF Stim]
DLC_s3 = videoIDs{5:8}; % set 2 [OFF Med, OFF Stim]
% DLC_5_ramp % OFF Med
DLC_s7 = videoIDs{9:12}; % set 1 [OFF Med, ON Stim]
DLC_s9 = videoIDs{13:16}; % set 2 [OFF Med, ON Stim]

% ON Med
DLC_13 = videoIDs{17:20}; % set 1 [ON Med, OFF Stim]
DLC_15 = videoIDs{21:24}; % set 2 [ON Med, OFF Stim]
% DLC_18_ramp % ON Med
DLC_20 = videoIDs{25:28}; % set 1 [ON Med, ON Stim]
DLC_22 = videoIDs{29:32}; % set 2 [ON Med, ON Stim]

% DLC sess IDs:
L_DLC_sessIDs_OffMed = {DLC_s1, DLC_s3, DLC_s7, DLC_s9};
L_DLC_sessIDs_OnMed = {DLC_13, DLC_15, DLC_20, DLC_22};


%% Run PSD code on HOC_LFP for each LFP stream
% outer loop for LFP streams
% inner loop for movement indices

L_streamsofInt_OffMed_simple = {L_streamOfInt_s1, L_streamOfInt_s3, L_streamOfInt_s7, L_streamOfInt_s9};

% Loop through each LFP stream in Left STN
for lfp_i = 1:length(L_streamsofInt_OffMed_simple)
    LFP_streamOfInt_L = L_streamsofInt_OffMed_simple{lfp_i};
    DLC_sessionIDs_L = L_DLC_sessIDs_OffMed{lfp_i};

    % Convert LFP data for processing
    ts_LFP = linspace(0, length(LFP_streamOfInt_L)/250, length(LFP_streamOfInt_L));
    ecgClean_L = perceive_ecg(transpose(LFP_streamOfInt_L), 250, 0);
    cleanLFP_L = ecgClean_L.cleandata;

    % Ensure the session array is correctly formatted as a cell array
    currentSessions = L_DLC_sessIDs_OffMed{lfp_i};
    if ~iscell(currentSessions) % Check if not a cell, make it a cell
        currentSessions = {currentSessions};
    end

    % Loop through each DLC session ID corresponding to the LFP stream
    for dlc_i = 1:length(currentSessions)
        DLC_session_id = currentSessions{dlc_i};

        % Ensure DLC_session_id is a string if it's a cell
        if iscell(DLC_session_id)
            DLC_session_id = DLC_session_id{1};
        end

        % Loop through moveCSV files and load corresponding dlcDAT MAT files
        for csv_i = 1:length(moveCSV)

            tmpCSV = moveCSV{csv_i};

            % Split file names to extract relevant parts (dateID, sessID)
            nameParts = split(tmpCSV,'_');
            dateID = nameParts{1};
            sessID = nameParts{3};
            matName_title = [dateID , '-' , sessID];

            % Find and load corresponding dlcDAT MAT file
            matTempfind = [dateID , '_' , sessID];
            matInd = contains(mainMAT2 , matTempfind);
            matName = mainMAT2{matInd};
            cd(DLC_mat_dir)
            load(matName , 'outDATA')

            dlc_labDAT = outDATA;

            %for dlc_labDAT_i = 1:length(L_DLC_sessIDs_OffMed)
            % Find and load corresponding dlcDAT MAT files in 'moveCSV'
            % that match the filenames in 'L_DLC_sessIDs_OffMed' and
            % loop through those mat files

            % Determine indices for 'Hand OC' movement
            HOC_start_index = eucINDICIES.StartInd(HandMov_indices(dlc_i));
            HOC_stop_index = eucINDICIES.StopInd(HandMov_indices(dlc_i));

            % Interpolate and trim DLC data
            ts_DLC = linspace(0, (HOC_stop_index - HOC_start_index + 1) / 60, HOC_stop_index - HOC_start_index + 1);
            dlc_interp = interpolateDLCData(dlc_labDAT(HOC_start_index:HOC_stop_index, :), ts_DLC, ts_LFP);

            % Find LFP indices for the interpolated DLC data
            [~, startIND_LFP] = min(abs(ts_LFP - ts_DLC(1)));
            [~, endIND_LFP] = min(abs(ts_LFP - ts_DLC(end)));

            if endIND_LFP > length(ts_LFP)  % Ensure end index does not exceed time vector length
                endIND_LFP = length(ts_LFP);
            end

            HOC_LFP = cleanLFP_L(startIND_LFP:endIND_LFP);

            % Ensure the time vector matches the data length
            ts_trimmed = ts_LFP(startIND_LFP:endIND_LFP);

            % Adjust the interpolated data length
            if length(dlc_interp.fTip1_x) > length(ts_trimmed)
                dlc_interp.fTip1_x = dlc_interp.fTip1_x(1:length(ts_trimmed));
            elseif length(dlc_interp.fTip1_x) < length(ts_trimmed)
                % Handle cases where interpolated data is shorter than the time vector
                dlc_interp.fTip1_x = [dlc_interp.fTip1_x; nan(length(ts_trimmed) - length(dlc_interp.fTip1_x), 1)];
            end

            % Check lengths before plotting
            fprintf('Length of ts: %d, Length of data: %d\n', length(ts_trimmed), length(dlc_interp.fTip1_x));

            % Plotting
            figure;
            tiledlayout(2, 1, "TileSpacing", "tight");
            nexttile;
            plot(ts_trimmed, smoothdata(dlc_interp.fTip1_x, 'gaussian', 70));  % Correct data field used
            xlabel('Time (s)');
            ylabel('fTip1, X deflection');
            nexttile;
            plot(ts_trimmed, HOC_LFP);  % Ensure LFP data is also matched in length
            xlabel('Time (s)');
            ylabel('LFP (uV)');
            title(sprintf('Kinematics and LFP for Stream %d' , lfp_i));
        end
    end

end

















%% (?) Trim / normalize the LFP stream durations by motor trial indices condition

% baseline - NA

% Off Stim
% Hand O/C
% Pronation/Supination
% Elbow Flexion/Extension

% Ramp
% finger tap

% On Stim
% Hand O/C
% Pronation/Supination
% Elbow Flexion/Extension


%% compute LFP power / instantaneous LFP beta power and plot PSDs per session (using pspectrum function)

% pspectrum
% fxx - freq
% pxx - power

% [p, f, t] = pspectrum(streamOfInt, 'spectrogram');


% % bin data
% overlap
% spectrogram
% bp filter beta
% convert power
% hilbert transform
% time-frequency decomp.
% retain temporal resolution

%%

% positional data
% constext-dependent kinematic behav: resting, initiation, braking
% feature extraction and weighting via unsupervised ML
% movement vigor - rest vs. move - freq. of movement


%% Subfunctions

function [outTAB] = getDAYtimes(inputTIMES , inputSAMPLES) % JAT helper function

dayTIMES = cell(height(inputTIMES),1);
dayS = cell(height(inputTIMES),1);
fullDtime = cell(height(inputTIMES),1);
durations = zeros(height(inputTIMES),1);
streamOffsets = nan(height(inputTIMES),1);

for ti = 1:height(inputTIMES)
    inputTi = inputTIMES{ti};

    % The input date-time string
    % dateTimeStr = '2023-09-08T17:47:31.000Z';

    % Convert the string to a datetime object
    dateTimeObj = datetime(inputTi, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z');

    % Set the time zone to UTC (Coordinated Universal Time) since the input is in UTC
    dateTimeObj.TimeZone = 'UTC';

    % Convert to Mountain Time
    dateTimeObj_Mountain = datetime(dateTimeObj, 'TimeZone', 'America/Denver');

    % Extract the time component in AM/PM format
    timeComponent_AMPM = datetime(dateTimeObj_Mountain,'Format','hh:mm:ss a z');
    timeComponent_DATE = datetime(dateTimeObj_Mountain,'Format','dd-MMM-yyyy');

    % Display the time component
    % disp(['Time in Mountain Time (AM/PM): ', timeComponent_AMPM]);
    dayTIMES{ti} = timeComponent_AMPM;
    dayS{ti} = timeComponent_DATE;
    fullDtime{ti} = dateTimeObj_Mountain;
    durations(ti) = round(length(inputSAMPLES{ti})/250);
end

for di = 1:height(fullDtime)
    if di < height(fullDtime)
        t1 = fullDtime{di} + seconds(durations(di));
        t2 = fullDtime{di + 1};

        streamOffsets(di) = seconds(time(between(t1,t2)));
    end
end

dayT2 = cellfun(@(x) char(x), dayTIMES, 'UniformOutput',false);
dayS2 = cellfun(@(x) char(x), dayS, 'UniformOutput',false);
fDt = cellfun(@(x) char(x), fullDtime, 'UniformOutput',false);

outTAB = table(dayT2 , dayS2  , fDt , durations , streamOffsets,...
    'VariableNames',{'TimeOccur','DayOccur','FullNAT','Duration','Offset'});
end


function [outTAB] = getBSLFPtimes(inputTIMES) % JAT helper function

dayTIMES = cell(height(inputTIMES),1);
dayS = cell(height(inputTIMES),1);
fullDtime = cell(height(inputTIMES),1);

for ti = 1:height(inputTIMES)
    inputTi = inputTIMES{ti};

    % Convert the string to a datetime object
    dateTimeObj = datetime(inputTi, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z');

    % Set the time zone to UTC (Coordinated Universal Time) since the input is in UTC
    dateTimeObj.TimeZone = 'UTC';

    % Convert to Mountain Time
    dateTimeObj_Mountain = datetime(dateTimeObj, 'TimeZone', 'America/Denver');

    % Extract the time component in AM/PM format
    timeComponent_AMPM = datetime(dateTimeObj_Mountain,'Format','hh:mm:ss a z');
    timeComponent_DATE = datetime(dateTimeObj_Mountain,'Format','dd-MMM-yyyy');

    % Display the time component
    % disp(['Time in Mountain Time (AM/PM): ', timeComponent_AMPM]);
    dayTIMES{ti} = timeComponent_AMPM;
    dayS{ti} = timeComponent_DATE;
    fullDtime{ti} = dateTimeObj_Mountain;
end

dayT2 = cellfun(@(x) char(x), dayTIMES, 'UniformOutput',false);
dayS2 = cellfun(@(x) char(x), dayS, 'UniformOutput',false);
fDt = cellfun(@(x) char(x), fullDtime, 'UniformOutput',false);

outTAB = table(dayT2 , dayS2  , fDt,...
    'VariableNames',{'TimeOccur','DayOccur','FullNAT'});
end


function [inst_freq_filtered, ts_LFP_tmp] = HilbertTransform(lfpData, fs)
% Calculate time vector
ts_LFP_tmp = (0:length(lfpData)-1) / fs;

% Apply Hilbert transform
hilbert_eeg = hilbert(lfpData);
inst_phase = angle(hilbert_eeg);
inst_freq = diff(unwrap(inst_phase)) / (2 * pi * (1 / fs));

% Filter instantaneous frequency for 13-30 Hz
bandpass_filt = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', 13, 'HalfPowerFrequency2', 30, ...
    'SampleRate', fs);
inst_freq_filtered = filtfilt(bandpass_filt, inst_freq);
inst_freq_filtered(end + 1) = 0;  % Padding the last value
end


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
[p, f] = pspectrum(lfpData, 250, 'FrequencyLimits', [0 50], 'FrequencyResolution', 3);
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

