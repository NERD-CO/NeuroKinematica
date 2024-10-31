% Clinical LFP analyses - LSTN, On Med

%% Combine Percept LFP with DLC Video

% Data ingredients:
% 1) LFP data - JSON Session Reports (multiple rows (stream recordings) per report [metadata informs ID of row])
% 2) Movement data (.mat files) and Movement Indices (.csvs)

clear; close all; clc;

%% Directory set-up - Navigate b/t machines
pcname = getenv('COMPUTERNAME');

switch pcname
    case 'DESKTOP-I5CPDO7'   %%% JAT Desktop

        % mainDir = '';

    case 'DSKTP-JTLAB-EMR'   %%% ER Desktop

        mainDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\DLC_LFP';

    case 'NSG-M-FQBPFK3'     %%% ER PC

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
JSON_name1 = 'Report_Json_Session_Report_20230912T115956.json';
JSON_name2 = 'Report_Json_Session_Report_20230912T115939.json';

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
    sessDate_field_2 = currentJSON.SessionDate;

    % Convert the string to a datetime object
    dateTimeObj_2 = datetime(sessDate_field_2, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

    % Set the time zone to UTC (Coordinated Universal Time) since the input is in UTC
    dateTimeObj_2.TimeZone = 'UTC';

    % Convert to Mountain Time
    dateTimeObj_Mountain_2 = datetime(dateTimeObj_2, 'TimeZone', 'America/Denver');

    % Extract the time component in AM/PM format
    timeComponent_AMPM = datetime(dateTimeObj_Mountain_2,'Format','hh:mm:ss a z');
    timeComponent_DATE = datetime(dateTimeObj_Mountain_2,'Format','dd-MMM-yyyy');

    % convert timedomain to table
    BSTD_2 = currentJSON.BrainSenseTimeDomain; % struct
    BSTD_2_table = struct2table(BSTD_2);

    % plot raw and ecg-filtered time domain data for each row of current JSON file
    for BSTD_i = 1:size(BSTD_2_table, 1) % Loop through each row in current JSON

        % filter out ECG for each row of time domain data in current JSON
        tempData_1 = transpose(BSTD_2_table.TimeDomainData{BSTD_i}); % Transpose raw data for current row
        ecg = perceive_ecg(tempData_1, fs, plotit);

    end

    % navigate to time domain data
    [outTAB] = getBSLFPtimes(BSTD_2_table.FirstPacketDateTime);

    % Determine session start time - display first row of outTAB
    disp(outTAB(1,:));

    % Store the first row of outTAB in the session_StartTimes structure
    session_StartTimes.(sprintf('File%d', json_i)) = outTAB(1,:);

end

% JSON_name1 (...5956.json), FullNAT {'12-Sep-2023 10:17:12'},  Off Med
% JSON_name2 (...5939.json), FullNAT {'12-Sep-2023 11:31:25'},  On Med


%% Process ON Med JSON Session Report

% load JSON Session Report
js_1_name = 'Report_Json_Session_Report_20230912T115956.json' % Off Med
js_1 = jsondecode(fileread(js_1_name));

js_2_name = 'Report_Json_Session_Report_20230912T115939.json'; % On Med
js_2 = jsondecode(fileread(js_2_name));

% calculate the date/time
sessDate_field_2 = js_2.SessionDate;

% Convert the string to a datetime object
dateTimeObj_2 = datetime(sessDate_field_2, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

% Set the time zone to UTC (Coordinated Universal Time) since the input is in UTC
dateTimeObj_2.TimeZone = 'UTC';

% Convert to Mountain Time
dateTimeObj_Mountain_2 = datetime(dateTimeObj_2, 'TimeZone', 'America/Denver');

% convert BrainSense timedomain data to table
BSTD_2 = js_2.BrainSenseTimeDomain; % struct
BSTD_2_table = struct2table(BSTD_2); % table

% navigate to time domain data
BSLFP_times_2 = getBSLFPtimes(BSTD_2_table.FirstPacketDateTime);
BSLFP_times_2_table = cell2table(BSLFP_times_2.FullNAT);
disp(BSLFP_times_2(1,3)); % ensure correct file by start-time info


%% test - filter out ecg for first LFP stream

plot(BSTD_2_table.TimeDomainData{1}); % plot raw time domain data for row 1

streamOfInt_transposed = transpose(BSTD_2_table.TimeDomainData{1}); % transpose raw data for row 1
ecg = perceive_ecg(streamOfInt_transposed,fs,plotit); % run perceive_ecg function

hold on
plot(ecg.cleandata) % plot ecg-filtered time domain data for row 1
hold off


%% Load Left & Right STN BrainSense LFP streaming sessions in ON Med JSON Session Report

% BrainSense timedomain data table
streaming_TAB_2 = BSTD_2_table;

% Trim by STN of interest - LSTN, RSTN
stream_LEFT_1 = streaming_TAB_2(contains(streaming_TAB_2.Channel,'LEFT'),:); % L STN, R body
stream_RIGHT_1 = streaming_TAB_2(contains(streaming_TAB_2.Channel,'RIGHT'),:); % R STN, L body

% Determine duration (in seconds) of each stream
stream_LEFT_1_times = getDAYtimes(stream_LEFT_1.FirstPacketDateTime, stream_LEFT_1.TimeDomainData);
LEFT_sessDurations_seconds_1 = stream_LEFT_1_times.Duration;

stream_RIGHT_1_times = getDAYtimes(stream_RIGHT_1.FirstPacketDateTime, stream_RIGHT_1.TimeDomainData);
RIGHT_sessDurations_seconds_1 = stream_RIGHT_1_times.Duration;


%% Load Left & Right STN BrainSense LFP streaming sessions in ON Med JSON Session Report

% BrainSense timedomain data table
streaming_TAB_2 = BSTD_2_table;

% Trim by STN of interest - LSTN, RSTN
stream_LEFT_2 = streaming_TAB_2(contains(streaming_TAB_2.Channel,'LEFT'),:); % L STN, R body
stream_RIGHT_2 = streaming_TAB_2(contains(streaming_TAB_2.Channel,'RIGHT'),:); % R STN, L body

% Determine duration (in seconds) of each stream
stream_LEFT_2_times = getDAYtimes(stream_LEFT_2.FirstPacketDateTime, stream_LEFT_2.TimeDomainData);
LEFT_sessDurations_seconds_2 = stream_LEFT_2_times.Duration;

stream_RIGHT_2_times = getDAYtimes(stream_RIGHT_2.FirstPacketDateTime, stream_RIGHT_2.TimeDomainData);
RIGHT_sessDurations_seconds_2 = stream_RIGHT_2_times.Duration;


%% Separate LFP streams and corresponding DLC sessions of interest by condition
% 9/12/2023 BrainSense LFP streams (14) in JSON Session Report 2 - On Med						

% baseline
% js2_Row	sessDur(s)	DLC_sessID	DLC_sessDur(s) 	Hemisphere	Notes
% 1	        32	        session011	    45	        L	        % • baseline_L, (On Med, Off Stim @ 0 mA)
% 2	        34	        session012	    42	        R	        % • baseline_R, (On Med, Off Stim @ 0 mA)

% Off Stim
% js2_Row	sessDur(s)	DLC_sessID	DLC_sessDur(s) 	Hemisphere	Notes
% 3	        38	        session013	    45	        L	        % • session013, set1 (On Med, Off Stim @ 0 mA) 
% 4	        29	        session014	    42	        R	        % • session014, set1 (On Med, Off Stim @ 0 mA) 
% 5	        32	        session015	    41	        L	        % • session015, set2 (On Med, Off Stim @ 0 mA)
    % 6                                                             % •  session016 interrupted by RC pt. compensation dropoff
% 7	        29	        session017	    37	        R	        % • session017, set2 (On Med, Off Stim @ 0 mA)

% Ramp
% js2_Row	sessDur(s)	DLC_sessID	DLC_sessDur(s) 	Hemisphere	Notes
% 8	        153	        session018	    189	        L	        % • session018 (On Med, Stim Ramping) 
    % 9                                                             % • session019 (On Med, Stim Ramping) - LFP not captured
    % 10                                                            % •  scrap 

% On Stim
% js2_Row	sessDur(s)	DLC_sessID	DLC_sessDur(s) 	Hemisphere	Notes   
% 11	    28	        session020	    36	        L	        % • session020, set1 (On Med, On Stim @ max mA)
% 12	    26	        session021	    34	        R	        % • session021, set1 (On Med, On Stim @ max mA)
% 13	    26	        session022	    34	        L	        % • session022, set2 (On Med, On Stim @ max mA)
% 14	    25	        session023	    33	        R	        % • session023, set2 (On Med, On Stim @ max mA)


%% Isolate specific rows of interest in JSON Session Report 1 - On Med

% Left STN, Right body
% Specific Rows of interest
L_rowfromTab_js2_bsln = 1;
L_streamOfInt_js2_bsln = stream_LEFT_2.TimeDomainData{L_rowfromTab_js2_bsln}; % L_baseline (On Med, Off Stim @ 0 mA)

L_rowfromTab_js2_s1_OffStim = 3;
L_streamOfInt_js2_s1_OffStim = stream_LEFT_2.TimeDomainData{L_rowfromTab_js2_s1_OffStim}; % L_set1 (On Med, Off Stim @ 0 mA)

L_rowfromTab_js2_s2_OffStim = 5;
L_streamOfInt_js2_s2_OffStim = stream_LEFT_2.TimeDomainData{L_rowfromTab_js2_s2_OffStim}; % L_set2 (On Med, Off Stim @ 0 mA)

L_rowfromTab_js2_ramp = 8;
L_streamOfInt_js2_ramp = stream_LEFT_2.TimeDomainData{L_rowfromTab_js2_ramp}; % L_ramp (On Med, Stim Ramping) 

L_rowfromTab_js2_s1_OnStim = 11;
L_streamOfInt_js2_s1_OnStim = stream_LEFT_2.TimeDomainData{L_rowfromTab_js2_s1_OnStim}; % L_set1 (On Med, On Stim @ max mA)

L_rowfromTab_js2_s2_OnStim = 13;
L_streamOfInt_js2_s2_OnStim = stream_LEFT_2.TimeDomainData{L_rowfromTab_js2_s2_OnStim}; % L_set2 (On Med, On Stim @ max mA)


% Right STN, Left body
% Specific Rows of interest
R_rowfromTab_js2_bsln = 2;
R_streamOfInt_js2_bsln = stream_RIGHT_2.TimeDomainData{R_rowfromTab_js2_bsln}; % R_baseline, (On Med, Off Stim @ 0 mA)

R_rowfromTab_js2_s1_OffStim = 4;
R_streamOfInt_js2_s1_OffStim = stream_RIGHT_2.TimeDomainData{R_rowfromTab_js2_s1_OffStim}; % R_set1 (On Med, Off Stim @ 0 mA)

R_rowfromTab_js2_s2_OffStim = 7;
R_streamOfInt_js2_s2_OffStim = stream_RIGHT_2.TimeDomainData{R_rowfromTab_js2_s2_OffStim}; % R_set2 (On Med, Off Stim @ 0 mA)

% R_rowfromTab_js2_ramp = 9;
% R_streamOfInt_js2_ramp = stream_RIGHT_2.TimeDomainData{R_rowfromTab_js2_ramp}; % R_ramp (On Med, Stim Ramping)

R_rowfromTab_js2_s1_OnStim = 12;
R_streamOfInt_js2_s1_OnStim = stream_RIGHT_2.TimeDomainData{R_rowfromTab_js2_s1_OnStim}; % R_set1 (On Med, On Stim @ max mA)

R_rowfromTab_js2_s2_OnStim = 14;
R_streamOfInt_js2_s2_OnStim = stream_RIGHT_2.TimeDomainData{R_rowfromTab_js2_s2_OnStim}; % R_set2 (On Med, On Stim @ max mA)


%% Sanity Check Plotting

% L_streamsofInt_OnMed = {L_streamOfInt_js2_bsln, L_streamOfInt_js2_s1_OffStim, L_streamOfInt_js2_s2_OffStim, L_streamOfInt_js2_ramp, L_streamOfInt_js2_s1_OnStim, L_streamOfInt_js2_s2_OnStim}; % w/ baseline
L_streamsofInt_OnMed = {L_streamOfInt_js2_s1_OffStim, L_streamOfInt_js2_s2_OffStim, L_streamOfInt_js2_s1_OnStim, L_streamOfInt_js2_s2_OnStim}; % w/o baseline and ramp

plotTitles = {
    'LSTN set1 (Off On, Off Stim @ 0 mA)', 
    'LSTN set2 (Off On, Off Stim @ 0 mA)', 
    
    'LSTN set1 (Off On, On Stim @ max mA)', 
    'LSTN set2 (Off On, On Stim @ max mA)'};

figure
for L_i = 1:length(L_streamsofInt_OnMed)
    subplot(4,1,L_i)
    plot(L_streamsofInt_OnMed{L_i})
    title(plotTitles{L_i}) 
    xlabel('Time (s)') 
    ylabel('LFP Amplitude (uV)') 
    grid on 
end


% R_streamsofInt_OnMed = {R_streamOfInt_js2_bsln, R_streamOfInt_js2_s1_OffStim, R_streamOfInt_js2_s2_OffStim, R_streamOfInt_js2_ramp, R_streamOfInt_js2_s1_OnStim, R_streamOfInt_js2_s2_OnStim}; % w/ baseline
R_streamsofInt_OnMed = {R_streamOfInt_js2_s1_OffStim, R_streamOfInt_js2_s2_OffStim, R_streamOfInt_js2_s1_OnStim, R_streamOfInt_js2_s2_OnStim}; % w/o baseline and ramp

plotTitles = {
    'RSTN set1 (On Med, Off Stim @ 0 mA)', 
    'RSTN set2 (On Med, Off Stim @ 0 mA)', 
    
    'RSTN set1 (On Med, On Stim @ max mA)', 
    'RSTN set2 (On Med, On Stim @ max mA)'};

figure
for R_i = 1:length(R_streamsofInt_OnMed)
    subplot(4,1,R_i)
    plot(R_streamsofInt_OnMed{R_i})
    title(plotTitles{R_i}) 
    xlabel('Time (s)') 
    ylabel('LFP Amplitude (uV)') 
    grid on 
end



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


%% Filter, Hilbert Transform (On Med, Off Stim @ 0 mA))

% L_streamOfInt_js2_s1_OffStim (LFP stream 3) corresponds with DLC session013 (movement video), set1 (On Med, Off Stim @ 0 mA)

% Define LFP stream of interest:
L_streamOfInt_OnOff = L_streamOfInt_js2_s1_OffStim;

fs = 250;
t = 1/250;

% Target signal (LFP) sampling rate at 250 Hz (samples per sec)
ts_LFP = 0:1/250:(height(L_streamOfInt_OnOff)-1)/250;

ecgClean = perceive_ecg(transpose(L_streamOfInt_OnOff),250,0);

cleanLFP_OnOff = ecgClean.cleandata;

% Calculate instantaneous phase and freq. using Hilbert transform
hilbert = hilbert(cleanLFP_OnOff);
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
LFP_notchfilt = filtfilt(b_notch, a_notch, cleanLFP_OnOff); % Apply notch filter to remove 60Hz line noise

% Comb Filter centered at 60Hz and its harmonics (120Hz, 180Hz, 240Hz, 300Hz, 360Hz, etc.) 
bw = f_notch / Q;  % Bandwidth, bw = (fo/(fs/2))/q;
n_harmonics = floor((fs/2) / f_notch);  % Number of harmonics within Nyquist frequency
[b_comb, a_comb] = iircomb(floor(fs/f_notch), bw/(fs/2), 'notch');  % Design comb filter
LFP_combfilt = filtfilt(b_comb, a_comb, cleanLFP_OnOff); % Apply comb filter to remove harmonics


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

%% Compute and Plot Power Spectral Density (PSD), (On Med, Off Stim @ 0 mA))
% Use pspectrum to compute the power spectral density
% [power, freq] = pspectrum(cleanLFP_s1, fs);

pspectrum(cleanLFP_OnOff,250,'FrequencyLimits',[0 50],'FrequencyResolution',3) % 2 or 3
[bPxx,bFxx] = pspectrum(cleanLFP_OnOff,250,'FrequencyLimits',[0 50],'FrequencyResolution',3); %
bPxxP = pow2db(bPxx);

% Plot the power spectral density
figure;
%plot(freq, 10*log10(power));
plot(bFxx, bPxxP);
xlim([0 50])  % Limit x-axis to 50 Hz for better visualization
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density of LFP (untrimmed)');


%% bursting analysis, based on Torrecillos et al., 2018 (On Med, Off Stim @ 0 mA)
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
title('Continuous Wavelet Transform (untrimmed)');

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

%% Trim beta power lfp and time data based on known/computed offset (On Med, Off Stim @ 0 mA))

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

%% continuous wavelet transform on trimmed LFP (On Med, Off Stim @ 0 mA))

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

%% Beta Burst Detection and Plotting for trimmed LFP (On Med, Off Stim @ 0 mA))

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

%% Compute and Plot PSD on trimmed LFP (On Med, Off Stim @ 0 mA))

pspectrum(lfp_data_L_trimmed,250,'FrequencyLimits',[0 50],'FrequencyResolution',3) % 2 or 3
[bPxx,bFxx] = pspectrum(lfp_data_L_trimmed,250,'FrequencyLimits',[0 50],'FrequencyResolution',3); %
bPxxP = pow2db(bPxx);

% Plot the power spectral density
figure;
%plot(freq, 10*log10(power));
plot(bFxx, bPxxP);
xlim([0 50])  % Limit x-axis to 50 Hz for better visualization
ylim([-10 25])  % Limit y-axis to 25 db for better visualization
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density of LFP');


%% Loop through each LFP stream in Left STN - Apply Filtering, PSD, CWT, and Beta Bursting code to all LFP streams (On Med, Off Stim @ 0 mA))
for i = 1:length(L_streamsofInt_OnMed)

    lfp_data_L = L_streamsofInt_OnMed{i};

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






%%  % Isolate LFP and Video of iterest (LFP stream 3, DLC session 13) - On Med, Off Stim
% Kinematics and LFP - On Med, Off Stim

% Off Stim
% js2_Row	sessDur(s)	DLC_sessID	DLC_sessDur(s) 	Hemisphere	Notes
% 3	        38	        session013	    45	        L	        % • session013, set1 (On Med, Off Stim @ 0 mA) 
% 4	        29	        session014	    42	        R	        % • session014, set1 (On Med, Off Stim @ 0 mA) 

% L_streamOfInt_js2_s1_OffStim (LFP stream 3) corresponds with DLC session013 (movement video), set1 (On Med, Off Stim @ 0 mA)

% Define LFP stream of interest:
L_streamOfInt_OnOff = L_streamOfInt_js2_s1_OffStim;

% Extract DLC video of interest
dlc_session_id_OnOff = eucINDICIES.videoID{18};

% Load mat file corresponding to DLC video session
cd(DLC_mat_dir)
load('dlcDAT_20230912_session013_idea08_resnet50_rightCam-0000DLC.mat','outDATA')
dlc_labDAT_OnOff = outDATA;

% Extract dlc Tablet start and stop frame index
tablet_start_index_OnOff = eucINDICIES.StartInd(tablet_indices(5));
tablet_stop_index_OnOff = eucINDICIES.StopInd(tablet_indices(5));

% Trim mat file corresponding to DLC video session (based on offset b/t Video start/stop to Tablet start/stop index)
dlc_trimmed_OnOff = dlc_labDAT_OnOff(tablet_start_index_OnOff:tablet_stop_index_OnOff,:); % 2821x40 table

%% Synchronize corresponding LFP stream (trimmed by Tablet start index) s.t. LFP_start index syncs with Tablet_start index (On Med, Off Stim @ 0 mA))

% Sampling Rate Conversions and Calculations
totalNumSamples_Vid_OnOff = height(dlc_trimmed_OnOff);
totalNumSecs_OnOff = totalNumSamples_Vid_OnOff/60; % 60 fps
totalNumSamples_LFP_OnOff = floor(totalNumSecs_OnOff*250); % 250 samples per second

% Original signal (video) sampling rate at 60 Hz (fps)
ts_DLC_OnOff = 0:1/60:(height(dlc_trimmed_OnOff)-1)/60;

% Target signal (LFP) sampling rate at 250 Hz (samples per sec)
ts_LFP_OnOff = 0:1/250:(height(L_streamOfInt_OnOff)-1)/250;

% Interpolate - upsample kinematic data to match LFP sampling rate
allColNs = dlc_trimmed_OnOff.Properties.VariableNames;
dlc_trimmed_OnOff_int = table;
for coli = 1:width(dlc_trimmed_OnOff)

    % Tmp col
    tmpCol = allColNs{coli};

    % % Interpolation
    % x_250 = interp1(ts_DLC, transpose(dlc_lab2use.(tmpCol)), ts_LFP, 'spline');
    x_250 = interp1(ts_DLC_OnOff, dlc_trimmed_OnOff.(tmpCol), ts_LFP_OnOff);

    dlc_trimmed_OnOff_int.(tmpCol) = transpose(x_250);
end

trimFrame_int_OnOff = linspace(tablet_start_index_OnOff,tablet_stop_index_OnOff, length(ts_LFP_OnOff));


% Kinematics and LFP plot (On Med, Off Stim @ 0 mA))
figure;
tiledlayout(2,1,"TileSpacing","tight")

xTime = ts_LFP_OnOff;

% kinematics
fTip1_X_smooth = smoothdata(dlc_trimmed_OnOff_int.fTip1_x,'gaussian',70);
nexttile
plot(xTime,fTip1_X_smooth);
xlim([0 round(max(ts_LFP_OnOff))])
xlabel('Time (s)')
ylabel('fTip1, X deflection (mm)')
title('Kinematics: Finger Tip X Position');

% LFP
ecg = perceive_ecg(transpose(L_streamOfInt_OnOff),250,0);
nexttile
plot(xTime,ecg.cleandata);
xlim([0 round(max(ts_LFP_OnOff))])
xlabel('Time (s)')
ylabel('LFP Amplitude (µV)');
title('Time Synced STN LFP Data'); 

% continuous wavelet transform
figure;
cwt(ecg.cleandata,fs,'FrequencyLimits',[1 125],'VoicesPerOctave',21); % try 14, mults of 7
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Continuous Wavelet Transform');


%% Hand Open/Close segment, LSTN On Med, Off Stim - HOC

% Extract start and stop frame indices for 'Hand OC' in DLC session data
HOC_start_index_OnOff = eucINDICIES.StartInd(HandMov_indices(5));
HOC_stop_index_OnOff = eucINDICIES.StopInd(HandMov_indices(5));

% Extract the kinematic data within the 'Hand OC' start and stop frame indices 
dlc_trimmed_HOC_OnOff = dlc_labDAT_OnOff(HOC_start_index_OnOff:HOC_stop_index_OnOff,:);

% Smooth changes in fTip1_x position data (euclidian distances) within frame indices
fTip1_X_smooth_HOC_OnOff = smoothdata(dlc_trimmed_HOC_OnOff.fTip1_x, 'gaussian', 70);

% Normalize kinematic data
fTip1_X_normalized = fTip1_X_smooth_HOC_OnOff - min(fTip1_X_smooth_HOC_OnOff) + (0.01 * fTip1_X_smooth_HOC_OnOff);

% Time vector for kinematic data
xTime_Kinematics = linspace(0, length(fTip1_X_smooth_HOC_OnOff)/60, length(fTip1_X_smooth_HOC_OnOff));

% Plot normalized kinematic data over time
figure;
plot(xTime_Kinematics, fTip1_X_normalized);
xlim([0, round(max(xTime_Kinematics)-1)])
ylim([0, round(max(fTip1_X_normalized)+1)]) 
xlabel('Time (s)')
ylabel('Normalized fTip1, X deflection (mm)')
title('Kinematics: Finger Tip X Position');

%% Create synchronized Kinematics and LFP plot for HOC task - On Med, Off Stim - HOC
% 'trimFrame_int_OffOff' has the indices for the LFP data after being interpolated

% Find the closest matching indices in the interpolated kinematic data
[~, startInd_Kin] = min(abs(trimFrame_int_OnOff - HOC_start_index_OnOff));
[~, stopInd_Kin] = min(abs(trimFrame_int_OnOff - HOC_stop_index_OnOff));

% Trim the interpolated kinematic data to match the 'Hand OC' segment
dlc_trimmed_HOC_OnOff_int = dlc_trimmed_OnOff_int(startInd_Kin:stopInd_Kin, :);

% Find the corresponding indices in the LFP time vector
[~, startIND_LFP] = min(abs(ts_LFP_OnOff - ts_DLC_OnOff(HOC_start_index_OnOff)));
[~, endIND_LFP] = min(abs(ts_LFP_OnOff - ts_DLC_OnOff(HOC_stop_index_OnOff)));

% Extract the LFP data for Hand OC
LFP_HOC_OnOff = L_streamOfInt_OnOff(startIND_LFP:endIND_LFP);

% Calculate the corresponding time vector for the trimmed kinematic data
xTime_Kinematics_HOC = linspace(0, (stopInd_Kin - startInd_Kin + 1) / fs, stopInd_Kin - startInd_Kin + 1);

% Normalize the trimmed kinematic data so that the min is 0 and max is 50
fTip1_X_smooth_HOC = smoothdata(dlc_trimmed_HOC_OnOff_int.fTip1_x, 'gaussian', 70);
fTip1_X_normalized_HOC = (fTip1_X_smooth_HOC - min(fTip1_X_smooth_HOC)) * 25 / (max(fTip1_X_smooth_HOC) - min(fTip1_X_smooth_HOC));

% HOC Kinematics and LFP plot (LFP stream 3, DLC session 1) - Off Med, Off Stim
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
xTime_LFP = linspace(0, length(LFP_HOC_OnOff)/250, length(LFP_HOC_OnOff));
ecg = perceive_ecg(LFP_HOC_OnOff', 250, 0);
nexttile
plot(xTime_LFP, ecg.cleandata);
xlim([0, round(max(xTime_LFP)-1)])
% ylim([round(min(ecg.cleandata)-1), round(max(ecg.cleandata)+1)])
ylim([-50, 25])
ylabel('LFP Amplitude (µV)');
title('Time Synced STN LFP Data');


%% Compute and Plot Power Spectral Density (PSD), (On Med, Off Stim @ 0 mA)) - HOC

% Filter ecg
ts_LFP_HOC_OnOff = 0:1/250:(height(LFP_HOC_OnOff)-1)/250;
ecgClean_HOC_OnOff = perceive_ecg(transpose(LFP_HOC_OnOff),250,0);
cleanLFP_HOC_OnOff = ecgClean_HOC_OnOff.cleandata;

pspectrum(cleanLFP_HOC_OnOff,250,'FrequencyLimits',[0 50],'FrequencyResolution',3) % 2 or 3
[bPxx_HOC_OnOff, bFxx_HOC_OnOff] = pspectrum(cleanLFP_HOC_OnOff,250,'FrequencyLimits',[0 50],'FrequencyResolution',3); %
bPxxP_HOC_OnOff = pow2db(bPxx_HOC_OnOff);

% Plot the power spectral density
figure;
%plot(freq, 10*log10(power));
plot(bFxx_HOC_OnOff, bPxxP_HOC_OnOff);
xlim([0 50])  % Limit x-axis to 50 Hz for better visualization
ylim([-10 25])  % Limit y-axis to 25 db for better visualization
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density');

%% bursting analysis, based on Torrecillos et al., 2018 (On Med, Off Stim @ 0 mA)) - HOC
% time-frq transform via complex Morlet wavelets (f0/σf = 7, f0 = 1-45 Hz, steps of 0.25 Hz).

% define 'lfp_data' as a LFP data matrix with dimensions [samples x trials]
lfp_data_HOC_OnOff = cleanLFP_HOC_OnOff;

fs = 250; % sampling rate
f0 = 1:0.25:45; % center frequency range from 1-45Hz in steps of 0.25Hz

% Frequency-Time Decomposition using complex Morlet wavelets
time = (1:length(lfp_data_HOC_OnOff)) / fs;
cwt_power = zeros(length(f0), length(lfp_data_HOC_OnOff));

for i = 1:length(f0)
    % create a sinusoidal wave modulated by a Gaussian window
    sinusoidal_wave = exp(2 * pi * 1i * f0(i) .* time); % complex cosine component of wavelet with an imaginary sine component
    Guassian_window = exp(-time.^2 / (2 * (f0(i) / 7)^2)); % tapers the sinusoid; (f0(i) / 7)^2 controls width of Gaussian
    wavelet = sinusoidal_wave .* Guassian_window; % applies wavelet across f0 range to create filters convolved with the LFP signal
    cwt_power(i, :) = abs(conv(lfp_data_HOC_OnOff, wavelet, 'same')).^2; % square to compute power
end

% time-domain bin size determined by 'VoicesPerOctave' (higher VPO = finer freq. res; lower time res?)
% continuous wavelet transform
figure;
cwt(lfp_data_HOC_OnOff,fs,'FrequencyLimits',[1 125],'VoicesPerOctave',14); % try 14, mults of 7
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
beta_power_time_courses_HOC_OnOff = mean(normalized_power(beta_range(peak_idx) + (-3:3), :), 1);

% Plot the beta power time course
figure('Renderer', 'opengl');
plot(time, beta_power_time_courses_HOC_OnOff);
xlim([0, round(max(time)-1)])
xlabel('Time (s)');
ylabel('Normalized Beta Power (µVp)');
title('LFP Beta Power')

%% Beta Burst Detection and Plotting for trimmed LFP - HOC On Med, Off Stim

% Compute mean beta power
mean_beta_power_HOC_OnOff = mean(beta_power_time_courses_HOC_OnOff);
std_beta_power_HOC_OnOff = std(beta_power_time_courses_HOC_OnOff);

% Set threshold at the 75th percentile of the mean beta power
threshold = prctile(mean_beta_power_HOC_OnOff, 75);

% Define beta bursts as time points exceeding the threshold for more than 3 oscillatory cycles (3 oscillatory cycles = ~100ms)
min_duration = round(0.1 * fs); % Minimum duration of a burst in samples (100 ms)
% min_duration = round(0.066 * fs);  % Reduce minimum duration to 2 osc. cycles(66 ms)

above_threshold = beta_power_time_courses_HOC_OnOff > threshold;
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
plot(time, beta_power_time_courses_HOC_OnOff);
xlim([0, round(max(time)-1)])
xlabel('Time (s)');
ylabel('Normalized Beta Power');
title('LFP Beta Power')


% Plot the beta power time course and beta burst overlays using xregion (test with L_streamOfInt_s1)
% https://www.mathworks.com/help/matlab/ref/xregion.html

figure('Renderer', 'opengl');
p = plot(time, beta_power_time_courses_HOC_OnOff); % store handle to the beta power plot
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







%%  % Isolate LFP and Video of iterest (LFP stream 11, DLC session 20) - LSTN On Med, On Stim
% Kinematics and LFP - LSTN On Med, On Stim

% On Stim
% js2_Row	sessDur(s)	DLC_sessID	DLC_sessDur(s) 	Hemisphere	Notes   
% 11	    28	        session020	    36	        L	        % • session020, set1 (On Med, On Stim @ max mA)
% 12	    26	        session021	    34	        R	        % • session021, set1 (On Med, On Stim @ max mA)

% Define LFP stream of interest:
L_streamOfInt_OnOn = L_streamOfInt_js2_s2_OnStim;

% Extract DLC video of interest
dlc_session_id_OnOn = eucINDICIES.videoID{30};

% Extract dlc Tablet start frame index
tablet_start_index_OnOn = eucINDICIES.StartInd(tablet_indices(8));
tablet_stop_index_OnOn = eucINDICIES.StopInd(tablet_indices(8));

% Load mat file corresponding to DLC video session
cd(DLC_mat_dir)
load('dlcDAT_20230912_session022_idea08_resnet50_rightCam-0000DLC.mat','outDATA')
dlc_labDAT_OnOn = outDATA;

% Trim mat file corresponding to DLC video session (based on offset b/t Video start/stop to Tablet start/stop index)
dlc_trimmed_OnOn = dlc_labDAT_OnOn(tablet_start_index_OnOn:tablet_stop_index_OnOn,:);


%% Synchronize corresponding LFP stream (trimmed by Tablet start index) s.t. LFP_start index syncs with Tablet_start index

% Sampling Rate Conversions and Calculations
totalNumSamples_Vid_OnOn = height(dlc_trimmed_OnOn);
totalNumSecs_OnOn = totalNumSamples_Vid_OnOn/60; % 60 fps
totalNumSamples_LFP_OnOn = floor(totalNumSecs_OnOn*250); % 250 samples per second

% Original signal (video) sampling rate at 60 Hz (fps)
ts_DLC_OnOn = 0:1/60:(height(dlc_trimmed_OnOn)-1)/60;

% Target signal (LFP) sampling rate at 250 Hz (samples per sec)
ts_LFP_OnOn = 0:1/250:(height(L_streamOfInt_OnOn)-1)/250;

% Interpolate - upsample kinematic data to match LFP sampling rate
allColNs = dlc_trimmed_OnOn.Properties.VariableNames;
dlc_trimmed_OnOn_int = table;
for coli = 1:width(dlc_trimmed_OnOn)

    % Tmp col
    tmpCol = allColNs{coli};

    % % Interpolation
    % x_250 = interp1(ts_DLC, transpose(dlc_lab2use.(tmpCol)), ts_LFP, 'spline');
    x_250 = interp1(ts_DLC_OnOn, dlc_trimmed_OnOn.(tmpCol), ts_LFP_OnOn);

    dlc_trimmed_OnOn_int.(tmpCol) = transpose(x_250);
end

trimFrame_int_OnOn = linspace(tablet_start_index_OnOn,tablet_stop_index_OnOn, length(ts_LFP_OnOn));


%% Kinematics and LFP plot - On Med, On Stim

figure;
tiledlayout(2,1,"TileSpacing","tight")

xTime = ts_LFP_OnOn;

% kinematics
fTip1_X_smooth = smoothdata(dlc_trimmed_OnOn_int.fTip1_x,'gaussian',70);
nexttile
plot(xTime,fTip1_X_smooth);
xlim([0 round(max(ts_LFP_OnOn))])
xlabel('Time (s)')
ylabel('fTip1, X deflection (mm)')
title('Kinematics: Finger Tip X Position');

% LFP
ecg = perceive_ecg(transpose(L_streamOfInt_OnOn),250,0);
nexttile
plot(xTime,ecg.cleandata);
xlim([0 round(max(ts_LFP_OnOn))])
xlabel('Time (s)')
ylabel('LFP Amplitude (µV)');
title('Time Synced STN LFP Data'); 

% continuous wavelet transform
figure;
cwt(ecg.cleandata,fs,'FrequencyLimits',[1 125],'VoicesPerOctave',21); % try 14, mults of 7
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Continuous Wavelet Transform');


%% HOC segments - On Med, On Stim

% Extract start and stop frame indices for 'Hand OC' (LFP stream 9, DLC session 7)
HOC_start_index_OnOn = eucINDICIES.StartInd(HandMov_indices(8));
HOC_stop_index_OnOn = eucINDICIES.StopInd(HandMov_indices(8));

% Extract the kinematic data within the 'Hand OC' start and stop frame indices 
dlc_trimmed_HOC_OnOn = dlc_labDAT_OnOn(HOC_start_index_OnOn:HOC_stop_index_OnOn,:);

% Smooth changes in fTip1_x position data (euclidian distances) within frame indices
fTip1_X_smooth_HOC_OnOn = smoothdata(dlc_trimmed_HOC_OnOn.fTip1_x, 'gaussian', 70);

% Normalize kinematic data
fTip1_X_normalized = fTip1_X_smooth_HOC_OnOn - min(fTip1_X_smooth_HOC_OnOn) + (0.01 * fTip1_X_smooth_HOC_OnOn);

% Time vector for kinematic data
xTime_Kinematics = linspace(0, length(fTip1_X_smooth_HOC_OnOn)/60, length(fTip1_X_smooth_HOC_OnOn));

% Plot normalized kinematic data over time
figure;
plot(xTime_Kinematics, fTip1_X_normalized);
xlim([0, round(max(xTime_Kinematics)-1)])
ylim([0, round(max(fTip1_X_normalized)+1)]) 
xlabel('Time (s)')
ylabel('Normalized fTip1, X deflection (mm)')
title('Kinematics: Finger Tip X Position');


%% Create synchronized Kinematics and LFP plot for HOC task - On Med, On Stim - HOC
% 'trimFrame_int_OffOn' has the indices for the LFP data after being interpolated

% Find the closest matching indices in the interpolated kinematic data
[~, startInd_Kin] = min(abs(trimFrame_int_OnOn - HOC_start_index_OnOn));
[~, stopInd_Kin] = min(abs(trimFrame_int_OnOn - HOC_stop_index_OnOn));

% Trim the interpolated kinematic data to match the 'Hand OC' segment
dlc_trimmed_HOC_OnOn_int = dlc_trimmed_OnOn_int(startInd_Kin:stopInd_Kin, :);

% Find the corresponding indices in the LFP time vector
[~, startIND_LFP] = min(abs(ts_LFP_OnOn - ts_DLC_OnOn(HOC_start_index_OnOn)));
[~, endIND_LFP] = min(abs(ts_LFP_OnOn - ts_DLC_OnOn(HOC_stop_index_OnOn)));

% Extract the LFP data for Hand OC
LFP_HOC_OnOn = L_streamOfInt_OnOn(startIND_LFP:endIND_LFP);

% Calculate the corresponding time vector for the trimmed kinematic data
xTime_Kinematics_HOC = linspace(0, (stopInd_Kin - startInd_Kin + 1) / fs, stopInd_Kin - startInd_Kin + 1);

% Normalize the trimmed kinematic data so that the min is 0 and max is 50
fTip1_X_smooth_HOC = smoothdata(dlc_trimmed_HOC_OnOn_int.fTip1_x, 'gaussian', 70);
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
xTime_LFP = linspace(0, length(LFP_HOC_OnOn)/250, length(LFP_HOC_OnOn));
ecg = perceive_ecg(LFP_HOC_OnOn', 250, 0);
nexttile
plot(xTime_LFP, ecg.cleandata);
xlim([0, round(max(xTime_LFP)-1)])
% ylim([round(min(ecg.cleandata)-1), round(max(ecg.cleandata)+1)])
ylim([-50, 35])
ylabel('LFP Amplitude (µV)');
title('Time Synced STN LFP Data');


%% Compute and Plot Power Spectral Density (PSD), On Med, On Stim - HOC

% Filter ecg
ts_LFP_HOC_OnOn = 0:1/250:(height(LFP_HOC_OnOn)-1)/250;
ecgClean_HOC_OnOn = perceive_ecg(transpose(LFP_HOC_OnOn),250,0);
cleanLFP_HOC_OnOn = ecgClean_HOC_OnOn.cleandata;

pspectrum(cleanLFP_HOC_OnOn,250,'FrequencyLimits',[0 50],'FrequencyResolution',3) % 2 or 3
[bPxx_HOC_OnOn, bFxx_HOC_OnOn] = pspectrum(cleanLFP_HOC_OnOn,250,'FrequencyLimits',[0 50],'FrequencyResolution',3); %
bPxxP_HOC_OnOn = pow2db(bPxx_HOC_OnOn);

% Plot the power spectral density
figure;
%plot(freq, 10*log10(power));
plot(bFxx_HOC_OnOn, bPxxP_HOC_OnOn);
xlim([0 50])  % Limit x-axis to 50 Hz for better visualization
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectral Density');

%% bursting analysis, based on Torrecillos et al., 2018 (LFP stream 9, DLC session 7) - On Med, On Stim - HOC
% time-frq transform via complex Morlet wavelets (f0/σf = 7, f0 = 1-45 Hz, steps of 0.25 Hz).

% define 'lfp_data' as a LFP data matrix with dimensions [samples x trials]
lfp_data_HOC_OnOn = cleanLFP_HOC_OnOn;

fs = 250; % sampling rate
f0 = 1:0.25:45; % center frequency range from 1-45Hz in steps of 0.25Hz

% Frequency-Time Decomposition using complex Morlet wavelets
time = (1:length(lfp_data_HOC_OnOn)) / fs;
cwt_power = zeros(length(f0), length(lfp_data_HOC_OnOn));

for i = 1:length(f0)
    % create a sinusoidal wave modulated by a Gaussian window
    sinusoidal_wave = exp(2 * pi * 1i * f0(i) .* time); % complex cosine component of wavelet with an imaginary sine component
    Guassian_window = exp(-time.^2 / (2 * (f0(i) / 7)^2)); % tapers the sinusoid; (f0(i) / 7)^2 controls width of Gaussian
    wavelet = sinusoidal_wave .* Guassian_window; % applies wavelet across f0 range to create filters convolved with the LFP signal
    cwt_power(i, :) = abs(conv(lfp_data_HOC_OnOn, wavelet, 'same')).^2; % square to compute power
end

% time-domain bin size determined by 'VoicesPerOctave' (higher VPO = finer freq. res; lower time res?)
% continuous wavelet transform
figure;
cwt(lfp_data_HOC_OnOn,fs,'FrequencyLimits',[1 125],'VoicesPerOctave',14); % try 14, mults of 7
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
beta_power_time_courses_HOC_OnOn = mean(normalized_power(beta_range(peak_idx) + (-3:3), :), 1);

% Plot the beta power time course
figure('Renderer', 'opengl');
plot(time, beta_power_time_courses_HOC_OnOn);
xlim([0, round(max(time)-1)])
xlabel('Time (s)');
ylabel('Normalized Beta Power (µVp)');
title('LFP Beta Power')

%% Beta Burst Detection and Plotting for trimmed LFP

% Compute mean beta power
mean_beta_power_HOC_OnOn = mean(beta_power_time_courses_HOC_OnOn);
std_beta_power_HOC_OnOn = std(beta_power_time_courses_HOC_OnOn);


% Set threshold at the 75th percentile of the mean beta power
threshold = prctile(mean_beta_power_HOC_OnOn, 75);

% Define beta bursts as time points exceeding the threshold for more than 3 oscillatory cycles (3 oscillatory cycles = ~100ms)
min_duration = round(0.1 * fs); % Minimum duration of a burst in samples (100 ms)
% min_duration = round(0.066 * fs);  % Reduce minimum duration to 2 osc. cycles(66 ms)

above_threshold = beta_power_time_courses_HOC_OnOn > threshold;
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
plot(time, beta_power_time_courses_HOC_OnOn);
xlim([0, round(max(time)-1)])
xlabel('Time (s)');
ylabel('Normalized Beta Power');
title('LFP Beta Power')


% Plot the beta power time course and beta burst overlays using xregion (test with L_streamOfInt_s1)
% https://www.mathworks.com/help/matlab/ref/xregion.html

figure('Renderer', 'opengl');
p = plot(time, beta_power_time_courses_HOC_OnOn); % store handle to the beta power plot
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

