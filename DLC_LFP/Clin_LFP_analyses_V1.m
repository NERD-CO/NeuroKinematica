% Clinical LFP analyses

%% Combine Percept LFP with DLC Video

% Data ingredients:
% 1) LFP - JSON Session Reports (1 report per hemisphere, multiple rows (stream recordings) per report [metadata informs ID of row])
% 2) Movement data (.mat files) and Movement Indices (.csvs)


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

        %     % % Optional: plot unfiltered data
        %     % figure; % Create new figure for each plot
        %     % plot(temp_BSTD_table.TimeDomainData{BSTD_i}); % blue
        %     % title(sprintf('File %d, Row %d', json_i, BSTD_i));
        %
        % filter out ECG for each row of time domain data in current JSON
        tempData_1 = transpose(BSTD_1_table.TimeDomainData{BSTD_i}); % Transpose raw data for current row
        ecg = perceive_ecg(tempData_1, fs, plotit);

        %     % % Optional: plot ecg-filtered data
        %     % hold on;
        %     % plot(ecg.cleandata); % orange
        %     % title(sprintf('ECG Filtered Data: File %d, Row %d', json_i, BSTD_i));
        %     % hold off;

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

% temp_JSON2_name = 'Report_Json_Session_Report_20230912T115939.json';
% temp_JSON2 = jsondecode(fileread(temp_JSON2_name));

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


%% filter out ecg

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


%% 9/12/2023 BrainSense LFP streams (12) in JSON Session Report 1 - Off Med

% js1_Row	sessDur(s)	DLC_sessID	DLC_sessDur(s) 	Hemisphere	Notes: case videos (10) for JSON Session Report 1 - Off Med
% 1	        32	        NA	            NA	        L	        % • L_baseline, (Off Med, Off Stim @ 0 mA)
% 2	        32	        NA		        NA          R	        % • R_baseline, (Off Med, Off Stim @ 0 mA)
% 3	        47	        session001	    56	        L	        % • session001, set1 (Off Med, Off Stim @ 0 mA)
% 4	        41	        session002		40          R	        % • session002, set1 (Off Med, Off Stim @ 0 mA)
% 5	        38	        session003	    48	        L	        % • session003, set2 (Off Med, Off Stim @ 0 mA)
% 6	        34	        session004		44          R	        % • session004, set2 (Off Med, Off Stim @ 0 mA)
% 7	        322	        session005	    336	        L	        % • session005 (Off Med, Stim Ramping)
% 8	        251	        session006		262         R	        % • session006 (Off Med, Stim Ramping)
% 9	        36	        session007	    45	        L	        % • session007, set1 (Off Med, On Stim @ max mA)
% 10	    29	        session008		37          R	        % • session008, set1 (Off Med, On Stim @ max mA)
% 11	    28	        session009	    35	        L	        % • session009, set2 (Off Med, On Stim @ max mA)
% 12	    26	        session010		36          R       	% • session010, set2 (Off Med, On Stim @ max mA)


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

L_streamsofInt_OffMed = {L_streamOfInt_bsln, L_streamOfInt_s1, L_streamOfInt_s3, L_streamOfInt_s5, L_streamOfInt_s7, L_streamOfInt_s9}; % w/ baseline
% L_streamsofInt_OffMed = {L_streamOfInt_s1, L_streamOfInt_s3, L_streamOfInt_s5, L_streamOfInt_s7, L_streamOfInt_s9}; % w/o baseline

% non-ecg-filtered
% maybe hp filter

% Titles for each plot in L_streamsofInt_OffMed
plotTitles = {
    'LSTN baseline (Off Med, Off Stim @ 0 mA)',
    'LSTN set1 (Off Med, Off Stim @ 0 mA)',
    'LSTN set2 (Off Med, Off Stim @ 0 mA)',
    'LSTN ramp (Off Med, Stim Ramping)',
    'LSTN set1 (Off Med, On Stim @ max mA)',
    'LSTN set2 (Off Med, On Stim @ max mA)'};

figure
for L_i = 1:length(L_streamsofInt_OffMed)
    ts_LFP = 0:1/250:(height(L_streamsofInt_OffMed{L_i})-1)/250;

    subplot(6,1,L_i)
    plot(ts_LFP, L_streamsofInt_OffMed{L_i})
    title(plotTitles{L_i})
    xlabel('Time (s)')
    ylabel('LFP Amplitude (uV)')
    grid on
end


R_streamsofInt_OffMed = {R_streamOfInt_bsln, R_streamOfInt_s2, R_streamOfInt_s4, R_streamOfInt_s6, R_streamOfInt_s8, R_streamOfInt_s10}; % w/ baseline
% R_streamsofInt_OffMed = {R_streamOfInt_s2, R_streamOfInt_s4, R_streamOfInt_s6, R_streamOfInt_s8, R_streamOfInt_s10}; % w/o baseline

% Titles for each plot in L_streamsofInt_OffMed
plotTitles = {
    'RSTN baseline (Off Med, Off Stim @ 0 mA)',
    'RSTN set1 (Off Med, Off Stim @ 0 mA)',
    'RSTN set2 (Off Med, Off Stim @ 0 mA)',
    'RSTN ramp (Off Med, Stim Ramping)',
    'RSTN set1 (Off Med, On Stim @ max mA)',
    'RSTN set2 (Off Med, On Stim @ max mA)'};

figure
for R_i = 1:length(R_streamsofInt_OffMed)
    ts_LFP = 0:1/250:(height(R_streamsofInt_OffMed{R_i})-1)/250;

    subplot(6,1,R_i)
    plot(ts_LFP, R_streamsofInt_OffMed{R_i})
    title(plotTitles{R_i})
    xlabel('Time (s)')
    ylabel('LFP Amplitude (uV)')
    grid on
end


%% Hilbert Transform (test w/ L_streamOfInt_s1)

fs = 250;
t = 1/250;

% Target signal (LFP) sampling rate at 250 Hz (samples per sec)
ts_LFP = 0:1/250:(height(L_streamOfInt_s1)-1)/250;

ecgClean_s1 = perceive_ecg(transpose(L_streamOfInt_s1),250,0);

cleanLFP_s1 = ecgClean_s1.cleandata;

% Calculate instantaneous phase and freq. using Hilbert transform
hilbert_s1 = hilbert(cleanLFP_s1);
inst_phase = angle(hilbert_s1);
inst_freq = diff(unwrap(inst_phase))/(2*pi*t);

% Filter instantaneous frequency for 13-30 Hz
bandpass_filt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',13,'HalfPowerFrequency2',30, ...
    'SampleRate',fs);
inst_freq_filtered = filtfilt(bandpass_filt, inst_freq);
inst_freq_filtered2 = [inst_freq_filtered , 0];

figure;
plot(ts_LFP,inst_freq_filtered2);
xlim([0 round(max(ts_LFP))])
ylabel('frequency (Hz)')
xlabel('time (seconds)')
title('Instantaneous Frequency of Filtered LFP Beta Band (13-30 Hz)');

% Compute and Plot Power Spectral Density (PSD)

% Use pspectrum to compute the power spectral density
[power, freq] = pspectrum(cleanLFP_s1, fs);

% Plot the power spectral density
figure;
plot(freq, 10*log10(power));
xlim([0 50])  % Limit x-axis to 50 Hz for better visualization
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density of LFP Signal');

%% bursting analysis (based on Torrecillos et al., 2018) test with L_streamOfInt_s1

% define 'lfp_data' as a LFP data matrix with dimensions [samples x trials]
lfp_data = L_streamOfInt_s1;

% Assuming 'movement_onset' is a vector indicating the time of movement onset for each trial

fs = 250;
f0 = 1:0.25:45; % Frequency range from 1 to 45 Hz in steps of 0.25 Hz

% Frequency-Time Decomposition using complex Morlet wavelets
time = (1:length(lfp_data)) / fs;
cwt_power = zeros(length(f0), length(lfp_data));
for i = 1:length(f0)
    wavelet = exp(2 * pi * 1i * f0(i) .* time) .* exp(-time.^2 / (2 * (f0(i) / 7)^2));
    cwt_power(i, :) = abs(conv(lfp_data, wavelet, 'same')).^2;
end

% time-domain bin size
% example to test
cwt(lfp_data,fs,'FrequencyLimits',[1 125],'VoicesPerOctave',7);

% sanity-check plot
imagesc(cwt_power)
set(gca, 'YDir','normal')
% remap y-axis to represent freq. steps
ytemp = yticks;
yticklabels(f0(ytemp))

% Normalize power for each frequency band
mean_power = mean(cwt_power, 2);
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
hold on;

% Overlay the detected bursts on the same plot
burst_indices = find(bursts);
for i = 1:length(burst_indices)
    burst_start = burst_indices(i);
    while i < length(burst_indices) && burst_indices(i + 1) - burst_indices(i) == 1
        i = i + 1;
    end
    burst_end = burst_indices(i);
    patch([time(burst_start), time(burst_end), time(burst_end), time(burst_start)], ...
          [min(beta_power_time_courses), min(beta_power_time_courses), ...
           max(beta_power_time_courses), max(beta_power_time_courses)], ...
          [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end

% https://www.mathworks.com/help/matlab/ref/xregion.html

hold off;
xlabel('Time (s)');
ylabel('Beta Power');
title('Beta Bursting Dynamics');
legend('Beta Power', 'Detected Bursts');


%% Loop through each LFP stream in Left STN
for i = 1:length(L_streamsofInt_OffMed)
    lfp_data = L_streamsofInt_OffMed{i};  

    % Hilbert Transform
    [inst_freq_filtered, ts_LFP] = HilbertTransform(lfp_data, fs);

    % Plotting
    figure;
    plot(ts_LFP, inst_freq_filtered);
    xlim([0 round(max(ts_LFP))]);
    ylabel('Frequency (Hz)');
    xlabel('Time (s)');
    title(sprintf('Instantaneous Frequency of filtered LFP Beta Band (13-30 Hz) for Stream %d', i));

    % Bursting Analysis
    [beta_power_time_courses, beta_peak_frequency, bursts, time] = burstingAnalysis(lfp_data, fs);

    % Plotting
    figure('Renderer', 'opengl');
    plot(time, beta_power_time_courses);
    hold on;
    burst_indices = find(bursts);
    for j = 1:length(burst_indices)
        burst_start = burst_indices(j);
        while j < length(burst_indices) && burst_indices(j + 1) - burst_indices(j) == 1
            j = j + 1;
        end
        burst_end = burst_indices(j);
        patch([time(burst_start), time(burst_end), time(burst_end), time(burst_start)], ...
              [min(beta_power_time_courses), min(beta_power_time_courses), ...
               max(beta_power_time_courses), max(beta_power_time_courses)], ...
              [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    end
    hold off;
    xlabel('Time (s)');
    ylabel('Beta Power');
    title(sprintf('Beta Bursting Dynamics for Stream %d', i));
    legend('Beta Power', 'Detected Bursts');
end

%% Loop through each LFP stream in Right STN
for i = 1:length(R_streamsofInt_OffMed)
    lfp_data = R_streamsofInt_OffMed{i};  

    % Hilbert Transform
    [inst_freq_filtered, ts_LFP] = HilbertTransform(lfp_data, fs);

    % Plotting
    figure;
    plot(ts_LFP, inst_freq_filtered);
    xlim([0 round(max(ts_LFP))]);
    ylabel('Frequency (Hz)');
    xlabel('Time (s)');
    title(sprintf('Instantaneous Frequency of filtered LFP Beta Band (13-30 Hz) for Stream %d', i));

    % Bursting Analysis
    [beta_power_time_courses, beta_peak_frequency, bursts, time] = burstingAnalysis(lfp_data, fs);

    % Plotting
    figure('Renderer', 'opengl');
    plot(time, beta_power_time_courses);
    hold on;
    burst_indices = find(bursts);
    for j = 1:length(burst_indices)
        burst_start = burst_indices(j);
        while j < length(burst_indices) && burst_indices(j + 1) - burst_indices(j) == 1
            j = j + 1;
        end
        burst_end = burst_indices(j);
        patch([time(burst_start), time(burst_end), time(burst_end), time(burst_start)], ...
              [min(beta_power_time_courses), min(beta_power_time_courses), ...
               max(beta_power_time_courses), max(beta_power_time_courses)], ...
              [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    end
    hold off;
    xlabel('Time (s)');
    ylabel('Beta Power');
    title(sprintf('Beta Bursting Dynamics for Stream %d', i));
    legend('Beta Power', 'Detected Bursts');
end

%% Separate LFP streams and corresponding DLC sessions of interest by condition

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

% EUC indicies
cd(mainDir2)
eucINDICIES = readtable("EUC_Indicies_1.xlsx");

% Find all rows in eucINDICIES where the 'moveID' is 'Tablet'
tablet_indices = find(strcmp(eucINDICIES.moveID, 'Tablet'));

% Define row indices in eucINDICIES for each 'moveID'
HandMov_indices = find(strcmp(eucINDICIES.moveID, 'HAND OC'));
PronSup_indices = find(strcmp(eucINDICIES.moveID, 'HAND PS'));
FlexExtend_indices = find(strcmp(eucINDICIES.moveID, 'ARM EF'));


%% test - Isolate LFP and Video of iterest (LFP stream 3, DLC session 1)

% L_streamOfInt_s1 (LFP stream 3) corresponds with DLC session001 (movement video 1), (Off Med, Off Stim @ 0 mA)

% Define LFP stream of interest: 
streamOfInt_1 = L_streamOfInt_s1; % 11688x1 double

% Extract DLC video of interest
session_id_1 = eucINDICIES.videoID{1};

% Extract dlc Tablet start frame index 
tablet_start_index_1 = eucINDICIES.StartInd(tablet_indices(1));
tablet_stop_index_1 = eucINDICIES.StopInd(tablet_indices(1));

% Load mat file corresponding to DLC video session
cd(DLC_mat_dir)
load('dlcDAT_20230912_session001_idea08_resnet50_rightCam-0000DLC.mat','outDATA')
dlc_labDAT_1 = outDATA;

% Trim mat file corresponding to DLC video session
dlc_trimmed_1 = dlc_labDAT_1(tablet_start_index_1:tablet_stop_index_1,:); % 2821x40 table

% Synchronize corresponding LFP stream (trim by dlc Tablet start frame index) 

% calculate offset from Video_start to Tablet_start 
% trim offset
% synchronize Lfp_start with Tablet_start


%% Sampling Rate Conversions and Calculations

totalNumSamples_Vid = height(dlc_trimmed_1);

totalNumSecs = totalNumSamples_Vid/60; % 60 fps

totalNumSamples_LFP = floor(totalNumSecs*250); % 250 samples per second
% y = floor(a) rounds fi object a to the nearest integer in the direction of negative infinity and returns the result in fi object y

% Original signal (video) sampling rate at 60 Hz (fps)
ts_DLC = 0:1/60:(height(dlc_trimmed_1)-1)/60;

% Target signal (LFP) sampling rate at 250 Hz (samples per sec)
ts_LFP = 0:1/250:(height(streamOfInt_1)-1)/250;

% Interpolate - upsample kinematic data to match LFP sampling rate
allColNs = dlc_trimmed_1.Properties.VariableNames;
dlc_trimmed_1_int = table;
for coli = 1:width(dlc_trimmed_1)

    % Tmp col
    tmpCol = allColNs{coli};

    % % Interpolation
    % x_250 = interp1(ts_DLC, transpose(dlc_lab2use.(tmpCol)), ts_LFP, 'spline');
    x_250 = interp1(ts_DLC, dlc_trimmed_1.(tmpCol), ts_LFP);

    dlc_trimmed_1_int.(tmpCol) = transpose(x_250);
end

%trimFrB_int = linspace(FRAMEstart,FRAMEend, length(ts_LFP));

%% Kinematics and LFP plot

tiledlayout(2,1,"TileSpacing","tight")
xTime = ts_LFP;

fTip1_X_smooth = smoothdata(dlc_trimmed_1_int.fTip1_x,'gaussian',70);
nexttile
plot(xTime,fTip1_X_smooth);
xlim([0 round(max(ts_LFP))])
xlabel('Time (s)')
ylabel('fTip1, X deflection')

ecg = perceive_ecg(transpose(streamOfInt_1),250,0);

nexttile
plot(xTime,ecg.cleandata);
xlim([0 round(max(ts_LFP))])
xlabel('Time (s)')
ylabel('LFP (uV)')

%%

% Extract start and stop frame indices for 'Hand OC'
HOC_start_index_1 = eucINDICIES.StartInd(HandMov_indices(1));
HOC_stop_index_1 = eucINDICIES.StopInd(HandMov_indices(1));


%% loop for test above




%%
% OFF Med
DLC_s1 = videoIDs(1:3); % set 1 [OFF Med, OFF Stim]
DLC_s3 = videoIDs(4:6); % set 2 [OFF Med, OFF Stim]
% DLC_5_ramp % OFF Med
DLC_s7 = videoIDs(7:9); % set 1 [OFF Med, ON Stim]
DLC_s9 = videoIDs(10:12); % set 2 [OFF Med, ON Stim]

% ON Med
DLC_13 = videoIDs(13:15); % set 1 [ON Med, OFF Stim]
DLC_15 = videoIDs(16:18); % set 2 [ON Med, OFF Stim]
% DLC_18_ramp % ON Med
DLC_20 = videoIDs(19:21); % set 1 [ON Med, ON Stim]
DLC_22 = videoIDs(22:24); % set 2 [ON Med, ON Stim]

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

%% subfunctions


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
    % Bursting Analysis: Time-Frequency Decomposition using complex Morlet wavelet transform
    f0 = 1:0.25:45;  % Frequency range from 1 to 45 Hz in steps of 0.25 Hz
    time = (1:size(lfp_data, 1)) / fs;
    cwt_power = zeros(length(f0), size(lfp_data, 1));
    for i = 1:length(f0)
        % wavelet = complex exponential (sine wave) modulated by a Gaussian envelope. 
        % (f0(i) / 7) controls the width of the Gaussian envelope, with a constant ratio of f0/σf = 7, 
        % σf = standard deviation of the Gaussian in the frequency domain
        wavelet = exp(2 * pi * 1i * f0(i) .* time) .* exp(-time.^2 / (2 * (f0(i) / 7)^2)); 
        cwt_power(i, :) = abs(conv(lfp_data, wavelet, 'same')).^2; % convolves the LFP data with the complex Morlet wavelet and then takes the absolute value squared to obtain the power
    end

    % Normalize power for each frequency band
    mean_power = mean(cwt_power, 2); % 2nd dim is the time dimension
    std_power = std(cwt_power, 0, 2);
    normalized_power = (cwt_power - mean_power) ./ std_power;

    % Beta Peak Selection
    beta_range = find(f0 >= 13 & f0 <= 30);
    [~, peak_idx] = max(mean(normalized_power(beta_range, :), 2));
    beta_peak_frequency = f0(beta_range(peak_idx));

    % Burst Detection
    % calculate average beta power time courses across a 6-Hz-wide frequency band centered on the peak beta frequency
    beta_power_time_courses = mean(normalized_power(beta_range(peak_idx) + (-3:3), :), 1);
    
    % compute mean beta power
    mean_beta_power = mean(beta_power_time_courses);
    
    % set threshold at the 75th percentile of the mean beta power across the session
    threshold = prctile(mean_beta_power, 75);

    % Define beta bursts as time points exceeding the threshold for more than 3 oscillatory cycles (3 oscillatory cycles = ~100ms)
    min_duration = round(0.1 * fs);  
    above_threshold = beta_power_time_courses > threshold;
    above_threshold = above_threshold(:);  % Ensure it's a column vector
    bursts = zeros(size(above_threshold));  % Initialize bursts array
    [start_indices, end_indices] = findConsecutiveIndices(above_threshold);  % Find start and end indices of consecutive regions above threshold
    for k = 1:length(start_indices)
        if (end_indices(k) - start_indices(k) + 1) >= min_duration
            bursts(start_indices(k):end_indices(k)) = 1;  % Mark burst regions
        end
    end
end

function [start_indices, end_indices] = findConsecutiveIndices(logical_array)
    % Find start and end indices of consecutive true regions in a logical array
    d = diff([0; logical_array(:); 0]);
    start_indices = find(d > 0);
    end_indices = find(d < 0) - 1;
end

