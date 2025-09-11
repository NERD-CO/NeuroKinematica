function [beta_power_time_courses, beta_peak_frequency, bursts, time] = burstingAnalysis_fun(lfp_data, fs)
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