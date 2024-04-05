











figure;
[cfs,frq] = cwt(cleanLFP,250,'FrequencyLimits',[2 100],'VoicesPerOctave',24);
t = 0:1/250:(length(cleanLFP)-1)/250;
absCFS = abs(cfs);
imagesc(t,frq,abs(cfs));
axis xy
% set(gca, 'YScale', 'log');
% imagesc(cfs)

betaIND = frq > 13 & frq < 30;
betaBAND = absCFS(betaIND,:)

betaBmean = mean(betaBAND,1)
betaBsm = smoothdata(betaBmean,'gaussian',30)
plot(betaBmean)
hold on
plot(betaBsm)



%%

% Assuming EEG_data is [channels x time]
signal1 = cleanLFP(1:end-1); % First signal
signal2 = pointX2sm; % Second signal

% Parameters for mscohere
window = 256; % Window size for the segments
noverlap = 128; % Number of points to overlap between segments
nfft = 512; % Number of FFT points
fs = 250; % Sampling frequency in Hz

% Calculate coherence
[Cxy, F] = mscohere(signal1, signal2, window, noverlap, nfft, fs);

% Plot coherence
figure;
plot(F, Cxy);
xlabel('Frequency (Hz)');
ylabel('Coherence');
title('Magnitude-Squared Coherence');

%%
[crossCorr, lags] = xcorr(signal1, signal2, 'coeff');

% Find lag with maximum correlation
[maxCorr, idx] = max(abs(crossCorr));
maxLag = lags(idx);

% Plot cross-correlation
figure;
plot(lags, crossCorr);
xlabel('Lag');
ylabel('Cross-Correlation Coefficient');
title('Cross-Correlation between Kinematic Data and Beta Instantaneous Frequency');
grid on;

% Highlight the lag with the maximum correlation
hold on;
plot(maxLag, maxCorr, 'ro');
legend('Cross-Correlation', 'Maximum Correlation Point');