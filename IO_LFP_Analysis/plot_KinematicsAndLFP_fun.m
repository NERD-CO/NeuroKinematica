function plot_KinematicsAndLFP_fun(ts, dlcData, lfpData, streamIndex, sessionID)
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