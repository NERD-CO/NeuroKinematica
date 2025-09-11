function bursts = betaBurstDetection_fun(cfs, freqs, time, beta_range, fs)
% Compute mean and std of power in the beta range
beta_power = mean(abs(cfs(beta_range,:)).^2, 1);
mean_beta_power = mean(beta_power);
std_beta_power = std(beta_power);

% Detect bursts as periods where power is significantly above mean
burst_threshold = mean_beta_power + 2 * std_beta_power; % Example threshold
bursts = beta_power > burst_threshold;
end