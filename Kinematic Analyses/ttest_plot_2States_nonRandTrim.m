function ttest_plot_2States_nonRandTrim(state1_name, state1_data, state2_name, state2_data, outputDir)

% Check for NaN values and remove them
state1_data.amplitudes = state1_data.amplitudes(~isnan(state1_data.amplitudes));
state2_data.amplitudes = state2_data.amplitudes(~isnan(state2_data.amplitudes));

state1_data.widths = state1_data.widths(~isnan(state1_data.widths));
state2_data.widths = state2_data.widths(~isnan(state2_data.widths));

state1_data.peakDists = state1_data.peakDists(~isnan(state1_data.peakDists));
state2_data.peakDists = state2_data.peakDists(~isnan(state2_data.peakDists));

% Ensure both datasets have the same length
minLength_amplitudes = min(length(state1_data.amplitudes), length(state2_data.amplitudes));
state1_data.amplitudes = state1_data.amplitudes(1:minLength_amplitudes);
state2_data.amplitudes = state2_data.amplitudes(1:minLength_amplitudes);

minLength_widths = min(length(state1_data.widths), length(state2_data.widths));
state1_data.widths = state1_data.widths(1:minLength_widths);
state2_data.widths = state2_data.widths(1:minLength_widths);

minLength_peakDists = min(length(state1_data.peakDists), length(state2_data.peakDists));
state1_data.peakDists = state1_data.peakDists(1:minLength_peakDists);
state2_data.peakDists = state2_data.peakDists(1:minLength_peakDists);

% Perform t-tests
[~, p_value_amplitude] = ttest2(state1_data.amplitudes, state2_data.amplitudes);
[~, p_value_width] = ttest2(state1_data.widths, state2_data.widths);
[~, p_value_peakDists] = ttest2(state1_data.peakDists, state2_data.peakDists);

% Display p-values with condition names
disp(['p-value for peak amplitude between ' state1_name ' and ' state2_name ': ' num2str(p_value_amplitude)]);
disp(['p-value for peak width between ' state1_name ' and ' state2_name ': ' num2str(p_value_width)]);
disp(['p-value for peak distance between ' state1_name ' and ' state2_name ': ' num2str(p_value_peakDists)]);

% Create p-value table
p_value_table = table(p_value_amplitude, p_value_width, p_value_peakDists, 'VariableNames', {'p_val_Amplitude', 'p_val_Width', 'p_val_PeakDistance'});

% Save the table to a CSV file
p_values_filename = sprintf('%sp_values_%s_vs_%s.csv', outputDir, state1_name, state2_name);
writetable(p_value_table, p_values_filename);
% writetable(p_value_table, [outputDir filesep p_values_filename]);


% Plotting stat comparisons (between 2 states)
state1_means = [mean(state1_data.amplitudes), mean(state1_data.widths), mean(state1_data.peakDists)];
state2_means = [mean(state2_data.amplitudes), mean(state2_data.widths), mean(state2_data.peakDists)];

state1_sd = [std(state1_data.amplitudes), std(state1_data.widths), std(state1_data.peakDists)];
state2_sd = [std(state2_data.amplitudes), std(state2_data.widths), std(state2_data.peakDists)];

% Data organization for plotting
means_matrix = [state1_means; state2_means];
sd_matrix = [state1_sd; state2_sd];

% Labels
group_labels = {state1_name, state2_name};
metric_labels = {'Amplitudes', 'Intra-movement Durations', 'Inter-movement Durations'};

% Plotting
figure;
colors = [0.4 0.2 0.6; 0.2 0.7 0.8]; % RGB triplet per condition bar

for i = 1:3
    subplot(1, 3, i);
    b = bar(means_matrix(:, i), 'FaceColor', 'flat');
    hold on;

    % Assign colors per condition bar
    for cond = 1:2
        b.CData(cond,:) = colors(cond, :);
    end

    % Add text annotations
    for cond = 1:2
        x = b.XEndPoints(cond); % get the x-coordinate of each bar
        y = b.YEndPoints(cond); % get the y-coordinate of each bar
        text(x, y + 0.02, sprintf('%.2f (SD=%.2f)', means_matrix(cond,i), sd_matrix(cond,i)),...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end

    % Use standard deviation (sd_matrix) for the error bars
    errorbar(1:2, means_matrix(:, i), sd_matrix(:, i), '.k', 'LineWidth', 1.5);

    set(gca, 'XTickLabel', group_labels);
    title(['Comparison of ', metric_labels{i}]);
    ylabel(['Mean ', metric_labels{i}]);
    xlabel('Condition');

    hold off;
end
end