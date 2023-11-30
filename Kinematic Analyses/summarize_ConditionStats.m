function summaryTable = summarize_ConditionStats(conditionName, amplitudes, widths, peakDists, outputDir)

% Compute overall mean, standard deviation, and variance
mean_amplitudes = mean(amplitudes);
std_amplitudes = std(amplitudes);
var_amplitudes = var(amplitudes);

mean_widths = mean(widths);
std_widths = std(widths);
var_widths = var(widths);

mean_peakDists = mean(peakDists);
std_peakDists = std(peakDists);
var_peakDists = var(peakDists);

% Create a table to store overall summary results
summaryTable = table(mean_amplitudes, std_amplitudes, var_amplitudes, ...
    mean_widths, std_widths, var_widths, ...
    mean_peakDists, std_peakDists, var_peakDists, ...
    'VariableNames', {'MeanAmplitude', 'StdAmplitude', 'VarAmplitude', ...
    'MeanWidth', 'StdWidth', 'VarWidth', ...
    'MeanPeakDist', 'StdPeakDist', 'VarPeakDist'});

% Assign row name to the summary table
summaryTable.Properties.RowNames = {sprintf('%s_Condition', conditionName)};

% Write summary stats table per condition to a CSV file
fileName = sprintf('%sfTipTracking_results_%s_summary_mm.csv', conditionName);
writetable(summaryTable, [outputDir filesep fileName],'WriteRowNames', true);

end