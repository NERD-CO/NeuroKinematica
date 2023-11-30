function ANOVA_plot_allStates(data, group_labels, measureName, yLabel, anovaTitle)

% Check if the lengths of data and group_labels are equal
if length(data) ~= length(group_labels)
    error('Data and group_labels must be the same length. Data length: %d, Group labels length: %d', length(data), length(group_labels));
end

% Clean the data to remove NaN values
validIndices = ~isnan(data);
data = data(validIndices);
group_labels = group_labels(validIndices);

% Make sure data and group labels are of the same length
assert(length(data) == length(group_labels), 'Data and group labels must be of the same length.');

% Perform ANOVA
[p_value_anova, tbl, stats] = anova1(data, group_labels, 'off');

% Display ANOVA results
disp([anovaTitle, ' ANOVA p-value: ', num2str(p_value_anova)]);
disp(tbl);

% Plot ANOVA Results
figure;
box_handle = boxplot(data, group_labels);
title([measureName, ' Comparison Across Conditions']);
ylabel(yLabel);
xlabel('Condition');

% Annotate with ANOVA p-value
x_limits = xlim;
y_limits = ylim;
text(x_limits(2) * 0.95, y_limits(2) * 0.95, ...
    sprintf('ANOVA p-value: %.3f', p_value_anova), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 10, 'FontWeight', 'bold');

% Calculate means and standard deviations for annotation
unique_groups = unique(group_labels);
for i = 1:length(unique_groups)
    groupData = data(strcmp(group_labels, unique_groups{i}));
    means(i) = mean(groupData);
    stds(i) = std(groupData);
end
% means = arrayfun(@(x) mean(data(strcmp(group_labels, x))), unique_groups);
% stds = arrayfun(@(x) std(data(strcmp(group_labels, x))), unique_groups);

% Get positions for annotations
boxes = findobj(gca, 'Tag', 'Box');
positions = arrayfun(@(x) x.XData(2), boxes);

% Reverse order since the boxes appear in reverse
positions = flip(positions);  % Annotate with means and standard deviations


for i = 1:length(means)
    text_position = [positions(i), y_limits(2) * 0.75];
    text(text_position(1), text_position(2), ...
        sprintf('Mean=%.2f\nSD=%.2f', means(i), stds(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontSize', 8, 'FontWeight', 'bold');
end

end
