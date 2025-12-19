%% scatter plot w/ median line

data = randn(100, 4);

colors = rand(4, 3);

figure; 
for mi = 1:4
    x_data = ones(size(data, 1), 1) * mi; % 4 groups, move along x-axis
    y_data = data(:, mi); % all rows per col
    hold on
    s1 = scatter(x_data, y_data, 60, colors(mi, :), 'filled');

    % use properties to create jitter
    s1.XJitter = "rand"; % within rec/sq space
    s1.XJitterWidth = 0.25; % adjust as needed

    medianValue = median(y_data); % Calculate the median of the current group
    medianLn_width = 0.2;
    line([mi-medianLn_width mi+medianLn_width], [medianValue medianValue], 'Color', colors(mi, :), 'LineWidth', 2); % Plot median line
end

xticks(1:4);
