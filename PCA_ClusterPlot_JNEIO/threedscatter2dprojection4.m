

%%
n_points = 100;
% 1. Generate Random 3D Data
rng(42); % For reproducibility
x = randn(100,1);
y = randn(100,1);
z = randn(100,1);
%%


x = score(:,1);
y = score(:,2);
x = score(:,3);

%%

% 2. Setup the Plot
figure;
% scatter3(x, y, z, 15, 'filled', 'k'); % Plot original points in black
hold on;
% grid on;
% axis equal;
% xlabel('X'); ylabel('Y'); zlabel('Z');

% 3. Define "Wall" Locations
% We push the walls slightly beyond the data limits so they don't overlap
x_wall = min(x) - 1;
y_wall = max(y) + 1;
z_wall = min(z) - 1;

% Set the axis limits to make room for the walls
xlim([x_wall, max(x)+1]);
ylim([min(y)-1, y_wall]);
zlim([z_wall, max(z)+1]);
%%

mtOrder = {'HAND OC','HAND PS','ARM EF'};
% threeClusts = randsample(3,100,1);
threeClusts = zeros(length(mtPerUnit),1);
for mii = 1:3
    mv = mtOrder{mii};
    idx = mtPerUnit == mv;
    threeClusts(idx) = mii; 
end

%% Cluster eval

[silvals] = silhouette(score(:,3),threeClusts);

clusEv = zeros(3,1);
for svi = 1:3
    mv = mtOrder{svi};
    idx = mtPerUnit == mv;
    medSV = median(silvals(idx));
    clusEv(svi) = medSV; 
end

moveClusts = strings(size(threeClusts));
for mii = 1:3
    moveClusts(threeClusts == mii) = mtOrder{mii};
end

dataTable = table(silvals,moveClusts,'VariableNames',{'SilVals','ClustIDs'});

factors = moveClusts; % Define factors for ANOVA

aov = anova(factors,silvals,FactorNames="ClustIDs")

m = multcompare(aov, "ClustIDs", "CriticalValueType", "bonferroni")

close all
figure;
tiledlayout(3,1)
for ppii = 1:3

    nexttile
    histogram(silvals(threeClusts == ppii),15)

    xline(clusEv(ppii),'-k','median')

    xlim([-1 1])
    title(mtOrder{ppii})
end



%%

close all
colorSS = ([ 
    217-15,240-15,211-15;     % HAND OC
    127,191,123;     % HAND PS
    27,120,55] ...   % ARM EF
    ./ 255);
hold on

for i = 1:3
    % --- PROJECTION 1: XY Plane (Bottom Wall) ---
    % We use x and y data, and calculate the 2D convex hull

    cInd = threeClusts == i;

    n_pointsi = sum(cInd);

    xi = x(cInd);
    yi = y(cInd);
    zi = z(cInd);

    k_xy = convhull(xi, yi);

    % Plot on the Z-wall (z = z_wall)
    % We create a vector of the wall position matching the size of the hull
    z_proj = repmat(z_wall, size(k_xy));
    p1 = patch(xi(k_xy), yi(k_xy), z_proj, colorSS(i,:), 'FaceAlpha',...
        0.1, 'EdgeColor', colorSS(i,:));

    p1.FaceColor = 'none';

    % --- D. Plot Projection on Z-Axis Floor (Black Crosses) ---
    % We fix z to 'z_floor' and plot x vs y
    % plot3(xi, yi, repmat(z_wall, n_pointsi, 1), [colorSS(i),'o'],...
    %     'LineWidth', 1);

    % --- PROJECTION 2: XZ Plane (Back Wall) ---
    % We use x and z data
    k_xz = convhull(xi, zi);

    % Plot on the Y-wall (y = y_wall)
    y_proj = repmat(y_wall, size(k_xz));
    p2 = patch(xi(k_xz), y_proj, zi(k_xz), colorSS(i,:), 'FaceAlpha',...
        0.1, 'EdgeColor', colorSS(i,:));

    p2.FaceColor = 'none';

    % --- C. Plot Projection on Y-Axis Wall (Red Crosses) ---
    % We fix y to 'y_wall' and plot x vs z
    % plot3(xi, repmat(y_wall, n_pointsi, 1), zi, [colorSS(i),'o'],...
    %     'LineWidth', 1);

    % --- PROJECTION 3: YZ Plane (Side Wall) ---
    % We use y and z data
    k_yz = convhull(yi, zi);

    % Plot on the X-wall (x = x_wall)
    x_proj = repmat(x_wall, size(k_yz));
    p3 = patch(x_proj, yi(k_yz), zi(k_yz), colorSS(i,:), 'FaceAlpha', 0.1,...
        'EdgeColor', colorSS(i,:));

    p3.FaceColor = 'none';

    % --- Plot Projection on X-Axis Wall 
    % fix x to 'x_wall' and plot y vs z
    % plot3(repmat(x_wall, n_pointsi, 1), yi, zi, [colorSS(i),'o'],...
    %     'LineWidth', 1);

    num_points = sum(cInd);
    original_indices = 1:n_pointsi;

    new_indices = linspace(1, n_pointsi, n_pointsi * 10);

    x_smooth = interp1(original_indices, xi, new_indices, 'spline');
    y_smooth = interp1(original_indices, yi, new_indices, 'spline');
    z_smooth = interp1(original_indices, zi, new_indices, 'spline');

    plot3(x_smooth, y_smooth, z_smooth, 'Color',colorSS(i,:),'LineWidth', 3);

end

title('');
view(3); 