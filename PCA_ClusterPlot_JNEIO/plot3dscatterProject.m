function [] = plot3dscatterProject(x,y,z,wallInfo,limInfo,color2use)

n_points = height(x);

% 2. Setup the Plot
hold on

% Set the axis limits to make room for the walls
xlim([limInfo(1,1), limInfo(1,2)]);
ylim([limInfo(2,1), limInfo(2,2)]);
zlim([limInfo(3,1), limInfo(3,2)]);

% --- PROJECTION 1: XY Plane (Bottom Wall) ---
% We use x and y data, and calculate the 2D convex hull
k_xy = convhull(x, y);

% Plot on the Z-wall (z = z_wall)
% We create a vector of the wall position matching the size of the hull
z_proj = repmat(wallInfo(3), size(k_xy));
patch(x(k_xy), y(k_xy), z_proj, 'FaceColor',color2use,...
    'FaceAlpha', 0.1, 'EdgeColor', color2use);

% --- D. Plot Projection on Z-Axis Floor (Black Crosses) ---
% We fix z to 'z_floor' and plot x vs y
plot3(x, y, repmat(wallInfo(3), n_points, 1),...
    'Color',color2use, 'LineWidth', 1);

% --- PROJECTION 2: XZ Plane (Back Wall) ---
% We use x and z data
k_xz = convhull(x, z);

% Plot on the Y-wall (y = y_wall)
y_proj = repmat(wallInfo(2), size(k_xz));
patch(x(k_xz), y_proj, z(k_xz), 'FaceColor',color2use,...
    'FaceAlpha', 0.1, 'EdgeColor', color2use);

% --- C. Plot Projection on Y-Axis Wall (Red Crosses) ---
% We fix y to 'y_wall' and plot x vs z
plot3(x, repmat(wallInfo(2), n_points, 1), z,...
    'Color',color2use, 'LineWidth', 1);

% --- PROJECTION 3: YZ Plane (Side Wall) ---
% We use y and z data
k_yz = convhull(y, z);

% Plot on the X-wall (x = x_wall)
x_proj = repmat(wallInfo(1), size(k_yz));
patch(x_proj, y(k_yz), z(k_yz), 'FaceColor',color2use,...
    'FaceAlpha', 0.1, 'EdgeColor', color2use);

% --- B. Plot Projection on X-Axis Wall (Green Crosses) ---
% We fix x to 'x_wall' and plot y vs z
plot3(repmat(wallInfo(1), n_points, 1), y, z,...
    'Color',color2use, 'LineWidth', 1);

title('3D Scatter with Convex Hull Projections');
view(3); % Set standard 3D view

end