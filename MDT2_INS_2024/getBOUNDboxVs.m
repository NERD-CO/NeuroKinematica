function [xB_min,yB_min,B_width,B_height] = getBOUNDboxVs(xVALUES,yVALUES,fracEXP)

% Find the minimum and maximum coordinates
x_MIN = min(xVALUES);
x_MAX = max(xVALUES);
y_MIN = min(yVALUES);
y_MAX = max(yVALUES);

% Create the rectangle
width_box = x_MAX - x_MIN;
height_box = y_MAX - y_MIN;

% Expand the boundary by 10%
exp_factor = fracEXP;
B_width = width_box * (1 + exp_factor);
B_height = height_box * (1 + exp_factor);

% Adjust the position to keep the expansion centered
xB_min = x_MIN - (B_width - width_box) / 2;
yB_min = y_MIN - (B_height - height_box) / 2;

end