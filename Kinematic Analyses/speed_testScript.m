% Define framerate of videos (time conversion factor)
fps = 100; % frames per second

% Convert distance units to mm (distance conversion factor)
pixels_to_mm = 2.109; % 232 mm / 110 pxl = 2.1091 mm per pixel

% compute euclidean distances between marker coord. across frames
euclidall = [];
for frame_i = 1:height(fTip1_x_dorsal)
    if frame_i ~= height(fTip1_x_dorsal)
        point1 = fTip1_x_dorsal(frame_i,:);
        point2 = fTip1_x_dorsal(frame_i + 1,:);
        euclidall(frame_i,1) = pdist2(point1 , point2);
    end
end

% Convert distance variables to mm usng conversion factor
euclidall = euclidall * pixels_to_mm; % converting euclidean distances to mm

frames = transpose([0:1:height(euclidall)]);
time = frames/fps; % seconds

time_step = 0.01; % step size
Df = diff(euclidall)/time_step; % velocity

speed = abs(Df); % speed

%[std_spd, mean_spd] = std(speed);





