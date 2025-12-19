% JNE Color scheme

% https://colorbrewer2.org/#type=diverging&scheme=PRGn&n=7
% color brewer MATLAB fileexchange

% 7 classes
% diverging nature of data
% colorblind safe

% STN Depth encoded by 1st 4 rows (shades of purple)
% MoveType encoded by last 3 rows (shades of green or teal)

JNE_colors = (...
    [118,42,131;    % dark purple      (dorsal STN)
    175,141,195;    % lavender         (central STN)
    231,212,232;    % light purple     (ventral STN)
    128,128,128;    % change white --> grey (REST)
    217,240,211;    % light green      (Hand OC)
    127,191,123;    % green            (Hand PS)
    27,120,55]) ... % dark green    (Arm EF)
    ./ 256; % RGB --> 8-bit fraction for MATLAB

% Lavender shades (STN depth)
purpleShades = (...
    [118,42,131;      % dark purple    (dorsal STN)
    175,141,195;      % lavender       (central STN)
    231,212,232]) ... % light purple   (ventral STN)])
    ./ 256; % RGB --> 8-bit fraction for MATLAB

% Green shades (Movement Contexts)
greenShades = (...
    [128,128,128;   % grey             (REST)
    217,240,211;    % light green      (Hand OC)
    127,191,123;    % green            (Hand PS)
    27,120,55]) ... % dark green       (Arm EF)
    ./ 256; % RGB --> 8-bit fraction for MATLAB

% Teal shades:
tealShades = (...
    [199,234,229;   % light teal
    90,180,172;     % teal
    1,102,94]) ...  % dark teal
    ./ 256; % RGB --> 8-bit fraction for MATLAB

% Create blend of green + teal shades (more green) for Movement Context
GreenTeal_Shades = (...
    [(217+199)./2, (240+234)./2, (211+299)./2;   % light green  (Hand OC)
    (127+90)./2, (191+180)./2, (123+172)./2;     % green        (Hand PS)
    (27+1)./1, (120+102)./2, (55+94)./2]) ...    % dark green   (Arm EF)
    ./ 256; % RGB --> 8-bit fraction for MATLAB
]