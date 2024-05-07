function [kinematicsOUT] = getKINematics(pointX2sm)


ampMEAN = mean(pointX2sm,'omitnan');
ampSTD = std(pointX2sm,'omitnan');
ampThresh = (ampMEAN + (ampSTD/2))*0.4; % 0.75

[peaks, locs, widths, prominences] = findpeaks(pointX2sm,...
    MinPeakHeight=mean(pointX2sm), MinPeakDistance=15,...
    MinPeakProminence=ampThresh, Annotate ='extents');

figure;

findpeaks(pointX2sm,...
    MinPeakHeight=mean(pointX2sm), MinPeakDistance=15,...
    MinPeakProminence=ampThresh, Annotate ='extents')

%%
fps = 60;
pixels_to_mm = 2.109;
amplitudes = peaks * pixels_to_mm; % converting amplitudes to mm

% Compute timepoints from locs (vector of integer indices corresponding to video frame number)
timepoints = locs / 60; % Convert frame numbers to time (in seconds) using video sampling rate (Fs) conversion factor

% Convert frame-relative variables to seconds using time conversion factor
widths_fps = widths / fps; % converting widths to seconds
halfWidths = widths_fps / 2;
% slope

% Define unique name for the findpeaks results based on the current CSV name
% findpeaks_output = [outputDir , filesep , 'findpeaks_output_' ,...
    % moveTypeIDs{mmi}, '_', tmpCSV(1:end-44) '.csv']; % assumes tmpCSV is a string ending in '.csv'

% Create a table to store results based on computed variables
kinematicsOUT = table(timepoints, locs, peaks, amplitudes, prominences, widths_fps, halfWidths, 'VariableNames',...
    {'Timepoints', 'Locations', 'Peaks', 'Amplitudes', 'Prominences', 'Widths', 'HalfWidths'});


end