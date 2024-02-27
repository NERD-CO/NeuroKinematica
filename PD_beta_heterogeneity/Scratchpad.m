%% PD beta peak feature heterogeneity

% Define data directories

mainDir = 'C:\Users\radclier\OneDrive\Documents\CU Denver - Anschutz Medical BioE PhD\Neuro_ThompsonLab\PD beta peak feature heterogeneity';

mainDir2 = 'C:\Users\radclier\OneDrive\Documents\CU Denver - Anschutz Medical BioE PhD\Neuro_ThompsonLab\PD beta peak feature heterogeneity\Subject Data';


cd(mainDir2)

mat_dir1 = dir('*.mat');
mat_dir2 = {mat_dir1.name}; % list of mat files

addpath('C:\MATLAB\GitHub\NeuroKinematica\PD_beta_heterogeneity') 
for i = 1:length(mat_dir2)
    outStruct = fixVecLength(mat_dir2{i})
    % overwrite orig. file with new mat file
    save(mat_dir2{i}, 'outStruct')
end

%% Sensing Session Data (15 streams per hemisphere, 30 streams per pt)

% outershell features
% 1. Pt. ID
% 2. Hem
% 3. PD symptom scores – get from Dr. Kern notes (choose clinPref)
% bsln hemibody symptoms
% optimal stim improvement


%% create Structured output (.mat file) of Sensing Session Data

% i.    Structured output (.mat file)

% addpath('C:\MATLAB\GitHub\MDT_Phase2\CreateSessionStructure')
% createFolderStructure_mdt2();


%% Assess sensing data isolated by patient ID and hemisphere via nested for loop
% Structured output (.mat file) contains 15 total contact configs per pt. hem

cd(mainDir2)

% initialize a summary table 'VariableNames'
pt_ID_all = cell(2000,1);
hem_all = cell(2000,1);
conditionName_all = cell(2000,1);
contactName_all = cell(2000,1);
fxx_all = cell(2000,1);
uVp_t_all = cell(2000,1);
PxxP_all = cell(2000,1);

% counter
row_counter = 1;

% outer loop: .mat_filenames
mat_dir1 = dir('*.mat');
mat_dir2 = {mat_dir1.name}; % list of mat files

% Initialize session_container as a cell array if the lengths can vary
% session_container = cell(:, 1);

for f_i = 1:length(mat_dir2)
    filename = mat_dir2{f_i}
    name_parts = split(filename, '_');
    pt_ID = name_parts{1}   % 1. Define participant ID
    hem = name_parts{2}     % 2. Define hemisphere

    load(filename)

    % ii. Tab-based structure (.mat file)
    struct_mat = outStruct; % 1x1 struct with 3-4 fields (that are also 1x1 structs)

    % iterate through each LFP recording condition
    conditionNames = fieldnames(struct_mat); % 4x1 cell
    for condition_i = 1:numel(conditionNames)                               % Layer 1: condition
        conditionName = conditionNames{condition_i};
        sensingSessions = struct_mat.(conditionName); % 1x1 struct with 3 fields

        % iterate through each field (session) within the condition
        sessionNames = fieldnames(sensingSessions); % 3x1 cell

        % Skip 'INFO' session names
        sessionNames = sessionNames(~matches(sessionNames, 'INFO'));

        % extract each contact configuration from the sensing configs struct
        contactConfigs = fieldnames(sensingSessions.S1);

        % iterate through each contact configuration within the condition
        for contact_i = 1:numel(contactConfigs)                             % Layer 2: contact per condition
            contactName = contactConfigs{contact_i};

            temp_RowNum = numel(sensingSessions.S1.(contactName));

            session_container = zeros(temp_RowNum,3);

            % iterate through each contact config per session per condition
            for session_i = 1:numel(sessionNames)                           % Layer 3: contact stream per session per condition
                sessionName = sessionNames{session_i};
                contact_stream = sensingSessions.(sessionName).(contactName);


                session_container(:, session_i) = contact_stream;
                % session_container{session_i} = contact_stream;

            end

            % compute average across sessions per contact config within the condition
            session_mean = mean(session_container, 2);

            % compute PSD on the average across sessions per contact config within the condition
            % [fxx, pxx, uVp_t, PxxP] = computePSD(stream_mean);
            window = hanning(250);
            [pxx, fxx] = pwelch(session_mean, window, 125, 256, 250, 'onesided');  % closest approx to what MDT tablet displays (125 overlap, 256 bits, 250 window)
            uVp_t = sqrt(pxx).*rms(window).*sqrt(2).*2.*250/256;  % method for scaling PSD --> microvolt peak output
            % PxxP(:,di2) = 10*log10(pxx);
            PxxP = pow2db(pxx); % method for scaling PSD --> decibil conversion output

            conditionName_all{row_counter} = conditionName;
            contactName_all{row_counter} = contactName;
            fxx_all{row_counter} = fxx;
            uVp_t_all{row_counter} = uVp_t;
            PxxP_all{row_counter} = PxxP;
            pt_ID_all{row_counter} = pt_ID;
            hem_all{row_counter} = hem;

            row_counter = row_counter + 1;

        end
    end

end


% find row where first empty cell in located
trim_vec = cellfun(@(x) ~isempty(x), conditionName_all, UniformOutput=true);

pt_ID_all2 = pt_ID_all(trim_vec);
hem_all2 = hem_all(trim_vec);
condition_all2 = conditionName_all(trim_vec);
contactName_all2 = contactName_all(trim_vec);
fxx_all2 = fxx_all(trim_vec);
uVp_t_all2 = uVp_t_all(trim_vec);
PxxP_all2 = PxxP_all(trim_vec);

% create a summary table with the following columns: pt_ID, hem, condition, contact, fxx, uVp_t, PxxP
summaryTable = table(pt_ID_all2, hem_all2, condition_all2, contactName_all2, fxx_all2, uVp_t_all2, PxxP_all2, 'VariableNames', {'Participant ID', 'Hemisphere', 'Condition', 'Contact(s)', 'fxx', 'uVp', 'PxxP'});

% save summaryTable as both a .mat and .csv
cd(mainDir)
save('summaryTable.mat', 'summaryTable');


%%

% do I need to save out session_mean per pt. hem?


%%

% load summaryTable
load("summaryTable.mat");

% initialize variables for updated summaryTable
ptID = zeros(height(summaryTable),1);
hemID = cell(height(summaryTable),1);
conditions_All = cell(height(summaryTable),1);
contacts_All = cell(height(summaryTable),1);
meanPOW = zeros(height(summaryTable),1);
peakValue = zeros(height(summaryTable),1);
maxPOW = zeros(height(summaryTable),1);
standDevPOW = zeros(height(summaryTable),1);
aucPOW = zeros(height(summaryTable),1);
skewPOW = zeros(height(summaryTable),1);
kurtPOW = zeros(height(summaryTable),1);
medianFreq = zeros(height(summaryTable),1);

% iterate through each row in the summaryTable
for pt_i = 1:height(summaryTable)

    % Extract necessary information for each measurement
    ptID(pt_i) = str2double(summaryTable.('Participant ID'){pt_i}); 
    hemID{pt_i} = summaryTable.('Hemisphere'){pt_i};
    conditions_All{pt_i} = summaryTable.('Condition'){pt_i};
    contacts_All{pt_i} = summaryTable.('Contact(s)'){pt_i};
    pxx = summaryTable.PxxP{pt_i}; % PSD values
    fxx = summaryTable.fxx{pt_i}; % Frequency values
    
    % for loop through freq. bands 
    % theta
    % beta

    % Compute statistics from pxx
    meanPOW(pt_i) = mean(pxx);
    [maxPOW(pt_i), peakIdx] = max(pxx);
    peakValue(pt_i) = fxx(peakIdx);
    standDevPOW(pt_i) = std(pxx);
    aucPOW(pt_i) = trapz(fxx, pxx);
    skewPOW(pt_i) = skewness(pxx);
    kurtPOW(pt_i) = kurtosis(pxx);
    
    % Compute median frequency
    totalPower = sum(pxx);
    cumulativePower = cumsum(pxx);
    halfTotalPower = totalPower / 2;
    medianIndex = find(cumulativePower >= halfTotalPower, 1, 'first');
    medianFreq(pt_i) = fxx(medianIndex);
end

% Append the calculated variables as new columns to 'summaryTable'
summaryTable.meanPOW = meanPOW;
summaryTable.peakValue = peakValue;
summaryTable.maxPOW = maxPOW;
summaryTable.standDevPOW = standDevPOW;
summaryTable.aucPOW = aucPOW;
summaryTable.skewPOW = skewPOW;
summaryTable.kurtPOW = kurtPOW;
summaryTable.medianFreq = medianFreq;

% save summaryTable as both a .mat and .csv
save('updatedSummaryTable.mat', 'summaryTable');


%%

%%

% Load summaryTable
load("summaryTable.mat");
dimdataTab = summaryTable;

% 5th column onwards in summaryTable contains the data for t-SNE (fxx, uVp, PxxP)
cellArray = table2array(dimdataTab(:,5:end));

% Initialize an empty cell array to store the numeric data
numericData = cell(size(cellArray, 1), size(cellArray, 2));

% Loop through each cell in the cell array and extract the numeric data
for i = 1:size(cellArray, 1)
    for j = 1:size(cellArray, 2)
        % Extract the numeric data from the cell
        numericData{i, j} = cellArray{i, j};
    end
end

% Convert the cell array of numeric data to a numeric matrix
dimDataNum = cell2mat(numericData);

% Set Perplexity to a value less than the number of rows in dimDataArray
perplexityValue = min(size(dimDataNum, 1) - 1, 15);

% t-SNE with adjusted Perplexity value
tsneOutdim = tsne(dimDataNum, 'Algorithm', 'exact', 'Distance', 'mahalanobis', 'Perplexity', perplexityValue);

% Define new color scheme
newCOLORS = [35, 61, 77;
             254, 127, 45;
             252, 202, 70;
             161, 193, 129;
             97, 155, 138];
newCOLORS2 = newCOLORS / 255;

% Overlay outlines for each subject
for subii = 1:5
    tsneSUB = tsneOutdim(ismember(dimdataTab.Subject,subii),:);
    scatter(tsneSUB(:,1),tsneSUB(:,2),20,newCOLORS2(subii,:),"filled")
    hold on
    tmpBound = boundary(tsneSUB,0.05);
    
    xVals = tsneSUB(tmpBound,1);
    yVals = tsneSUB(tmpBound,2);
    
    centroid_x = mean(xVals);
    centroid_y = mean(yVals);
    
    scale_factor = 1.2; % Adjust as necessary
    
    % Move each point away from the centroid
    new_x = centroid_x + scale_factor * (xVals - centroid_x);
    new_y = centroid_y + scale_factor * (yVals - centroid_y);
    
    plot(new_x,new_y,'Color',newCOLORS2(subii,:))
    
    ftmp = fill(new_x,new_y,newCOLORS2(subii,:),'FaceAlpha',0.3);
    ftmp.EdgeColor = newCOLORS2(subii,:);
    ftmp.LineWidth = 1.2;
end

legend({'','','1','','','2','','','3','','','4','','','5'},"Location","southeast")

xticks([-25 0 25])
xlabel('tSNE 1')
yticks([-20 0 20])
ylabel('tSNE 2')

axis square
ax1 = gca;
ax1.TitleHorizontalAlignment = 'left';
title('Subject separation - beta power parameters')

% Compute silhouette values for each subject
figure;
[silhouVals_Subject,~] = silhouette(tsneOutdim,dimdataTab.Subject);

% Plot histograms of silhouette values for each subject
figure;
tiledlayout(5,1,'TileSpacing','tight','Padding','compact')
for subii = 1:5
    nexttile
    h = histogram(silhouVals_Subject(dimdataTab.Subject == subii));
    h.Normalization = 'probability';
    h.BinWidth = 0.1;
    h.FaceColor = newCOLORS2(subii,:);
    h.EdgeColor = "none";
    h.FaceAlpha = 0.5;
    
    meanLOC = mean(silhouVals_Subject(dimdataTab.Subject == subii));
    meanTxt = ['Ave. ' , num2str(round(meanLOC,2))];
    meanColr = newCOLORS2(subii,:);
    xline(meanLOC,'-',meanTxt,'Color',meanColr,'LabelHorizontalAlignment','left')
    
    xlim([-1 1])
    ylim([0 0.6])
    xlabel('Silhouette value')
    ylabel('Fraction of values')
end

set(gcf, 'Position', [1160 266 360 911])

































































%%

%% load (50 second recordings)
% Run PSD
% Use MDT variation of PSD computation


%% 3. PD symptom scores – get from Dr. Kern notes (choose clinPref)

% bsln hemibody symptoms
% optimal stim improvement


%% Filter LFP Timeseries Data

% fs = 250;
% nfft = 250;
% window = 250;
% overlap = 150;
% lfpsamplerate = 2;
%
%
% % field of int = average of sessions per contact
%
% streamOfInt = ex_stream1;
%
% % Define perceive_ecg function params and add function paths
% fs = 250;
% plotit = 0; % 0 = don't plot, 1 = plot
% addpath 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\perceive-master'
%
% % ecg = perceive_ecg(transpose(streamOfInt),250,0);
% % plot(ecg.cleandata);
% % yyaxis right
%
% [px,fx,tx] = pspectrum(streamOfInt, fs, 'spectrogram');
% % [pxx,fxx,txx] = pspectrum(ecg.cleandata, fs, 'spectrogram');
% % outputs: power matrix, frequency vec, time vec



%% compute LFP power / instantaneous LFP beta power and plot PSDs per session (using pspectrum function)

% % bin data - bin size
% % narrow, 2Hz windows
%
%
% % Hilbert transform - obtain the amplitude envelope of LFP
% % study the instantaneous power of the LFP data
% hilbertTransformed_LFP = hilbert(streamOfInt, fs);
% % hilbertTransformed_LFP = hilbert(ecg.cleandata, fs);
% amplitudeEnvelope = abs(hilbertTransformed_LFP);
%
% % bandpass filter for beta
% betaBand = [13 35]; % Beta frequency range
% bpFiltered_LFP = bandpass(streamOfInt, betaBand, fs);
% % bpFiltered_LFP = bandpass(ecg.cleandata, betaBand, fs);
%
% % plot PSD to vizualize strength of the variations (energy) as a function of frequency
% [pxx, fxx, txx] = pspectrum(bpFiltered_LFP, fs, 'spectrogram');
% plot(fxx, pxx);
% xlabel('Frequency (Hz)');
% ylabel('Power');
% title('Power Spectral Density');
%
%
% % Hilbert transform - obtain the amplitude envelope of LFP in the beta band.
% % study the instantaneous power of the LFP data
%
% hilbertTransformed_beta = hilbert(bpFiltered_LFP);
% betaEnvelope = abs(hilbertTransformed_beta);



%% functions

% % compute PSD on the average of all streams per contact per sensing config (condition)
% function [fxx, pxx, uVp_t, PxxP] = computePSD(data)
% [pxx, fxx] = pwelch(data, hanning(250), 125, 256, 250, 'onesided');         % closest approx to what MDT tablet displays (125 overlap, 256 bits, 250 window)
% uVp_t = sqrt(pxx) .* rms(hanning(250)) .* sqrt(2) .* 2 .* 250 / 256;        % method for scaling PSD --> microvolt peak output
% % PxxP(:,di2) = 10*log10(pxx);
% PxxP = pow2db(pxx);                                                         % method for scaling PSD --> decibil conversion output
% end
