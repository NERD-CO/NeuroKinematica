%% PD beta peak feature heterogeneity

% Define data directory

mainDir = 'C:\Users\radclier\OneDrive\Documents\CU Denver - Anschutz Medical BioE PhD\Neuro_ThompsonLab\PD beta peak feature heterogeneity';

cd(mainDir)

%% Analyze data isolated by Pt.ID and hemisphere

% outershell features
% 1. Pt. ID
% 2. Hem
% 3. PD symptom scores – get from Dr. Kern notes (choose clinPref)
% bsln hemibody symptoms
% optimal stim improvement


% 1. Participants: 'MPR8', 'MPR9', 'MPR10'

pt_ID = 'MPR8';
% pt_ID = 'MPR9';
% pt_ID = 'MPR10';

mainDir2 = [mainDir , filesep , pt_ID , filesep , 'Monopolar Sensing'];


% 2. Define hemisphere

hemisphere = 'L';
% hemisphere = 'R';

switch hemisphere
    case 'L'

        mainDir3 = [mainDir2 , filesep , 'Left STN'];   % , filesep , 'Sessions'];

    case 'R'

        mainDir3 = [mainDir2 , filesep , 'Right STN'];  % , filesep , 'Sessions'];
end

cd(mainDir3);

%% Sensing Session Data (15 streams per hemisphere, 30 streams per pt)

% i.    Structured output (.mat file)

% addpath('C:\MATLAB\GitHub\MDT_Phase2\CreateSessionStructure')
% createFolderStructure_mdt2();


%% load Structured output (.mat file) in Current Folder

% mat_name = pt_ID_hemisphere_*.mat
% load("mat_name")
% load("MPR8_L_08222023.mat")

% define pattern for structured .mat filename
filename_pattern = sprintf('%s_%s_*.mat', pt_ID, hemisphere);

% list all files in the current directory that match the pattern
pattern_match = dir(fullfile(mainDir3, filename_pattern));

% load Structured output (.mat file) in Current folder
mat_filename = fullfile(pattern_match.folder, pattern_match.name);
load(mat_filename);


%% unpack by hem, 15 total contact configs

% ii.	Tab-based structure
struct_mat = outStruct; % 1x1 struct with 4 fields (that are also 1x1 structs)
table_mat = struct2table(struct_mat); % 1x4 table with 4 fields (that are also 1x1 structs)

% % isoloate by sensing config (each contains 3-4 sub-stucts with fields)
% bipol_lev = struct2table(table_mat.BL);
% bipol_seg = struct2table(table_mat.BS);
% mono_lev = struct2table(table_mat.ML);
% mono_seg = struct2table(table_mat.MS);
%
% % S1, Info, S2, S3 (each sub-struct field contains doubles (LFP stream per contact))
% ex_stream1 = mono_seg.S1.E1;
% plot(ex_stream1)


% nested for loop
for condition_i = 1:width(table_mat) % Layer 1: condition
    sensing_config = table_mat(:,condition_i); % 1x1 table containing 1x1 stuct with 3-4 fields
    % extract each field (stream) from sensing config struct
    %
    for contact_i = 1:width() % Layer 2: contacts per condition
        % lev_or_seg =
        for stream_avg_i = 1:width() % Layer 3: average across streams
            % extract and compute average of all fields (streams) per contact per sensing config (condition)
            %
            % store average w/in that sensing config (as new field)
            %
            % compute PSD on the average of all streams per contact per sensing config (condition)
            [pxx, fxx] = pwelch(tmpMean, hanning(250), 125, 256, 250, 'onesided');  % closest approx to what MDT tablet displays (125 overlap, 256 bits, 250 window)
            uVp_t = sqrt(pxx).*rms(hanning(250)).*sqrt(2).*2.*250/256;  % method for scaling PSD --> microvolt peak output
            % PxxP(:,di2) = 10*log10(pxx);
            PxxP = pow2db(pxx); % method for scaling PSD --> decibil conversion output
        end
    end
end



% resave struct_mat with new fields (average of streams per contact per condition, fxx, pxx, uVp_t, PxxP)


%%



% save out raw avg, fxx, uVp_t, PxxP

% unsupervised ML dataset for clustering

% create a summary table
% columns
% pt_ID
% hem
% condition
% contact
% uVp_t
% PxxP




%%
% load (50 second recordings)
% Run PSD
% Use MDT variation of PSD computation


%% 3. PD symptom scores – get from Dr. Kern notes (choose clinPref)
% bsln hemibody symptoms
% optimal stim improvement


%% Filter LFP Timeseries Data

fs = 250;
nfft = 250;
window = 250;
overlap = 150;
lfpsamplerate = 2;


% field of int = average of sessions per contact

streamOfInt = ex_stream1;

% ecg = perceive_ecg(transpose(streamOfInt),250,0);
% plot(ecg.cleandata);
% yyaxis right

[px,fx,tx] = pspectrum(streamOfInt, fs, 'spectrogram');
% [pxx,fxx,txx] = pspectrum(ecg.cleandata, fs, 'spectrogram');
% outputs: power matrix, frequency vec, time vec



%% compute LFP power / instantaneous LFP beta power and plot PSDs per session (using pspectrum function)

% bin data - bin size
% narrow, 2Hz windows


% AO recordings - LFP data
% dorsal STN - 2.5mm above target
% central STN
% ventral STN
% use TTLs to pull out LFP period
% use 500ms + TTL1 as baseline
% timescale for analysis
% isolate LFP recording segments by movement context per trial
% rest/pre-move
% Hand Open/Close
% rest/transition
% Arm Pron/Sup
% rest/transition
% Elbow Flex/Extend
% rest/post-move
% kinematically derive epochs


% Hilbert transform - obtain the amplitude envelope of LFP
% study the instantaneous power of the LFP data
hilbertTransformed_LFP = hilbert(streamOfInt, fs);
% hilbertTransformed_LFP = hilbert(ecg.cleandata, fs);
amplitudeEnvelope = abs(hilbertTransformed_LFP);

% bandpass filter for beta
betaBand = [13 35]; % Beta frequency range
bpFiltered_LFP = bandpass(streamOfInt, betaBand, fs);
% bpFiltered_LFP = bandpass(ecg.cleandata, betaBand, fs);

% plot PSD to vizualize strength of the variations (energy) as a function of frequency
[pxx, fxx, txx] = pspectrum(bpFiltered_LFP, fs, 'spectrogram');
plot(fxx, pxx);
xlabel('Frequency (Hz)');
ylabel('Power');
title('Power Spectral Density');


%% Hilbert transform - obtain the amplitude envelope of LFP in the beta band.
% study the instantaneous power of the LFP data

hilbertTransformed_beta = hilbert(bpFiltered_LFP);
betaEnvelope = abs(hilbertTransformed_beta);


%%
