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

mainDir2 = [mainDir , filesep , pt_ID , 'MonopolarSensing'];

% 2. Define hemisphere

hemisphere = 'Left';
%hemisphere = 'Right';

switch hemisphere
    case 'Left'

        mainDir3 = [mainDir2 , filesep , 'Left' , filesep, 'Sessions'];

    case 'Right'

        mainDir3 = [mainDir2 , filesep , 'Right' , filesep, 'Sessions'];
end


%% Sensing Session Data (15 streams per hemisphere, 30 streams per pt)

% i.    Structured output (.mat file)

addpath('C:\MATLAB\GitHub\MDT_Phase2\CreateSessionStructure')
createFolderStructure_mdt2();

% ii.	Tab-based structure


%% unpack by hem, 15 total contact config

% load("MPR8_L_08222023.mat")

% Rename Session######### in 'Sessions' to 'Session'
% struct_mat = [mainDir3 , filesep, 'Sessions' , filesep, 'Session' , filesep, 'StructureFile'];

struct_mat = outStruct;

% isoloate by sensing config (each contains 3-4 sub-stucts with fields)
bipol_lev = struct2table(struct_mat.BL);
bipol_seg = struct2table(struct_mat.BS);
mono_lev = struct2table(struct_mat.ML);
mono_seg = struct2table(struct_mat.MS);

% S1, Info, S2, S3 (each sub-struct field contains doubles (LFP stream per contact))
ex_stream1 = mono_seg.S1.E1;

plot(ex_stream1)


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
