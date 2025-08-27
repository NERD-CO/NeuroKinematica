%% IO data analysis outline


%% Analog Front-End Signals

% 1) MER (using AO Neuro Omega electrode channels)

  % Micro Contact
    % f_micro_spike = 300 - 1000 Hz
    % f_micro_LFP = 0-200 Hz
  % Macro Contact
    % f_macro_LFP = 0-200 Hz

    SR_MER_Raw = 44 kHz
    SR_MER_LFP = 1.375 kHz
    Gain_MER = 20

% 2) Video (using two cameras (Blackfly USB 3.0) + FLIR-Spinnaker-Python software + custom python GUI)

    % f_movement = 0-10 Hz

    SR_movement = 100 fps
    Exposure = 5000

% 3) ACC (using AO Neuro Omega Analog Inputs - SensorBox)

    % f_movement = 0-10 Hz

    SR_ACC = 2.75 kHz
    Gain_ACC = 0.25

% 4) sEMG (using AO Neuro Omega Headbox module)

    % f_EMG = 0-1000 Hz

    SR_EMG = 44 kHz
    Gain_EMG = 55

%% Clock - TTL

% CDIG (Hz, Down, TimeBegin, TimeEnd)
f_CDIG = 44 kHz

%% 

% > trial_ID_determination
% > save_IO_EphysProcFiles
    % contains functions (2)
    % > [mat_filelist, ACC_check] = save_IO_mat_filenames(studyID, IO_DataDir); (1)
    % > save_IO_mat_ProcFiles(mat_filelist, Case_DataDir, ACC_check); (2)
% save_IO_mat_ProcFiles
    % 'CSPK','CLFP','CMacro_LFP', 'CEMG', 'CACC' (+subfields)

%% Preprocessing

% Filter
% Downsample? No 
% Segment (into epochs)
    % by trial type
% Feature Extraction
    % PSD - quant. power distrib
    % STFT (time-frequency analysis)
    % Hilbert
    % envelope function

%% beta bursting pipeline


%% AO recordings - LFP data
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