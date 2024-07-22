%% Environment / Data Directory Inputs

% specify directory where case-specific data files are located 
curPCname = getenv('COMPUTERNAME'); 

switch curPCname
    case 'DESKTOP-I5CPDO7'  % PC_1
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';  
    case 'DSKTP-JTLAB-EMR'  % Lab Desktop
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative'; 
    case 'NSG-M-H8J3X34'    % PC_2
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
end

% isolate a specific CaseDate / studyID (in Subject_AO csv)

% CaseDate = '03_09_2023'; % studyID = 1 
CaseDate = '03_23_2023'; % studyID = 2 

Case_DataDir = [IO_DataDir, filesep, CaseDate];

%% Hardcode directories

% directory where all IO data is located
RawDataDir = [Case_DataDir, filesep, 'Raw Electrophysiology MATLAB'];       % directory where raw MATLAB data files are located (case-specific)
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology'];       % directory where processed MATLAB data should be saved (case-specific)
ClustSpkTimesDir = [ProcDataDir, filesep, 'ClusteredSpikeTimes'];           % directory where clustered spike times should be saved (case-specific)

% define datastream sampling rates (Alpha Omega and Video sampling fs)
AO_fs = 44000; % Hz
DLC_fs = 1000; % fps

%% go to ProcDataDir
cd(ProcDataDir)

% turn into loop later
% load Processed .mat files of interest
load("Processed_LT1D1.064F0002.mat")

% isolate fields of interest
TTL_Down = ProcEphys.TTL.Down; % number of TTLs in frames
TTL_startTime = ProcEphys.TTL.startTime;

% define TTL of interest based on task start time
TTL_of_int = TTL_Down(:,101); % samplele # in AO


%% go to ClustSpkTimesDir
cd(ClustSpkTimesDir)

% turn into loop later
% load spike cluster(s) of interest
load("SpCl_E1_LT1D1.064F0002.mat")
AO_startTime = spikeClInfo.AOstartTS;

% index of spike time that corresponds with timepoint of TTL
time_offset = TTL_startTime - AO_startTime; % in seconds

% convert time offset to samples
sample_offset = round(time_offset*AO_fs); % in samples wrt to AO clock

% 
TTL_spike_idx = TTL_of_int + sample_offset; % in samples wrt to AO clock

% Movement block example
TTL_of_int_Start = TTL_Down(:,501);
TTL_of_int_End = TTL_Down(:,1000);
TTL_spike_idx_Start = TTL_of_int_Start + sample_offset;
TTL_spike_idx_End = TTL_of_int_End + sample_offset;

SpikeTS = spikeClInfo.SpikeTSindex; % spikes in samp # wrt AO clock
spikes_in_move1 = SpikeTS > TTL_spike_idx_Start & SpikeTS < TTL_spike_idx_End; % logical of spikes w/in moveblock in AO_time



%%
% check # of clusters

%
lastSpikeTS = spikeClInfo.SpikeTSindex(end);
numSpikes = spikeClInfo.SpikeTSindex(end) - spikeClInfo.SpikeTSindex(1);
num_seconds = numSpikes/AO_fs;

% waveforms of spikes only necessary for figures
waves = spikeClInfo.waveForms.waves; 

