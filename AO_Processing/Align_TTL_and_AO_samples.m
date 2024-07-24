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

% load Processed .mat files of interest
matfiles = dir('*.mat');
matnames = {matfiles.name};
for mat_names = 1:length(matnames)
    load(matnames{mat_names},'ProcEphys'); 
end

% test code within loop below with example .mat file
load("Processed_LT1D1.064F0002.mat") % 'ProcEphys'

for mat_names = 1:length(matnames)
    % isolate fields of interest
    TTL_Down = ProcEphys.TTL.Down; % number of TTLs (frames)
    TTL_clockStart = ProcEphys.TTL.startTime; % in time wrt TTL clock

    % define AO recording times                                        % Q: get these from movement index csv? (time or frame index?)
    rec_startTime = 101; % index of AO recording initiation
    rec_endTime = 1001; % index of AO recording termination

    % define TTL (frame) of interest based on task context
    TTL_taskStart = TTL_Down(:, rec_startTime); % number of samples wrt TTL clock
    TTL_taskEnd = TTL_Down(:, rec_endTime); 

    % index of task time that corresponds with timepoint of TTL
    time_offset = TTL_clockStart - rec_startTime; % time (seconds .. or ms) wrt TTL clock

    % convert time offset to sample offset
    sample_offset = round(time_offset*AO_fs); % number of samples wrt AO clock

    % Q: what's really happening here / why?
    TTL_task_idx_Start = TTL_taskStart + sample_offset; % number of samples wrt AO clock
    TTL_task_idx_End = TTL_taskEnd + sample_offset;
end


% rec_offset = 500; % duration in ms, example
% rec_endTime = rec_startTime + rec_offset; % frame index for AO recording start-time



%% go to ClustSpkTimesDir
cd(ClustSpkTimesDir)

% list of filename
SPKmatfiles = dir('*.mat');
SPKmatnames = {SPKmatfiles.name};


% test code within loop below with example .mat file
load("SpCl_E1_LT1D1.064F0002.mat")  % spikeClInfo

% storage container
all_clusts = zeros(size(SPKmatnames));

% check # of clusters
for SpCl_names = 1:length(SPKmatnames)
    load(SPKmatnames{SpCl_names},'spikeClInfo');
    temp_clust = numel(unique(spikeClInfo.clusterIDS));
    all_clusts(SpCl_names) = temp_clust;
    
    if numel(unique(spikeClInfo.clusterIDS)) ~= 1
    % determine clust IDs
    clust_IDs = unique(spikeClInfo.clusterIDS);
    for cii = 1:length(clust_IDs)
        % find unique cluster ID indices
        clust_index = ismember(spikeClInfo.clusterIDS,clust_IDs(cii));
        clust_time = spikeClInfo.SpikeTSindex(clust_index);
    end
    

    %%% NOT sure how to proceed from here %%%

    end 

end

%% use distinct clust_time vectors as input in subsequent spike cluster analysis

% load distinct spike clusters of interest


% define AO recording times
AO_startTime = spikeClInfo.AOstartTS;

% index of spike time that corresponds with timepoint of TTL
spk_time_offset = TTL_clockStart - AO_startTime; % time (in seconds .. or ms) wrt TTL clock

% convert time offset to sample offset
spk_sample_offset = round(spk_time_offset*AO_fs); % number of samples wrt AO clock


% Movement block example
TTL_spk_Start = TTL_Down(:,501);
TTL_spk_End = TTL_Down(:,1000);
TTL_spk_idx_Start = TTL_spk_Start + spk_sample_offset; % number of samples wrt AO clock
TTL_spk_idx_End = TTL_spk_End + spk_sample_offset;


SpikeTS = spikeClInfo.SpikeTSindex; % spikes in samp # wrt AO clock
spikes_in_move1 = SpikeTS > TTL_spk_idx_Start & SpikeTS < TTL_spk_idx_End; % logical of spikes w/in moveblock in AO_time


%% movement block structure
% define start frame of event
% define end frame of event
% define offset duration before & after event


%% initial code

lastSpikeTS = spikeClInfo.SpikeTSindex(end);
numSpikes = spikeClInfo.SpikeTSindex(end) - spikeClInfo.SpikeTSindex(1);
num_seconds = numSpikes/AO_fs;

% waveforms of spikes only necessary for figures
waves = spikeClInfo.waveForms.waves; 

