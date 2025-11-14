%% Spike Cluster Quality Reviewer ppt output

% Run this in pre-2025 MATLAB version

%% Environment / Directory Set-up

% specify directory where case-specific data files are located
curPCname = getenv('COMPUTERNAME');

switch curPCname
    case 'DESKTOP-I5CPDO7'  % PC_1
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
    case 'DSKTP-JTLAB-EMR'  % Lab Desktop
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        GitHub_RepoDir = 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis';
        addpath 'C:\Users\erinr\OneDrive - The University of Colorado Denver\Documents 1\GitHub\NeuroKinematica\IO_FR_Analysis'
    case 'NSG-M-H8J3X34'    % PC_2
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';
        GitHub_RepoDir = 'C:\GitHub\NeuroKinematica\IO_FR_Analysis';
        addpath 'C:\GitHub\NeuroKinematica\IO_FR_Analysis'
end

cd(IO_DataDir)
Subject_AO = readtable('Subject_AO.xlsx');
Subject_Hem_CaseMap = readtable('Subject_Hem_MetaSummary.xlsx'); % Subject_Hem_MetaSummary.xlsx


%% Next step (replace CaseDate input block below):

% Loop through the 'Subject_AO' or 'Subject_Hem_CaseMap' CaseFolder rows
% (correspond to CaseDate) and perform this processing, i.e., 
% run the BatchSpikeClusterQuality function on each 
% Spike Cluster Directory per CaseFolder 


%% Input: Case-specific Ephys data dir

CaseDate = '03_23_2023';

% '03_23_2023';             % NER 2025, NANS 2026, INS 2026      % 1
% '04_05_2023';             % NER 2025, NANS 2026, INS 2026      % 1
% '04_13_2023_bilateral';                        % GRC 2026      % 2 (LSTN)
% '05_18_2023_b_bilateral'; % NER 2025, NANS 2026, INS 2026      % 3 (LSTN) %% error in indice bounds --> check movement indices
% '05_31_2023';                                  % INS 2026      % 4
% '06_08_2023_bilateral';   % NER 2025, NANS 2026, INS 2026      % 5 (LSTN) %% check movement indices
% '07_06_2023_bilateral';                        % INS 2026      % 6 (LSTN)
% '07_13_2023_bilateral';                        % INS 2026      % 7 (LSTN, RSTN)
% '08_23_2023';                       % NANS 2026, INS 2026      % 8
% '11_30_2023_bilateral';                        % GRC 2026      % 9 (LSTN, RSTN)


%% Define case-specific data input and outputs directories

Ephys_CaseDir = [IO_DataDir, filesep, CaseDate];       % case-specific ephys data input directory 

% Define directories where case-specific IO ephys data are located (inputs)
ProcDataDir = [Ephys_CaseDir, filesep, 'Processed Electrophysiology'];       % directory where processed ephys data should be saved (case-specific)


%% Handle bilateral cases and hemisphere selection

isBilateral = contains(CaseDate, 'bilateral', 'IgnoreCase', true);

if isBilateral
    fprintf('[INFO] Bilateral case detected: %s\n', CaseDate);

    % Prompt user for hemisphere choice (LSTN or RSTN)
    CaseDate_hem = input('Enter hemisphere (LSTN or RSTN): ', 's');

    % Validate input
    validHems = {'LSTN','RSTN'};
    if ~ismember(CaseDate_hem, validHems)
        error('Invalid input. Please enter LSTN or RSTN.');
    end
else
    CaseDate_hem = ''; % No hemisphere for unilateral cases
end

% Append hemisphere folder if needed
if ~isempty(CaseDate_hem)
    ProcDataDir = fullfile(ProcDataDir, CaseDate_hem);
    fprintf('[INFO] Hemisphere-specific input ephys directory set: %s\n', ProcDataDir);
else
    fprintf('[INFO] Using base input directory: %s\n', ProcDataDir);   
end


%% Define Clustered Spike Times directory

Cluster_SpkTimesDir = fullfile(ProcDataDir, 'ClusteredSpikeTimes'); % directory where clustered spike times should be saved (case-specific)
if ~isfolder(Cluster_SpkTimesDir)
    error('[ERROR] ClustSpkTimesDir does not exist: %s', Cluster_SpkTimesDir);
end
cd(Cluster_SpkTimesDir);


%% Run BatchSpikeClusterQuality review function 
% outputs 

BatchSpikeClusterQuality(Cluster_SpkTimesDir)