%% Kinematic pipeline notes

%% 1) Compress raw videos (on flashdrive)
> compressVideos (RW-provided function in 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\DLC_Processing')
- or -
> compressVideosJAT (JAT-modified function in 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\DLC_VideoCheck_GUI')
Q: What is the difference?


%% Add step: 
% 2) Ensure all study vids are present, compressed, and converted to MP4
** check for missing videos or videos that are uncompressed AVI vs. compressed MP4
** create generic function from one-off script 
    - video_filetype_conversion (in 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\DLC_Processing')


%% Data transfer & organization
- upload compressed videos to 'Compressed Raw Video' subfolder within case folder in Box
- after case folder is populated with AO / LFP data, download full case folder in Box to lab desktop 
[in 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative' or 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical']


%% DLC video labelling & model training
[in 'C:\Users\erinr\Desktop\DLC']
- create case video folder for dlc model and copy case videos (must be compressed & converted to MP4; if not, return to steps 1-3)
* (TIME): Label case videos/video pairs in dlc                                      
- train dlc model, evaluate network, analyze videos
- create dlc-labelled videos


%% DLC Processing
[in 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative' or 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical']
- create case folder in outer 'post_DLCproc' folder 
    - contains 'video folder', 'csv folder', 'mat folder'; (video folder contains 'Converted for dual GUI' subfolder)
        - populate 'video folder' with dlc-labelled case video pairs
        - populate 'csv folder' with dlc timeseries labelling outputs
> run_dlc_processCSV2 - convert dlc timeseries labelling outputs per video from .csv to .mat


%% Movement Indexing 
[in 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative' or 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical']
[within 'post_DLCproc' > [case folder] > 'video folder' subfolder]
- convert case video pairs (for dual movecheck GUI)
- save [within 'post_DLCproc' > [case folder] > 'video folder' > 'Converted for dual GUI' subfolder]
* (TIME): Movement index dlc-labelled videos/video pairs (via dual movecheck GUI)


** create generic functions from one-off script
> preGUI_vidConversion (in 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\DLC_VideoCheck_GUI')

** (?) Unfamiliar functions / scripts in 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\DLC_VideoCheck_GUI'
> compressVideosJAT
> fillDroppedFrames_JT_v2
> preGUI_vidConversion



%% Kinematic Analyses

% Scripts to run necessary functions in order
> Order_of_Operations__Clin_kinematics
-or
> Order_of_Operations__IO_kinematics

    % Process + vizualize kinematic data
    > run_MovementProcessing_Clin_v1
    -or-
    > run_MovementProcessing_IO_v1
    o	Iterative sanity checking per trial
    o	artifactRejection

    % Outputs [in 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\[Clin or IO]\Kinematic Analyses\[casedate]\[hem]\processedMovement']
        % findpeaks results + plot
        % cleaned dlcDAT data


    
    % Statistically summarize & compare movement timeseries data across conditions (IO depths -or- Clin states)
    > run_MovementStats_Clin_v1
    -or-
    > run_MovementStats_IO_v1
    
    o	Random data trimming for t-test and ANOVA 
    o   Descriptive stats: mean, stdev, variance
    o   kruskalwallis vs. anova1 - non-parametric comparison
    
    functions
    - compute_ConditionStats
    - summarize_ConditionStats
    - ttest_plot_2States_nonRandTrim
    - ttest_plot_2States_RandTrim
    - rand_trim_dataset
    - ANOVA_plot_allStates (**multiple variations - which is appropriate?)

    o	Outputs [in 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\[Clin or IO]\Kinematic Analyses\[casedate]\[hem]\movementStats']


** create generic function
%% generic function goal

% inputs:
    % case date
    % hemisphere
    % movement type (Hand O/C, Pron/Sup, Elbow F/E, Finger Tap)
    % other variables: conditiion (clinical) or depth (IO)
    % raw dlc label data - csv or mat (outDATA)
    % artifact flag: use or don't use artifactRejection function

% outputs:
    % plot of raw data per dlc marker of interest
    % plot of interpolated / cleaned data: outData interp
    % decision per frame table - binary status per marker per frame (accept/reject): outData Index
    % interpolated / cleaned data table: outData interp
    % movement analyses based on function inputs


% Develop algorthmic functions based on movment type and dlc-labelled markers of interest
    % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9883391/
    % Hand Movement (Hand O/C): area of the convex hull (ACH) of the four finger tips key-points and the palm key-point, measured in units of estimated-standing-height squared (H2).
    % Pronation-Supination (Pron/Sup): angular velocity of the vector from the thumb-tip key-point to the little-finger-tip key-points, measured in degrees per frame.
    % Finger Tapping (Finger Tap): Euclidean distance between the thumb-tip key-point and the index finger tip key-point, measured in units of estimated-standing height.
    % Find/develop definition/method/algorthm for Elbow Flex/Extend (Elbow F/E)





%% Analysis goals/feasibility check

For 12/18:

-   Clinical Movement analysis - Hand Open/Close  (Bonus: FingerTap) 
o	'Clin_2023-09-12_LSTN'
o	'Clin_2023-09-12_RSTN'


-	Clinical LFP analysis 
o	'Clin_2023-09-12_LSTN' % target first
o	'Clin_2023-09-12_RSTN'


-	IO Movement analysis - Hand Open/Close
o	'IO_2023-03-09' (RSTN): 3 depths  % target first
o	'IO_2023-05-11' (LSTN): 2 depths

-	IO LFP analysis
o	'IO_2023-03-09' (RSTN): 3 depths  % target first
o	'IO_2023-05-11' (LSTN): 2 depths





%% lab meeting notes_20231204
lfp power
beta filter
- Hilbert Transform
bursting detection
bursting threshold
burst frequency
burst duration

spike/MER data
envelope function?
waveform shape

artifact/outlier rejection
smoothing
filter

naive bayes classifier, labels