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



%% Kinematic Analysis
> run_MovementProcessing_v1 (Processing + Viz)
o	artifactRejection
o	Iterative sanity checking per trial

> run_MovementStats_v1
o	Random data trimming for t-test and ANOVA 
o   Descriptive stats: mean, stdev, variance
o   kruskalwallis vs. anova1 - non-parametric comparison
functions
- compute_ConditionStats
- summarize_ConditionStats
- ttest_plot_2States_nonRandTrim
- ttest_plot_2States_RandTrim
- rand_trim_dataset
- ANOVA_plot_allStates (multiple variations - which is appropriate?)


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