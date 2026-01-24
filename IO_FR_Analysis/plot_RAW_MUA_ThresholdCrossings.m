

ProcEphys_dir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\07_13_2023_bilateral\Processed Electrophysiology\LSTN';
cd(ProcEphys_dir)
load('Processed_LT1D4.071F0002.mat')

ClustSpk_dir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\07_13_2023_bilateral\Processed Electrophysiology\LSTN\ClusteredSpikeTimes';
cd(ClustSpk_dir)
load('SpCl_E1_LT1D4.071F0002.mat')
plot(ProcEphys.Spike.E1.rawData)
plot(ProcEphys.Spike.E1.rawData(1:88000))
 
spkTrim = spikeClInfo.SpikeTSindex(spikeClInfo.SpikeTSindex < 88000);

hold on;
plot(spkTrim, ProcEphys.Spike.E1.rawData(spkTrim), '*r')

%%


% loop through first 3 seconds 
ProcE_list = dir('*.mat');
ProcE_list2 = struct2table(ProcE_list);
ProcE_list3 = ProcE_list2.name;

close all;

for file_i = 1:length(ProcE_list3)
    load(ProcE_list3{file_i}) 
    plot(ProcEphys.Spike.E1.rawData)
    plot(ProcEphys.Spike.E1.rawData(1:88000))
    title(ProcE_list3{file_i});

    pause
    close all;

end

