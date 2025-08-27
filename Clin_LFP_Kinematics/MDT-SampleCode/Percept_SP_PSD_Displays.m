clear variables;clc; close all

d = uigetdir();
files = dir(fullfile(d,'*.json'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for i = 1:length(files)
    fname = [files(i).folder,'\',files(i).name];
    if contains(fname,'BiLevL1')
        val = jsondecode(fileread(fname));
        plot_chan_data(val,1,1,'Bipolar')
    elseif contains(fname,'BiSet1L1')
        val = jsondecode(fileread(fname));
        plot_chan_data(val,1,2,'Bipolar')
    elseif contains(fname,'BiSet2')
        val = jsondecode(fileread(fname));
        plot_chan_data2(val,1,3,1)
    end
    if contains(fname,'BiLevL2')
        val = jsondecode(fileread(fname));
        plot_chan_data(val,2,1,'Bipolar')
    elseif contains(fname,'BiSet1L2')
        val = jsondecode(fileread(fname));
        plot_chan_data(val,2,2,'Bipolar')
    elseif contains(fname,'BiSet2')
        val = jsondecode(fileread(fname));
        plot_chan_data2(val,2,3,2)
    end
    if contains(fname,'MonoLevL1')
        val = jsondecode(fileread(fname));
        plot_chan_data(val,3,1,'Monopolar')
    elseif contains(fname,'MonoSegL1')
        val = jsondecode(fileread(fname));
        plot_chan_data(val,3,2,'Monopolar')
    end
    if contains(fname,'MonoLevL2')
        val = jsondecode(fileread(fname));
        plot_chan_data(val,4,1,'Monopolar')
    elseif contains(fname,'MonoSegL2')
        val = jsondecode(fileread(fname));
        plot_chan_data(val,4,2,'Monopolar')
    end
end

figure(1)
subplot(3,1,1)
title('Lead 1:Bipolar Configuration')
figure(2)
subplot(3,1,1)
title('Lead 2: Bipolar Configuration')
figure(3)
subplot(2,1,1)
title('Lead 1: Monopolar Configuration')
figure(4)
subplot(2,1,1)
title('Lead 2: Monopolar Configuration')

function plot_chan_data(val,fignum,subplotnum,config)
    channels = {};
    for x = 1:length(val.BrainSenseTimeDomain)
        data = val.BrainSenseTimeDomain(x).TimeDomainData;
        channels{x} = val.BrainSenseTimeDomain(x).Channel;
        figure(fignum)
        if contains(config,'Bipolar')
            subplot(3,1,subplotnum)
        elseif contains(config,'Monopolar')
            subplot(2,1,subplotnum)
        end
        [pxx,fxx] = pwelch(data, hanning(250), 125, 256, 250, 'onesided');   % psd estimate
        semilogy(fxx,sqrt(pxx).*rms(hanning(250)).*sqrt(2).*2.*250/256, 'LineWidth',3) % markervec{markervecindex},
        hold on
        ylabel('LFP Magnitude (uVp)')
        xlabel('Frequency (Hz)')
        legend(channels,'Interpreter','none')
        xlim([0 40])
    end
end

function plot_chan_data2(val,fignum,subplotnum,half)
    channels = {};
    for x = 1:3
        if half == 1
            xx = x;
        else
            xx = x + 3;
        end
        data = val.BrainSenseTimeDomain(xx).TimeDomainData;
        channels{x} = val.BrainSenseTimeDomain(xx).Channel;
        figure(fignum)
        subplot(3,1,subplotnum)
        [pxx,fxx] = pwelch(data, hanning(250), 125, 256, 250, 'onesided');   % psd estimate
        semilogy(fxx,sqrt(pxx).*rms(hanning(250)).*sqrt(2).*2.*250/256, 'LineWidth',3) % markervec{markervecindex},
        hold on
        ylabel('LFP Magnitude (uVp)')
        xlabel('Frequency (Hz)')
        legend(channels,'Interpreter','none')
        xlim([0 40])
    end
end
