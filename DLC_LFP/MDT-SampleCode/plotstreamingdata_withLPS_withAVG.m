clear variables; clc; close all
fs=250;
nfft=250;
window=250;
overlap=200;

% Import and read file
[val1,path] = Import_Read_first;
[val2,path] = Import_Read_second(path);
[val3,path] = Import_Read_third(path);

%remove first 15 seconds
figure(1)
figure(2)
figure(3)

i = length(val1.BrainSenseTimeDomain);
channels = {};


%% Plot into figure 1-3
for x = 1:i
    data = val1.BrainSenseTimeDomain(x).TimeDomainData;
    channels{x} = val1.BrainSenseTimeDomain(x).Channel;
    
    figure(1)
    subplot(i,1,x)
    plot(data);
    
    figure(2)
    subplot(i,1,x)
    [s_L,f_L,t_L, P_L]=spectrogram(data*.015,window,overlap,nfft,fs);
    imagesc(t_L,f_L,10*log10(abs(P_L)+eps));
    set(gca,'YDir','normal');
    ylabel('Freq (Hz)')
    xlabel('Time (s)')
    colormap(gca,'jet')
    caxis([-100 35])
    colorbar off
    ylim([0,100])
    
    figure(3)
    [pxx,fxx] = pwelch(data, hanning(250), 125, 256, 250, 'onesided');   % psd estimate
    semilogy(fxx,sqrt(pxx).*rms(hanning(250)).*sqrt(2).*2.*250/256, 'LineWidth',3) % markervec{markervecindex},
    hold on
    ylabel('LFP Magnitude (uVp)')
    xlabel('Frequency (Hz)')
    legend(channels)
        
end

%% Plot into figure 4-6
for x = 1:i
    data = val2.BrainSenseTimeDomain(x).TimeDomainData;
    channels{x} = val2.BrainSenseTimeDomain(x).Channel;
    
    figure(4)
    subplot(i,1,x)
    plot(data);
    
    figure(5)
    subplot(i,1,x)
    [s_L,f_L,t_L, P_L]=spectrogram(data*.015,window,overlap,nfft,fs);
    imagesc(t_L,f_L,10*log10(abs(P_L)+eps));
    set(gca,'YDir','normal');
    ylabel('Freq (Hz)')
    xlabel('Time (s)')
    colormap(gca,'jet')
    caxis([-100 35])
    colorbar off
    ylim([0,100])
    
    figure(6)
    [pxx,fxx] = pwelch(data, hanning(250), 125, 256, 250, 'onesided');   % psd estimate
    semilogy(fxx,sqrt(pxx).*rms(hanning(250)).*sqrt(2).*2.*250/256, 'LineWidth',3) % markervec{markervecindex},
    hold on
    ylabel('LFP Magnitude (uVp)')
    xlabel('Frequency (Hz)')
    legend(channels)
        
end

%% Plot into figure 7-9
for x = 1:i
    data = val3.BrainSenseTimeDomain(x).TimeDomainData;
    channels{x} = val3.BrainSenseTimeDomain(x).Channel;
    
    figure(7)
    subplot(i,1,x)
    plot(data);
    
    figure(8)
    subplot(i,1,x)
    [s_L,f_L,t_L, P_L]=spectrogram(data*.015,window,overlap,nfft,fs);
    imagesc(t_L,f_L,10*log10(abs(P_L)+eps));
    set(gca,'YDir','normal');
    ylabel('Freq (Hz)')
    xlabel('Time (s)')
    colormap(gca,'jet')
    caxis([-100 35])
    colorbar off
    ylim([0,100])
    
    figure(9)
    [pxx,fxx] = pwelch(data, hanning(250), 125, 256, 250, 'onesided');   % psd estimate
    semilogy(fxx,sqrt(pxx).*rms(hanning(250)).*sqrt(2).*2.*250/256, 'LineWidth',3) % markervec{markervecindex},
    hold on
    ylabel('LFP Magnitude (uVp)')
    xlabel('Frequency (Hz)')
    legend(channels)
        
end

%% Average PSD

%Array for colors of lines
color1 = [0.00,0.45,0.74; 0.85,0.33,0.10; 0.93,0.69,0.13; 0.49,0.18,0.56; 0.47,0.67,0.19; 0.30,0.75,0.93];
for x = 1:i
    channels{x} = val1.BrainSenseTimeDomain(x).Channel;
    data1 = val1.BrainSenseTimeDomain(x).TimeDomainData;
    data2 = val2.BrainSenseTimeDomain(x).TimeDomainData;
    data3 = val3.BrainSenseTimeDomain(x).TimeDomainData;
    
    figure(10)
    [pxx1,fxx1] = pwelch(data1, hanning(250), 125, 256, 250, 'onesided');   % psd estimate
    [pxx2,fxx2] = pwelch(data2, hanning(250), 125, 256, 250, 'onesided');   % psd estimate
    [pxx3,fxx3] = pwelch(data3, hanning(250), 125, 256, 250, 'onesided');   % psd estimate
    p1 = sqrt(pxx1).*rms(hanning(250)).*sqrt(2).*2.*250/256;
    p2 = sqrt(pxx2).*rms(hanning(250)).*sqrt(2).*2.*250/256;
    p3 = sqrt(pxx3).*rms(hanning(250)).*sqrt(2).*2.*250/256;
    semilogy(fxx,mean(vertcat(p1',p2',p3')), 'LineWidth',3) % markervec{markervecindex},
    hold on
    
    
    curve1 = max(vertcat(p1',p2',p3'));
    curve2 = min(vertcat(p1',p2',p3'));
    s = patch([fxx1' fliplr(fxx1')], [curve2 fliplr(curve1)], color1(x,:));
    alpha(s,.2);
    
    
    ylabel('LFP Magnitude (uVp)')
    xlabel('Frequency (Hz)')
%     legend(channels)
end

%% Import Read function
function [val,path] = Import_Read_first()
%Import_Read Lets user select .json file for use in the Field Support
%Dashboard
%   This function allows the user to select a .json file and decodes it
%   into a structure labeled 'val' which is then used to populate the Field
%   Support Dashboard
% [file,path] = uigetfile('*.json','Select .json file',...

[file,path] = uigetfile('*.json','Select first .json file','C:\Medtronic\Medtronic CORTEX\spstudio\perceptsp\controlapp\Subjects');
fname = [path,file];
val = jsondecode(fileread(fname));

end

%% Import Read function (follow-up)
function [val,path] = Import_Read_second(path)
%Import_Read Lets user select .json file for use in the Field Support
%Dashboard
%   This function allows the user to select a .json file and decodes it
%   into a structure labeled 'val' which is then used to populate the Field
%   Support Dashboard
% [file,path] = uigetfile('*.json','Select .json file',...

[file,path] = uigetfile('*.json','Select second .json file',path);
fname = [path,file];
val = jsondecode(fileread(fname));

end

%% Import Read function (follow-up)
function [val,path] = Import_Read_third(path)
%Import_Read Lets user select .json file for use in the Field Support
%Dashboard
%   This function allows the user to select a .json file and decodes it
%   into a structure labeled 'val' which is then used to populate the Field
%   Support Dashboard
% [file,path] = uigetfile('*.json','Select .json file',...

[file,path] = uigetfile('*.json','Select third .json file',path);
fname = [path,file];
val = jsondecode(fileread(fname));

end