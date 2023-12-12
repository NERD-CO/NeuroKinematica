% MEDTRONIC CONFIDENTIAL
% 
% DO NOT DISTRIBUTE
% 
% RESEARCH USE ONLY

% To be shared only within Medtronic and with UCH researchers
% NOT VALIDATED, FOR RESEARCH USE IN MDT-UCH COLLABORATIVE STUDY A 1718256

function UCH_PowerSnapTD_short(val,path)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Support Variables Identifying Channels
num_stream = length(val.BrainSenseTimeDomain);
conc_ch = cell(1,num_stream);

for k = 1:num_stream
    conc_ch{k} = val.BrainSenseTimeDomain(k).Channel;
end

order_ch = {'ZERO_TWO_LEFT','ONE_THREE_LEFT','ZERO_THREE_LEFT',...
    'ZERO_TWO_RIGHT','ONE_THREE_RIGHT','ZERO_THREE_RIGHT'};
for k = 1:6
    index{k} = find(contains(conc_ch,order_ch{k}));
end

%% Testing TicksInMses update (2.0.4584-2.0.4594)
isoldformat = false;

for chan = 1:6
    if length(index{chan})>=2
        TickS1 = (str2num([val.BrainSenseTimeDomain(index{chan}(1)).TicksInMses]));
        TickS2 = (str2num([val.BrainSenseTimeDomain(index{chan}(2)).TicksInMses]));
        if TickS1(1)==TickS2(1)
            isoldformat = true;
            break
        end
    end
end

%% Time array (accounting dropped packets)

if isoldformat % Old format requires special handling due to TicksInMses duplication
    init_tick = min(str2num([val.BrainSenseTimeDomain(1).TicksInMses])); % intial streaming tick for session
    for k = 1:6 % For each possible sensing configuration
        if not(isempty(index{k}))
            GPsize = str2num([val.BrainSenseTimeDomain(index{k}).GlobalPacketSizes]);
            TickS = (str2num([val.BrainSenseTimeDomain(index{k}(1)).TicksInMses])-init_tick)/1000;

            TDtime = TickS(end)-GPsize(end)/250:1/250:TickS(end)-1/250;
            for i = length(GPsize):-1:2
                if TickS(i)-TickS(i-1) > GPsize(i)/250 + 1/250
                    Prev_packet = TickS(i-1)-GPsize(i-1)/250:1/250:TickS(i-1)-1/250;
                    TDtime = [Prev_packet,TDtime];
                else
                    Prev_packet = TDtime(1)-GPsize(i-1)/250:1/250:TDtime(1)-1/250;
                    TDtime = [Prev_packet,TDtime];
                end
            end
            index_TDtime{k} = TDtime;
            clear TDtime
        else
            index_TDtime{k} = [];
        end
    end
else % New format does not contain TicksInMses duplication
    index_TDtime = cell(1,6);
    for k = 1:num_stream % For each streaming sample
        GPsize = str2num([val.BrainSenseTimeDomain(k).GlobalPacketSizes]);
        init_tick = min(str2num([val.BrainSenseTimeDomain(k).TicksInMses])); % intial streaming tick for streaming sample
        if k == 1
            TickS = (str2num([val.BrainSenseTimeDomain(k).TicksInMses])-init_tick)/1000;
        else
            StartTime(1) = datetime(val.BrainSenseTimeDomain(1).FirstPacketDateTime,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z');
            StartTime(2) = datetime(val.BrainSenseTimeDomain(k).FirstPacketDateTime,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z');
            diffStartTime = seconds(diff(StartTime));
            TickS = (str2num([val.BrainSenseTimeDomain(k).TicksInMses])-init_tick+1000*diffStartTime)/1000;
        end
        TDtime = TickS(end)-GPsize(end)/250:1/250:TickS(end)-1/250;
        for i = length(GPsize):-1:2
            if TickS(i)-TickS(i-1) > GPsize(i)/250 + 1/250
                Prev_packet = TickS(i-1)-GPsize(i-1)/250:1/250:TickS(i-1)-1/250;
                TDtime = [Prev_packet,TDtime];
            else
                Prev_packet = TDtime(1)-GPsize(i-1)/250:1/250:TDtime(1)-1/250;
                TDtime = [Prev_packet,TDtime];
            end
        end
        placement = find(contains(order_ch,val.BrainSenseTimeDomain(k).Channel));
        index_TDtime{placement} = [index_TDtime{placement},TDtime];
        clear TDtime
    end
end

    
%% Channel TD data array

for k = 1:6 % For each possible sensing configuration
    data_cat{k}=vertcat(val.BrainSenseTimeDomain(index{k}).TimeDomainData);
end

%% Plotting TD

channel = {['0-2L'],['1-3L'],['0-3L'],['0-2R'],['1-3R'],['0-3R']};

color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};

if not(isempty(vertcat(data_cat{1:3})))
    figure(301)
    subplot(3,1,1)
    for i = 1:3
        if not(isempty(data_cat{i}))
            hold on
            plot(index_TDtime{i},data_cat{i},'Linestyle',':',...
                'Marker','.','MarkerSize',3,'DisplayName',channel{i},...
                'color',color{i})
        end
    end
    xlim([min([index_TDtime{1:6}]),max([index_TDtime{1:6}])])
    xlabel('Time (s)')
    ylabel('Volt (uV)')
    legend
    title('Left Hemisphere')
    
%     figure(1)
%     subplot(6,3,[10,11])
%     for i = 1:3
%         if not(isempty(data_cat{i}))
%             hold on
%             plot(index_TDtime{i},data_cat{i},'Linestyle',':',...
%                 'Marker','.','MarkerSize',3,'DisplayName',channel{i},...
%                 'color',color{i})
%         end
%     end
%     xlim([min([index_TDtime{1:6}]),max([index_TDtime{1:6}])])
%     xlabel('Time (s)')
%     ylabel('Volt (uV)')
%     legend
%     title('Streaming (Left)')
end

if not(isempty(vertcat(data_cat{4:6})))
    figure(302)
    subplot(3,1,1)
    for i = 4:6
        if not(isempty(data_cat{i}))
            hold on
            plot(index_TDtime{i},data_cat{i},'Linestyle',':',...
                'Marker','.','MarkerSize',3,'DisplayName',channel{i},...
                'color',color{i-3})
        end
    end
    xlim([min([index_TDtime{1:6}]),max([index_TDtime{1:6}])])
    xlabel('Time (s)')
    ylabel('Volt (uV)')
    legend
    title('Right Hemisphere')
    
%     figure(2)
%     subplot(6,3,[10,11])
%     for i = 4:6
%         if not(isempty(data_cat{i}))
%             hold on
%             plot(index_TDtime{i},data_cat{i},'Linestyle',':',...
%                 'Marker','.','MarkerSize',3,'DisplayName',channel{i},...
%                 'color',color{i-3})
%         end
%     end
%     xlim([min([index_TDtime{1:6}]),max([index_TDtime{1:6}])])
%     xlabel('Time (s)')
%     ylabel('Volt (uV)')
%     legend
%     title('Streaming (Right)')
end

%% Spectrogram & Plotting

% fs=250;
% nfft=64;
% window=250;
% overlap=60;
fs=250;
nfft=250;
window=250;
overlap=150;

index_TDtimeL = [index_TDtime{1:3}];
index_TDtimeR = [index_TDtime{4:6}];

data_catL = vertcat(data_cat{1:3});
data_catR = vertcat(data_cat{4:6});

[BL,IL] = sort(index_TDtimeL);
[BR,IR] = sort(index_TDtimeR);


spectro_catL = data_catL(IL);
for i = length(spectro_catL):-1:2
    if BL(i)-BL(i-1)>1/249
        add_sampl = zeros(round((BL(i)-BL(i-1))*250),1);
        spectro_catL(i:end+length(add_sampl)) = ...
            vertcat(add_sampl,spectro_catL(i:end));
        clear add_sampl
    end
end
if not(isempty(BL))&BL(1) >= 0
    add_sampl = zeros(round((BL(1)-0)*250),1);
    spectro_catL = vertcat(add_sampl,spectro_catL);
    clear add_sampl
end
spectro_catR = data_catR(IR);
for i = length(BR):-1:2
    if BR(i)-BR(i-1)>1/249
        add_sampl = zeros(round((BR(i)-BR(i-1))*250),1);
        spectro_catR(i:end+length(add_sampl)) = ...
            vertcat(add_sampl,spectro_catR(i:end));
        clear add_sampl
    end
end
if not(isempty(BR))& BR(1) >= 0
    add_sampl = zeros(round((BR(1)-0)*250),1);
    spectro_catR = vertcat(add_sampl,spectro_catR);
    clear add_sampl
end

endfill = abs(length(spectro_catL)-length(spectro_catR));
if endfill>0
    if length(spectro_catL)<length(spectro_catR)
        spectro_catL = vertcat(spectro_catL, zeros(endfill,1));
    else
        spectro_catR = vertcat(spectro_catR, zeros(endfill,1));
    end
end


% Spectrogram of Time Domain Data
if not(isempty(data_catL))
    figure(301)
    subplot(3,1,2)
    % spectrogram(spectro_catL,window,overlap,nfft,fs,'yaxis')
    % ylabel('Norm Freq')
    [s_L,f_L,t_L, P_L]=spectrogram(spectro_catL,window,overlap,nfft,fs);
    imagesc(t_L,f_L,10*log10(abs(P_L)+eps));
    set(gca,'YDir','normal');
    ylabel('Freq (Hz)')
    xlabel('Time (s)')
    colormap(gca,'jet')
    caxis([-100 35])
    colorbar off

end
if not(isempty(data_catR))
    figure(302)
    subplot(3,1,2)
    %spectrogram(spectro_catR,window,overlap,nfft,fs,'yaxis')
    %ylabel('Norm Freq')
    [s_R,f_R,t_R, P_R]=spectrogram(spectro_catR,window,overlap,nfft,fs);
    imagesc(t_R,f_R,10*log10(abs(P_R)+eps));
    set(gca,'YDir','normal');
    ylabel('Freq (Hz)')
    xlabel('Time (s)')
    colormap(gca,'jet')
    caxis([-100 35])
    colorbar off
end


end