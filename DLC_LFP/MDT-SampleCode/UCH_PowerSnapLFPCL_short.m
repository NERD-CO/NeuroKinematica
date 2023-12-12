% MEDTRONIC CONFIDENTIAL
% 
% DO NOT DISTRIBUTE
% 
% RESEARCH USE ONLY

% To be shared only within Medtronic and with UCH researchers
% NOT VALIDATED, FOR RESEARCH USE IN MDT-UCH COLLABORATIVE STUDY A 1718256

function PowerSnapLFPCL_short(val,path)
% Function to create figure containing sense configuration specific
% subplots
%   updated 8/9/2020

%% Support Variables Identifying Channels
num_stream = length(val.BrainSenseLfp);
conc_ch = cell(1,num_stream);

for k = 1:num_stream
    conc_ch{k} = val.BrainSenseLfp(k).Channel;
end

order_ch = {'ZERO_TWO_LEFT','ONE_THREE_LEFT','ZERO_THREE_LEFT',...
    'ZERO_TWO_RIGHT','ONE_THREE_RIGHT','ZERO_THREE_RIGHT'};
index = cell(1,6);
for k = 1:6
    index{k} = find(contains(conc_ch,order_ch{k}));
end


%% Time array (accounting dropped packets)

TDLFP=[];LDoutL=[];LDoutR=[];StimoutL=[];StimoutR=[];
init_tick = val.BrainSenseLfp(1).LfpData(1).TicksInMs; % intial streaming tick for session

for k = 1:num_stream % For each streaming sample
    
    TDx = [];
    LFPL = [];
    LFPR = [];
    StimL = [];
    StimR = [];
    
    TDx = horzcat(TDx,[val.BrainSenseLfp(k).LfpData.TicksInMs]);
    TDx = TDx - TDx(1);
    for j = 1:length(val.BrainSenseLfp(k).LfpData)
        LFPL = horzcat(LFPL,[val.BrainSenseLfp(k).LfpData(j).Left.LFP]);
        LFPR = horzcat(LFPR,[val.BrainSenseLfp(k).LfpData(j).Right.LFP]);
        StimL = horzcat(StimL,[val.BrainSenseLfp(k).LfpData(j).Left.mA]);
        StimR = horzcat(StimR,[val.BrainSenseLfp(k).LfpData(j).Right.mA]);
    end
    
    
    %         Check for rollover
    %         if not(isempty(TDLFP))
    %             while TDLFP(end)>(TDx(1)-init_tick)/1000
    %                 TDx = TDx + 3276750;
    %             end
    %         end
    %         Check for rollover
    if k~=1
        StartTime(1) = datetime(val.BrainSenseLfp(1).FirstPacketDateTime,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z');
        StartTime(2) = datetime(val.BrainSenseLfp(k).FirstPacketDateTime,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z');
        diffStartTime = seconds(diff(StartTime));
        difftickStartTime = (TDx(1) - val.BrainSenseLfp(1).LfpData(1).TicksInMs)/1000;
        
        TDx = TDx+(diffStartTime*1000);
        % old code handle rollover. Added adjustments to use datetime for
        % sync each stream
%         while diffStartTime-15 > difftickStartTime % When labeled datetime values are larger than expectation from tick add rollover;15 is to allow some inaccuracy in non-rollover events. inaccuracies should be investigated
%             TDx = TDx + 3276750;
%             difftickStartTime = (TDx(1) - val.BrainSenseLfp(1).LfpData(1).TicksInMs)/1000;
%         end

    end
    TDLFP = [TDLFP,(TDx)/1000];
    
    LDoutL = [LDoutL,LFPL];
    LDoutR = [LDoutR,LFPR];
    StimoutL = [StimoutL,StimL];
    StimoutR = [StimoutR,StimR];
    clear TDx LFPL LFPR StimL StimR
    
end


%% Cut non-continuous sections

for i = length(TDLFP):-1:2
    if TDLFP(i)-TDLFP(i-1)>0.6
        TDLFP(i:end+1) = [nan,TDLFP(i:end)];
        LDoutL(i:end+1) = horzcat(nan,LDoutL(i:end));
        LDoutR(i:end+1) = horzcat(nan,LDoutR(i:end));
        StimoutL(i:end+1) = horzcat(nan,StimoutL(i:end));
        StimoutR(i:end+1) = horzcat(nan,StimoutR(i:end));
    end
end

%% Plotting
figure(301)
subplot(3,1,1)
xlimacrossL = gca;
figure(302)
subplot(3,1,1)
xlimacrossR = gca;

color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
stimcolor = {[0 0.4470 0.7410]*3/4,[0.8500 0.3250 0.0980]*3/4};

figure(301)
if not(isempty([LDoutL]))
    figure(301)
    subplot(3,1,3)
    yyaxis left
    hold on
    plot(TDLFP, LDoutL','Linestyle',':','Marker','.',...
        'MarkerSize',3,'color',color{1})
    ylabel('LFP (LSB)')
    try
    ylim([0, 5*mean(LDoutL,'omitnan')])
    catch
        ylim([-1,1])
    end
    
    yyaxis right
    hold on
    plot(TDLFP, StimoutL','Linestyle',':','Marker','.',...
        'MarkerSize',3,'color',stimcolor{2})
    ylabel('Stim (mA)')
    xlabel('Time (s)')
    ylim([-1,6])
    if xlimacrossL.XLim ~= [0 1]
        xlim(xlimacrossL.XLim)
    end
end

% figure(1)
% subplot(6,3,[16,17])
% if not(isempty([LDoutL]))
%     yyaxis left
%     hold on
%     plot(TDLFP, LDoutL','Linestyle',':','Marker','.',...
%         'MarkerSize',3,'color',color{1})
%     ylabel('LFP (LSB)')
%     try
%     ylim([0, max(LDoutL)])
%     catch
%         ylim([-1,1])
%     end
%     
%     yyaxis right
%     hold on
%     plot(TDLFP, StimoutL','Linestyle',':','Marker','.',...
%         'MarkerSize',3,'color',stimcolor{2})
%     ylabel('Stim (mA)')
%     xlabel('Time (s)')
%     ylim([-1,6])
%     if xlimacrossL.XLim ~= [0 1]
%         xlim(xlimacrossL.XLim)
%     end
% end


if not(isempty([LDoutR]))
    figure(302)
    subplot(3,1,3)
    yyaxis left
    hold on
    plot(TDLFP, LDoutR','Linestyle',':','Marker','.',...
        'MarkerSize',3,'color',color{1})
    ylabel('LFP (LSB)')
    try
    ylim([0, 5*mean(LDoutR,'omitnan')])
    catch
        ylim([0,1])
    end
    
    yyaxis right
    hold on
    plot(TDLFP, StimoutR','Linestyle',':','Marker','.',...
        'MarkerSize',3,'color',stimcolor{2})
    ylabel('Stim (mA)')
    xlabel('Time (s)')
    ylim([-1,6])
    if xlimacrossR.XLim ~= [0 1]
        xlim(xlimacrossR.XLim)
    end
end

% figure(2)
% subplot(6,3,[16,17])
% if not(isempty([LDoutR]))
%     yyaxis left
%     hold on
%     plot(TDLFP, LDoutR','Linestyle',':','Marker','.',...
%         'MarkerSize',3,'color',color{1})
%     ylabel('LFP (LSB)')
%     try
%         ylim([0, max(LDoutR)])
%     catch
%         ylim([0,1])
%     end
%     
%     yyaxis right
%     hold on
%     plot(TDLFP, StimoutR','Linestyle',':','Marker','.',...
%         'MarkerSize',3,'color',stimcolor{2})
%     ylabel('Stim (mA)')
%     xlabel('Time (s)')
%     ylim([-1,6])
%     if xlimacrossR.XLim ~= [0 1]
%         xlim(xlimacrossR.XLim)
%     end
% end


%% linkaxes
fig_axes_L=zeros(1,3);
fig_axes_R=zeros(1,3);

figure(301);
for sub_ind=1:3
    subplot(3,1,sub_ind);
    fig_axes_L(sub_ind)=gca;
end

figure(302);
for sub_ind=1:3
    subplot(3,1,sub_ind);
    fig_axes_R(sub_ind)=gca;
end

try
linkaxes(fig_axes_L, 'x')
catch
end
try
linkaxes(fig_axes_R, 'x')
catch
end
end