%% Combine Percept LFP with DLC Video

% Determine frames for video based on recording start

mainDIR = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\DLC_LFP';

casedate_hem = '09_12_2023_LSTN';
mainDIR2 = [mainDIR, filesep , casedate_hem];
cd(mainDIR2)

%% Converted left and right cam case videos

vidLoc = [mainDIR2, filesep , 'ConvertedVideos'];
cd(vidLoc)

%% Tablet video

% Assuming 9_12_2023 LSTN, R body Tablet Frames were captured via Right Cam
TabletVideos = {'20230912_idea08_session001_rightCam-0000DLC_resnet50_Clin_2023-09-12_LSTN_v3Oct3shuffle1_100000_labeled-converted.mp4',
                '20230912_idea08_session003_rightCam-0000DLC_resnet50_Clin_2023-09-12_LSTN_v3Oct3shuffle1_100000_labeled-converted.mp4',
                '20230912_idea08_session005_rightCam-0000DLC_resnet50_Clin_2023-09-12_LSTN_v3Oct3shuffle1_100000_labeled-converted.mp4',
                '20230912_idea08_session007_rightCam-0000DLC_resnet50_Clin_2023-09-12_LSTN_v3Oct3shuffle1_100000_labeled-converted.mp4',
                '20230912_idea08_session009_rightCam-0000DLC_resnet50_Clin_2023-09-12_LSTN_v3Oct3shuffle1_100000_labeled-converted.mp4',
                '20230912_idea08_session013_rightCam-0000DLC_resnet50_Clin_2023-09-12_LSTN_v3Oct3shuffle1_100000_labeled-converted.mp4',
                '20230912_idea08_session015_rightCam-0000DLC_resnet50_Clin_2023-09-12_LSTN_v3Oct3shuffle1_100000_labeled-converted.mp4',
                '20230912_idea08_session018_rightCam-0000DLC_resnet50_Clin_2023-09-12_LSTN_v3Oct3shuffle1_100000_labeled-converted.mp4',
                '20230912_idea08_session020_rightCam-0000DLC_resnet50_Clin_2023-09-12_LSTN_v3Oct3shuffle1_100000_labeled-converted.mp4',
                '20230912_idea08_session022_rightCam-0000DLC_resnet50_Clin_2023-09-12_LSTN_v3Oct3shuffle1_100000_labeled-converted.mp4'};

% for tab_vid_i = 1:length(TabletVideos)
%     tab_vidObj = VideoReader(TabletVideos{tab_vid_i});
%     tab_vid = struct('cdata',zeros(tab_vidObj.Height,tab_vidObj.Width,3,'uint8'),'colormap',[]);
% end


tab_vidObj = VideoReader('20230912_idea08_session001_rightCam-0000DLC_resnet50_Clin_2023-09-12_LSTN_v3Oct3shuffle1_100000_labeled-converted.mp4');
tab_vid = struct('cdata',zeros(tab_vidObj.Height,tab_vidObj.Width,3,'uint8'),'colormap',[]);

%% Load in experimental camera 2

% dlcVideos 

dlc_vidObj = VideoReader('20230912_idea08_session001_leftCam-0000-converted.mp4');
dlc_vid = struct('cdata',zeros(dlc_vidObj.Height,dlc_vidObj.Width,3,'uint8'),'colormap',[]);


%% Convert Tablet Video to dataframe

frami = 1;
while hasFrame(tab_vidObj)
   tab_vid(frami).cdata = readFrame(tab_vidObj);
%    imshow(frame)
   frami = frami+1;
end
whos tab_vid
v1who = whos('tab_vid');
round(v1who.bytes/1000000/1000,2)
disp('Video1 done!')

%% Determine which frames to trim on based presence of tablet
% Loop through / Plot / Title with Frame #

tab_vid_table = struct2table(tab_vid);

for fi = 1:height(tab_vid_table)
    imshow(tab_vid_table.cdata{fi})
    title(num2str(fi))
    pause
end

% session001 
% StartFrame =  196; 
% StopFrame = 3060;


%% Convert Experimental Video 2 to dataframe

frami = 1;
while hasFrame(dlc_vidObj)
   dlc_vid(frami).cdata = readFrame(dlc_vidObj);
%    imshow(frame)
   frami = frami+1;
end
dlcTab_vid = struct2table(dlc_vid);
disp('Video1 done!')

% session001 
StartFrame =  206; % 206
StopFrame = 3060;  % 210

%% Trim DLC video by start and stop from Tablet video

dlcTab_vid2use = dlcTab_vid(StartFrame:StopFrame,:);

%% Load and Process CSV file

post_DLCproc_Dir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\post_DLCproc\09_12_2023_LSTN_v2';
cd(post_DLCproc_Dir)

csv_vidLoc  = [post_DLCproc_Dir , filesep , 'csv folder'];
save_matLoc = [post_DLCproc_Dir , filesep , 'mat folder'];
cd(csv_vidLoc)

addpath 'C:\Users\erinr\OneDrive\Documents\GitHub\NeuroKinematica\DLC_Processing'

dlcIO_processCSV('dirLOC',1,'userLOC',csv_vidLoc,'saveLOC',save_matLoc)

%% Load and Trim MAT file

% cd(save_matLoc)
% load('dlcDAT_20230912_idea08_R1.mat','outDATA')
% % d1_labDAT = outDATA.labelTab.leftCam;
% d1_labDAT = outDATA.labelTab.rightCam;

% Load MAT file corresponding to session
load('dlcDAT_20230912_session001_idea08_resnet50_rightCam-0000DLC.mat','outDATA')
d1_labDAT2 = outDATA;

% Trim MAT file
dlc_lab2use = d1_labDAT2(StartFrame:StopFrame,:);

% offsetSamples = 195;
% offsetSecs = 195/60;
totalNumSamplesV = height(dlc_lab2use);
totalNumSecs = totalNumSamplesV/60; % 60 fps

totalNumSampsLFP = floor(totalNumSecs*250);


%% Load and find first streaming from LFP

% lfpLoc = 'C:\Users\erinr\Downloads\erin_clinic_video\erin_clinic_video';

lfpLoc = mainDIR2;
cd(lfpLoc)
lfpFname = 'Report_Json_Session_Report_20230912T115956.json';

js = jsondecode(fileread(lfpFname));
streaming = js.BrainSenseTimeDomain;

% Convert to table
streamTAB = struct2table(streaming);

% Trim by STN of interest
streamLEFT = streamTAB(contains(streamTAB.Channel,'LEFT'),:);


% Right body / Left STN
% Specific Row of interest
rowfromTab = 3;
streamOfInt = streamLEFT.TimeDomainData{rowfromTab};


%% Sampling Rates

% Original signal at 60 Hz
ts_DLC = 0:1/60:(height(dlc_lab2use)-1)/60;

% Target sampling rate at 250 Hz
ts_LFP = 0:1/250:(height(streamOfInt)-1)/250;

allColNs = dlc_lab2use.Properties.VariableNames;
dlc_lab2use2int = table;
for coli = 1:width(dlc_lab2use)

    % Tmp col
    tmpCol = allColNs{coli};
    
    % Interpolation
    x_250 = interp1(ts_DLC, transpose(dlc_lab2use.(tmpCol)),...
        ts_LFP, 'spline');

    dlc_lab2use2int.(tmpCol) = transpose(x_250);
end

%% Kinematic and LFP plot
close all

tiledlayout(2,1,"TileSpacing","tight")
xTime = ts_LFP;
% plot(dlc_lab2use2int.fTip1_x)
% hold on
% yyaxis left
fTip1_X_smooth = smoothdata(dlc_lab2use2int.fTip1_x,'gaussian',70);
nexttile
plot(xTime,fTip1_X_smooth);
xlim([0 round(max(ts_LFP))])
ylabel('fTip1 X deflection')
xlabel('Time in seconds')
ecg = perceive_ecg(transpose(streamOfInt),250,0);
% plot(ecg.cleandata);
% yyaxis right
nexttile
plot(xTime,ecg.cleandata);
xlim([0 round(max(ts_LFP))])
ylabel('LFP (uV)')
xlabel('Time in seconds')


% %% 
% close all
% for di = 1:height(dlcTab_vid)
%      imshow(dlcTab_vid.cdata{di})
%      pause
% 
% end


%% Check DLC labeled video
dlcLab_vidLoc = vidLoc;
cd(dlcLab_vidLoc)

dlc_lab_vidObj = VideoReader('20230912_idea08_session001_rightCam-0000DLC_resnet50_Clin_2023-09-12_LSTN_v3Oct3shuffle1_100000_labeled-converted.mp4');
dlc_lab_vid = struct('cdata',zeros(dlc_lab_vidObj.Height,dlc_lab_vidObj.Width,3,'uint8'),'colormap',[]);

frami = 1;
while hasFrame(dlc_lab_vidObj)
   dlc_lab_vid(frami).cdata = readFrame(dlc_lab_vidObj);
%    imshow(frame)
   frami = frami+1;
end
dlcLabTab_vid = struct2table(dlc_lab_vid);
disp('Video1 done!')


dlc_lablab2use = dlcLabTab_vid(StartFrame:StopFrame,:);


%% Animated plot

plotXaxis = xTime;
kinematicY = fTip1_X_smooth;
lfpY = ecg.cleandata;

% dlc_lablab2use


% 4.1 samples per frame
stePS = round(linspace(1,11688,2805));

for fi = 1:height(dlc_lablab2use)
    subplot(4,1,1:2)

    imshow(dlc_lablab2use.cdata{fi})

    subplot(4,1,3)
    plot(plotXaxis(1:stePS(fi)),kinematicY(1:stePS(fi)))
    ylim([min(kinematicY) max(kinematicY)])
    xlim([0 max(plotXaxis)])

    subplot(4,1,4)
    plot(plotXaxis(1:stePS(fi)),lfpY(1:stePS(fi)))
    ylim([min(lfpY) max(lfpY)])
    xlim([0 max(plotXaxis)])

    if fi == 1
        pause
    end

    pause(0.01)

end




%%
for ffi = 1:height(dlc_lablab2use.cdata)
    imshow(dlc_lablab2use.cdata{ffi})
    hold on
    % plot(dlc_lab2use.PalmBase_x(ffi),dlc_lab2use.PalmBase_y(ffi),'ro')
    pause(0.01)
    cla
end



