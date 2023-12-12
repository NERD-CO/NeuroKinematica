%% Combine Percept LFP with DLC Video

% Determine frames for video based on recording start

% mainDIR = 'D:\Dropbox\erin_clinic_video\Clin_2023-09-12_session1_soffoffLi_Lcam-emr-2023-09-16\videos';
mainDIR = 'C:\Users\erinr\Downloads\erin_clinic_video\erin_clinic_video\Clin_2023-09-12_session1_soffoffLi_Lcam-emr-2023-09-16\videos';

%% Tablet video

tab_vidLoc = [mainDIR, filesep , 'TabletVideo'];
cd(tab_vidLoc)

tab_vidObj = VideoReader('20230912_idea08_session001_rightCam-0000.mp4');
tab_vid = struct('cdata',zeros(tab_vidObj.Height,tab_vidObj.Width,3,'uint8'),'colormap',[]);

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

tabTab_vid = struct2table(tab_vid);

for fi = 2900:height(tabTab_vid)
    imshow(tabTab_vid.cdata{fi})
    title(num2str(fi))
    pause
end

% Start 196
% Stop 3059

%% Load in experimental camera
dlc_vidLoc = [mainDIR, filesep , 'DLCVideo'];
cd(dlc_vidLoc)

dlc_vidObj = VideoReader('20230912_idea08_session001_leftCam-0000.mp4');
dlc_vid = struct('cdata',zeros(dlc_vidObj.Height,dlc_vidObj.Width,3,'uint8'),'colormap',[]);

frami = 1;
while hasFrame(dlc_vidObj)
   dlc_vid(frami).cdata = readFrame(dlc_vidObj);
%    imshow(frame)
   frami = frami+1;
end
dlcTab_vid = struct2table(dlc_vid);
disp('Video1 done!')

%% Trim DLC video by start and stop from Tablet video
dlcTab_vid2use = dlcTab_vid(206:3010,:);

%% Load and Process CSV file
csv_vidLoc  = [mainDIR , filesep , 'CSVdata'];
save_matLoc = [mainDIR , filesep , 'MATdata'];
cd(csv_vidLoc)

dlcIO_processCSV('dirLOC',1,'userLOC',csv_vidLoc,'saveLOC',save_matLoc)

%% Load and Trim MAT file
cd(save_matLoc)
load('dlcDAT_20230912_idea08_R1.mat','outDATA')
d1_labDAT = outDATA.labelTab.leftCam;


%%

dlc_lab2use = d1_labDAT(206:3010,:);

% offsetSamples = 195;
% offsetSecs = 195/60;
totalNumSamplesV = height(dlc_lab2use);
totalNumSecs = totalNumSamplesV/60; % 60 fps

totalNumSampsLFP = floor(totalNumSecs*250);


%% Load and find first streaming from LFP
lfpLoc = 'D:\Dropbox\erin_clinic_video';
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


%%
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
% plot(dlc_lab2use2int.Elbow_x)
% hold on
% yyaxis left
elbowXsm = smoothdata(dlc_lab2use2int.Elbow_x,'gaussian',70);
nexttile
plot(xTime,elbowXsm);
xlim([0 round(max(ts_LFP))])
ylabel('X deflection')
xlabel('Time in seconds')
ecg = perceive_ecg(transpose(streamOfInt),250,0);
% plot(ecg.cleandata);
% yyaxis right
nexttile
plot(xTime,ecg.cleandata);
xlim([0 round(max(ts_LFP))])
ylabel('uV')
xlabel('Time in seconds')



%% 
close all
for di = 1:height(dlcTab_vid)
     imshow(dlcTab_vid.cdata{di})
     pause

end


%%

for ffi = 1:height(dlcTab_vid2use.cdata)
    imshow(dlcTab_vid2use.cdata{ffi})
    hold on
    plot(dlc_lab2use.PalmBase_x(ffi),dlc_lab2use.PalmBase_y(ffi),'ro')
    pause(0.01)
    cla
end

%% To do

% 1. clean up interpolated point of interest

% 2. clean up LFP with artifact removal

% 3. instanteous theta , beta and gamma


%% Check DLC labeled video
dlcLab_vidLoc = mainDIR;
cd(dlcLab_vidLoc)

dlc_lab_vidObj = VideoReader('20230912_idea08_session001_leftCam-0000DLC_resnet50_Clin_2023-09-12_session1_labeled.mp4');
dlc_lab_vid = struct('cdata',zeros(dlc_lab_vidObj.Height,dlc_lab_vidObj.Width,3,'uint8'),'colormap',[]);

frami = 1;
while hasFrame(dlc_lab_vidObj)
   dlc_lab_vid(frami).cdata = readFrame(dlc_lab_vidObj);
%    imshow(frame)
   frami = frami+1;
end
dlcLabTab_vid = struct2table(dlc_lab_vid);
disp('Video1 done!')


dlc_lablab2use = dlcLabTab_vid(206:3010,:);


%% Animated plot

plotXaxis = xTime;
kinematicY = elbowXsm;
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



