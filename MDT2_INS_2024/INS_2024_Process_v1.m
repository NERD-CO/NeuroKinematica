pcName = getenv("COMPUTERNAME");

switch pcName
    case 'DESKTOP-I5CPDO7'
        mainDir = 'D:\Dropbox\ConferencePresentationsMeta\INS_2024\MDT2';
    case 'DESKTOP-EGFQKAI'
        mainDir = 'C:\Users\johna\Dropbox\ConferencePresentationsMeta\INS_2024\MDT2';

    otherwise


end

% loop through and store as structs
leftSTNLoc = [mainDir,'\Streaming\Left STN'];
leftSTNfile = 'Report_Json_Session_Report_20230908T131441.json';

rightSTNLoc = [mainDir,'\Streaming\Right STN'];
rightSTNfile = 'Report_Json_Session_Report_20230911T133035.json';

% cd(leftSTNLoc)
leftFF = [leftSTNLoc , filesep , leftSTNfile];
rightFF = [rightSTNLoc , filesep , rightSTNfile];


% Left STN / Right Body videos
leftSTN_rightBD_dir = [mainDir , filesep , 'RightBody'];
search_path = fullfile(leftSTN_rightBD_dir, '*.csv');
leftSTN_rightBD_csvLt = dir(search_path );
leftSTN_rightBD_csvL = {leftSTN_rightBD_csvLt.name};

leftSTN_rightBD_csvTable = getCSVinfo(leftSTN_rightBD_dir , leftSTN_rightBD_csvL);
leftSTN_rightBD_Vinfo_fname = fullfile(leftSTN_rightBD_dir ,"VideoINFO.xlsx");

% Right STN / Left Body videos

stnLEFTjson  = jsondecode(fileread(leftFF));
stnRIGHTjson = jsondecode(fileread(rightFF));

%
stnLEFTstream = stnLEFTjson.BrainSenseTimeDomain;
stnLSt_tab = struct2table(stnLEFTstream);
% Get relevant rows for hemi recorded
stnLSt_tab_LH = stnLSt_tab(contains(stnLSt_tab.Channel,'LEFT'),:);

% process with a function
stnLEFTsTtimes = getDAYtimes(stnLSt_tab_LH.FirstPacketDateTime ,...
    stnLSt_tab_LH.TimeDomainData);

% LEFT brain LFP + stimulation
stnLEFTBSlfp = struct2table(stnLEFTjson.BrainSenseLfp);

% Get times from BS LFP
[BSlfptimes] = getBSLFPtimes(stnLEFTBSlfp.FirstPacketDateTime);
%

% first test '11:56:27 AM MDT'
% second test '12:11:18 PM MDT'
%
% FIND GROUP ID FROM ----- 
% Therapy SnapShot for GROUP ID

%% Get Group A (first group)
groupIDchk = stnLEFTBSlfp.TherapySnapshot;
groupIDs = cell(height(groupIDchk),1);
for gi = 1:height(groupIDchk)

    tmpGrpR = groupIDchk{gi}.ActiveGroup;
    tmpGrpR1 = extractAfter(tmpGrpR,'.');
    groupIDs{gi} = tmpGrpR1;

end

group2use333 = 'b';

[GROUP_Video_Name , GROUP_Video_MoveIND_Name, GROUP_Video_CSV,...
    stnL_LFP_St_GrpA_KinLFP]...
    = getGroupInfo(group2use333 , groupIDs, stnLEFTBSlfp , stnLEFTsTtimes ,...
    stnLSt_tab_LH , leftSTN_rightBD_Vinfo_fname);




% DLC CVS convert video frame CSV
% Load and Process CSV file
csv_vidLoc  = leftSTN_rightBD_dir;
save_matLoc = leftSTN_rightBD_dir;
cd(csv_vidLoc)

dlcIO_processCSV('dirLOC',1,'userLOC',csv_vidLoc,'saveLOC',save_matLoc,...
    'selCSV',GROUP_Video_CSV,'USERid','JAT','hemiS','L','bodyS','R')

cd(save_matLoc)

% DLC name
matLIST = dir('*.mat');
matLISTu = {matLIST.name};

stnSIDE = 'L';

csvNparts = split(GROUP_Video_Name , {'_','-'});
csvMatSearch = [csvNparts{1},'_',csvNparts{3},'_',csvNparts{4},'_',stnSIDE,'STN'];
dlcMATname = matLISTu{contains(matLISTu,csvMatSearch)};

load(dlcMATname)
camFields = fieldnames(outDATA.labelTab);
d1_labDAT = outDATA.labelTab.(camFields{1});

% VIDEO WILL BE LONGER THAN LFP ---- TRIM FRAMES by LFP --- KEEP TRACK of
% INDEX
% STEP 1 ---- TRIM FRAMES by LFP
%%%%%%%%%%%%%% HARD CODED

switch group2use333
    case 'a'
        tabletFRAMEstart = 180;
        tabletFRAMEend = 6445;
    case 'b'
        tabletFRAMEstart = 63;
        tabletFRAMEend = 4624;
end

dlc_lab2use = d1_labDAT(tabletFRAMEstart:tabletFRAMEend,:);

totalNumSamplesV = height(dlc_lab2use);
totalNumSecs = totalNumSamplesV/60; % 60 fps

totalNumSampsLFP = floor(totalNumSecs*250);

% LFP
streamOfInt = stnL_LFP_St_GrpA_KinLFP.TimeDomainData{1};
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

trimFrB_int = linspace(tabletFRAMEstart,tabletFRAMEend, length(ts_LFP));


%%
% Pull in 
moveIND = readtable(GROUP_Video_MoveIND_Name);
moveINDc = moveIND(~moveIND.BeginF == 0,:);
if matches(group2use333,'a')
    moveINDc = moveINDc(1:8,:); % remove last too short in time
end

groupLFPAmp = zeros(height(moveINDc),1); 
groupKinAmp = zeros(height(moveINDc),1);
for mii = 1:height(moveINDc)

    tmpMOve = moveINDc(mii,:);
    startIND = tmpMOve.BeginF;
    endIND = tmpMOve.EndF;

    [~ , startINDi] = min(abs(startIND - trimFrB_int));
    [~ , endINDi] = min(abs(endIND - trimFrB_int));

    movKIN = [dlc_lab2use2int.fTip1_x(startINDi:endINDi) ,...
        dlc_lab2use2int.fTip1_y(startINDi:endINDi)];
    movLFP = streamOfInt(startINDi:endINDi);

    ts_LFPtmp = 0:1/250:(height(movLFP)-1)/250;

    tiledlayout(2,1,"TileSpacing","tight")
    xTime = ts_LFPtmp;
    pointXsm = smoothdata(movKIN(:,1),'gaussian',5);
    pointYsm = smoothdata(movKIN(:,2),'gaussian',5);
    nexttile
    plot(xTime,pointXsm,'k');
    % hold on
    % plot(xTime,pointYsm,'r');
    xlim([0 round(max(ts_LFPtmp))])
    ylabel('X deflection')
    xlabel('Time in seconds')

    Fs = 250;
    T = 1/250;

    ecgClean = perceive_ecg(transpose(movLFP),250,0);
    nexttile

    cleanLFP = ecgClean.cleandata;

    % Calculate instantaneous phase using Hilbert transform
    hilbert_eeg = hilbert(cleanLFP);
    inst_phase = angle(hilbert_eeg);
    inst_freq = diff(unwrap(inst_phase))/(2*pi*T);

    % Filter instantaneous frequency for 13-30 Hz
    bandpass_filt = designfilt('bandpassiir','FilterOrder',4, ...
        'HalfPowerFrequency1',13,'HalfPowerFrequency2',30, ...
        'SampleRate',Fs);
    inst_freq_filtered = filtfilt(bandpass_filt, inst_freq);
    inst_freq_filtered2 = [inst_freq_filtered , 0];

    plot(xTime,inst_freq_filtered2);
    xlim([0 round(max(ts_LFPtmp))])
    ylabel('voltage')
    xlabel('Time in seconds')

    groupLFPAmp(mii) = max(abs(inst_freq_filtered2));
    groupKinAmp(mii) = max(abs(pointXsm));

    % pause 
    close all
    clc

 

end


%%%% SOME CURSORY ANALYSIS between BETA and MOVE KIN


%%%% SOME CURSORY ANALYSIS between GROUP A and GROUP B



%%
% close all
% 
% tiledlayout(2,1,"TileSpacing","tight")
% xTime = ts_LFP;
% % plot(dlc_lab2use2int.Elbow_x)
% % hold on
% % yyaxis left
% elbowXsm = smoothdata(dlc_lab2use2int.Elbow_x,'gaussian',70);
% nexttile
% plot(xTime,elbowXsm);
% xlim([0 round(max(ts_LFP))])
% ylabel('X deflection')
% xlabel('Time in seconds')
% ecg = perceive_ecg(transpose(streamOfInt),250,0);
% % plot(ecg.cleandata);
% % yyaxis right
% nexttile
% plot(xTime,ecg.cleandata);
% xlim([0 round(max(ts_LFP))])
% ylabel('uV')
% xlabel('Time in seconds')
% 









% contact list
% BSLFP - LFPDATA = 2 Hz (0.5 s) steps - FFT (POWER) / mA
% Groups.Final(1).ProgramSettings.SensingChannel(1).ElectrodeState{1,1}.
% :::: Electrode , mA , fraction????? 
% EACH cell of ElectrodeState is a different contact - LOOK for everything
% BUT CASE 
% Sensing channel is rows for each hemisphere
% 



%%



function [GROUPA_Video_Name ,GROUPA_Video_MoveIND_Name ,...
    GROUPA_Video_CSV, stnL_LFP_St_GrpA_KinLFP]...
    = getGroupInfo(groupType , allGroups, BS_LFP , LFPstreamTimes ,...
    LFPstreamTab , videoFname)

switch groupType
    case 'a'
        group2use = 'GROUP_A';
        group2use2 = 'A';


    case 'b'
        group2use = 'GROUP_B';
        group2use2 = 'B';

    case 'c'
        group2use = 'GROUP_C';
        group2use2 = 'C';
end

% Get Group A rows
groupA_inds = matches(allGroups,group2use);
stnL_LFP_BS_grA = BS_LFP(groupA_inds,:); % BS LFP 
% stnL_LFP_BStms_grA = BSlfptimes(groupA_inds,:);
stnL_Ttims_grA = LFPstreamTimes(groupA_inds,:);
stnLSt_tab_LH_grA = LFPstreamTab(groupA_inds,:); % Stream LFP

% Pull in video durations

% Load Video info XslX
stnL_BD_R_vinfo = readtable(videoFname);
% Pull out Group A and Group B
groupA_Vinfo = stnL_BD_R_vinfo(contains(stnL_BD_R_vinfo.ftype,group2use2),:);

[~ , stnL_BD_BS_row] = min(abs(groupA_Vinfo.duration - stnL_Ttims_grA.Duration));

% GROUP A - Condition TEST
stnL_BD_BS_GrpA_testInfo = stnL_Ttims_grA(stnL_BD_BS_row,:);
stnL_LFP_BS_GrpA_KinLFP = stnL_LFP_BS_grA(stnL_BD_BS_row,:);
stnL_LFP_St_GrpA_KinLFP = stnLSt_tab_LH_grA(stnL_BD_BS_row,:);

GROUPA_Video_CSV = groupA_Vinfo.csvName;
allVideoNames = dir('*.mp4');
allVideoNames2 = {allVideoNames.name};
allVideoNames3 = cellfun(@(x) [extractBefore(x,'_labeled'),'.mp4'],...
    allVideoNames2, 'UniformOutput',false);
GROUPA_Video_Name = allVideoNames2(contains(allVideoNames3,...
    extractBefore(GROUPA_Video_CSV{1},'.')));
GROUPA_Video_MoveIND_Name = [extractBefore(GROUPA_Video_Name{1},'.'),'_MoveIndex.csv'];









end



function [outTAB] = getDAYtimes(inputTIMES , inputSAMPLES)

dayTIMES = cell(height(inputTIMES),1);
dayS = cell(height(inputTIMES),1);
fullDtime = cell(height(inputTIMES),1);
durations = zeros(height(inputTIMES),1);
streamOffsets = nan(height(inputTIMES),1);
for ti = 1:height(inputTIMES)

    inputTi = inputTIMES{ti};

    % The input date-time string
    % dateTimeStr = '2023-09-08T17:47:31.000Z';

    % Convert the string to a datetime object
    dateTimeObj = datetime(inputTi, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z');

    % Set the time zone to UTC (Coordinated Universal Time) since the input is in UTC
    dateTimeObj.TimeZone = 'UTC';

    % Convert to Mountain Time
    dateTimeObj_Mountain = datetime(dateTimeObj, 'TimeZone', 'America/Denver');

    % Extract the time component in AM/PM format
    timeComponent_AMPM = datetime(dateTimeObj_Mountain,'Format','hh:mm:ss a z');
    timeComponent_DATE = datetime(dateTimeObj_Mountain,'Format','dd-MMM-yyyy');

    % Display the time component
    % disp(['Time in Mountain Time (AM/PM): ', timeComponent_AMPM]);
    dayTIMES{ti} = timeComponent_AMPM;
    dayS{ti} = timeComponent_DATE;
    fullDtime{ti} = dateTimeObj_Mountain;
    durations(ti) = round(length(inputSAMPLES{ti})/250);

end


for di = 1:height(fullDtime)

    if di < height(fullDtime)
        t1 = fullDtime{di};
        t2 = fullDtime{di + 1};

        streamOffsets(di) = seconds(time(between(t1,t2)));

    end

end

dayT2 = cellfun(@(x) char(x), dayTIMES, 'UniformOutput',false);
dayS2 = cellfun(@(x) char(x), dayS, 'UniformOutput',false);
fDt = cellfun(@(x) char(x), fullDtime, 'UniformOutput',false);

outTAB = table(dayT2 , dayS2  , fDt , durations , streamOffsets,...
    'VariableNames',{'TimeOccur','DayOccur','FullNAT','Duration','Offset'});

% end

end



function [outTAB] = getBSLFPtimes(inputTIMES)

dayTIMES = cell(height(inputTIMES),1);
dayS = cell(height(inputTIMES),1);
fullDtime = cell(height(inputTIMES),1);

for ti = 1:height(inputTIMES)

    inputTi = inputTIMES{ti};

    % The input date-time string
    % dateTimeStr = '2023-09-08T17:47:31.000Z';

    % Convert the string to a datetime object
    dateTimeObj = datetime(inputTi, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z');

    % Set the time zone to UTC (Coordinated Universal Time) since the input is in UTC
    dateTimeObj.TimeZone = 'UTC';

    % Convert to Mountain Time
    dateTimeObj_Mountain = datetime(dateTimeObj, 'TimeZone', 'America/Denver');

    % Extract the time component in AM/PM format
    timeComponent_AMPM = datetime(dateTimeObj_Mountain,'Format','hh:mm:ss a z');
    timeComponent_DATE = datetime(dateTimeObj_Mountain,'Format','dd-MMM-yyyy');

    % Display the time component
    % disp(['Time in Mountain Time (AM/PM): ', timeComponent_AMPM]);
    dayTIMES{ti} = timeComponent_AMPM;
    dayS{ti} = timeComponent_DATE;
    fullDtime{ti} = dateTimeObj_Mountain;

end

dayT2 = cellfun(@(x) char(x), dayTIMES, 'UniformOutput',false);
dayS2 = cellfun(@(x) char(x), dayS, 'UniformOutput',false);
fDt = cellfun(@(x) char(x), fullDtime, 'UniformOutput',false);

outTAB = table(dayT2 , dayS2  , fDt,...
    'VariableNames',{'TimeOccur','DayOccur','FullNAT'});

end









function [outTABLE] = getCSVinfo(dirLOC , csvLIST)

cd(dirLOC)

[rowS , colS] = size(csvLIST);
if colS > rowS
    csvLIST = transpose(csvLIST);
end

numFrames = zeros(length(csvLIST),1);
for ci = 1:length(csvLIST)
    tmpCSV = readtable(csvLIST{ci});
    numFrames(ci) = height(tmpCSV);
end

outTABLE = table(csvLIST , numFrames,'VariableNames',{'CSVname','NumFrames'});


end












