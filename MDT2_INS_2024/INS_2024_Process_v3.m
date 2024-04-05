% function [] = generateINS_DatasetS(subject2u , bodySIDE , LFPSIDE , group2use333)
% INPUT ARGUMENT for SUBJECT
% LOAD CSV file with:
% JSON left and right files


subjectID = 'MDT9';

pcName = getenv("COMPUTERNAME");

switch pcName
    case 'DESKTOP-I5CPDO7'
        mainDir = 'D:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data';
    case 'DESKTOP-EGFQKAI'
        mainDir = 'C:\Users\johna\Dropbox\ConferencePresentationsMeta\INS_2024\MDT2';
    case 'DESKTOP-FAGRV5G'
        sumMetaTLoc = 'E:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data';
        mainDir = 'E:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data';

    otherwise


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%        MDT9
% first test '11:56:27 AM MDT'
% second test '12:11:18 PM MDT'
%
% FIND GROUP ID FROM ----- 
% Therapy SnapShot for GROUP ID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mainFOLDloc = [mainDir , filesep , subjectID];

% --------- LOAD BOTH LEFT and RIGHT STN lfp
% LEFT STN LFP
STNLoc.Left = [mainFOLDloc,'\Streaming\Left STN'];
STNfile.Left = 'Report_Json_Session_Report_20230908T131441.json';
STN_fulFile.Left = [STNLoc.Left , filesep , STNfile.Left];
% RIGHT STN LFP
STNLoc.Right = [mainDir,filesep,subjectID,'\Streaming\Right STN'];
STNfile.Right = 'Report_Json_Session_Report_20230911T133035.json';
STN_fulFile.Right = [STNLoc.Right , filesep , STNfile.Right];

% Left STN / Right Body videos
cd(mainFOLDloc)

% SWITCH CASE FOR bodySIDE [input argument]

BodyVidDir.Right = [mainFOLDloc , '\RightBody'];
BodyCSVsrch.Right = fullfile(BodyVidDir.Right, '*.csv');
BodyCSVList.Right = dir(BodyCSVsrch.Right);
BodyCSVListCA.Right = {BodyCSVList.Right.name};

% leftSTN_rightBD_csvTable = getCSVinfo(BodyVidDir.Right , BodyCSVListCA.Right);
BodyVideoIndex.Right = fullfile(BodyVidDir.Right ,"VideoINFO.xlsx");

% Left and Right STN JSON Files
stnLEFTjson  = jsondecode(fileread(STN_fulFile.Left));
stnRIGHTjson = jsondecode(fileread(STN_fulFile.Right));

% Extract BrainSenseTimeDomain
stnLEFTstream = stnLEFTjson.BrainSenseTimeDomain;
stnRIGHTstream = stnRIGHTjson.BrainSenseTimeDomain;
stnLeftSt_tab = struct2table(stnLEFTstream);
stnRightSt_tab = struct2table(stnRIGHTstream);

% Left JSON is Left hemi primary (contralateral always collected)
% Select LEFT STN fields from LEFT JSON file
stnLeftSt_tabU.Left = stnLeftSt_tab(contains(stnLeftSt_tab.Channel,'LEFT'),:);
stnLeftSt_tabU.Right = stnLeftSt_tab(contains(stnLeftSt_tab.Channel,'RIGHT'),:);

stnRightSt_tabU.Left = stnRightSt_tab(contains(stnRightSt_tab.Channel,'LEFT'),:);
stnRightSt_tabU.Right = stnRightSt_tab(contains(stnRightSt_tab.Channel,'RIGHT'),:);

% process with a function
stnLeftSt_times.Left = getDAYtimes(stnLeftSt_tabU.Left.FirstPacketDateTime ,...
    stnLeftSt_tabU.Left.TimeDomainData);

% LEFT brain LFP + stimulation
stnBSLFP.Left = struct2table(stnLEFTjson.BrainSenseLfp);

% Get times from BS LFP
% BSlfptimes = getBSLFPtimes(stnBSLFP.Left.FirstPacketDateTime);




%% Get Group A (first group)

% EXTRACT ALL Groups from Left [primary] JSON 
groupIDchk = stnBSLFP.Left.TherapySnapshot;
groupIDs = cell(height(groupIDchk),1);
for gi = 1:height(groupIDchk)
    tmpGrpR = groupIDchk{gi}.ActiveGroup;
    tmpGrpR1 = extractAfter(tmpGrpR,'.');
    groupIDs{gi} = tmpGrpR1;
end

% ADD as input argument
group2use333 = 'b';
% stnLEFT_Stream.Left
[GROUP_Video_Names , GROUP_Video_MoveIND_Names, ~,...
    stnLEFT_Stream.Left]...
    = getGroupInfo(groupIDs, stnBSLFP.Left , stnLeftSt_times.Left ,...
    stnLeftSt_tabU.Left , BodyVideoIndex.Right , BodyVidDir.Right);

%%% SWITCH Case for Body Side

% DLC CVS convert video frame CSV
% Load and Process CSV file
csv_vidLoc  = BodyVidDir.Right;
save_matLoc = BodyVidDir.Right;
cd(csv_vidLoc)

% DO NOT NEED TO RE-RUN ------------------
% dlcIO_processCSV('dirLOC',1,'userLOC',csv_vidLoc,'saveLOC',save_matLoc,...
%     'selCSV',GROUP_Video_CSV,'USERid','JAT','hemiS','L','bodyS','R')

cd(save_matLoc)

% DLC name
matLIST = dir('*.mat');
matLISTu = {matLIST.name};

stnSIDE = 'L';

% FIX per GROUP
% csvNparts = split(GROUP_Video_Names.B , {'_','-'});
% csvMatSearch = [csvNparts{1},'_',csvNparts{3},'_',csvNparts{4},'_',stnSIDE,'STN'];
% dlcMATname = matLISTu{contains(matLISTu,csvMatSearch)};

[dlcMATnames] = getMATDLCnames(stnSIDE , matLISTu , GROUP_Video_Names);

load(dlcMATnames.A)
camFields.A = fieldnames(outDATA.labelTab);
d1_labDAT.A = outDATA.labelTab.(camFields.A{1});

load(dlcMATnames.B)
camFields.B = fieldnames(outDATA.labelTab);
d1_labDAT.B = outDATA.labelTab.(camFields.B{1});

% VIDEO WILL BE LONGER THAN LFP ---- TRIM FRAMES by LFP --- KEEP TRACK of
% INDEX
% STEP 1 ---- TRIM FRAMES by LFP
%%%%%%%%%%%%%% HARD CODED -------------------------------------------------
% ADD TO CSV

tabletFRAMEstart.A = 180;
tabletFRAMEend.A = 6445;
tabletFRAMEstart.B = 63;
tabletFRAMEend.B = 4624;

dlc_lab2use.A = d1_labDAT.A(tabletFRAMEstart.A:tabletFRAMEend.A,:);
dlc_lab2use.B = d1_labDAT.B(tabletFRAMEstart.B:tabletFRAMEend.B,:);

totalNumSamplesV.A = height(dlc_lab2use.A);
totalNumSecs.A = totalNumSamplesV.A/60; % 60 fps
totalNumSamplesV.B = height(dlc_lab2use.B);
totalNumSecs.B = totalNumSamplesV.B/60; % 60 fps

% totalNumSampsLFP = floor(totalNumSecs*250);

%%%%%%%%%%%%%%%%%%%% ------ LFP ---------------------%%%%%%%%%%%%%%%%%%%%%%
streamOfInt.Left.Left.GroupA = stnLEFT_Stream.Left.GroupA.TimeDomainData{1};
streamOfInt.Left.Left.GroupB = stnLEFT_Stream.Left.GroupB.TimeDomainData{1};

% Original signal at 60 Hz
ts_DLC.A = 0:1/60:(height(dlc_lab2use.A)-1)/60;
ts_DLC.B = 0:1/60:(height(dlc_lab2use.B)-1)/60;

% Target sampling rate at 250 Hz
ts_LFP.A = 0:1/250:(height(streamOfInt.Left.Left.GroupA)-1)/250;
ts_LFP.B = 0:1/250:(height(streamOfInt.Left.Left.GroupB)-1)/250;


allColNs = dlc_lab2use.B.Properties.VariableNames;
dlc_lab2use2int.A = table;
dlc_lab2use2int.B = table;
for coli = 1:width(dlc_lab2use.B)
    % Tmp col
    tmpCol = allColNs{coli};

    % Interpolation
    x_250A = interp1(ts_DLC.A, transpose(dlc_lab2use.A.(tmpCol)),...
        ts_LFP.A, 'spline');

    x_250B = interp1(ts_DLC.B, transpose(dlc_lab2use.B.(tmpCol)),...
        ts_LFP.B, 'spline');

    dlc_lab2use2int.A.(tmpCol) = transpose(x_250A);
    dlc_lab2use2int.B.(tmpCol) = transpose(x_250B);
end

trimFrame_int.A = linspace(tabletFRAMEstart.A,tabletFRAMEend.A, length(ts_LFP.A));
trimFrame_int.B = linspace(tabletFRAMEstart.B,tabletFRAMEend.B, length(ts_LFP.B));

%%
% GROUP Move indices ---- 
moveIND = readtable(GROUP_Video_MoveIND_Names.B);
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

    [~ , startINDi] = min(abs(startIND - trimFrame_int.B));
    [~ , endINDi] = min(abs(endIND - trimFrame_int.B));

    movKIN = [dlc_lab2use2int.B.fTip1_x(startINDi:endINDi) ,...
        dlc_lab2use2int.B.fTip1_y(startINDi:endINDi)];
    movLFP = streamOfInt.Left.Left.GroupB(startINDi:endINDi);

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

test = 1;







% contact list
% BSLFP - LFPDATA = 2 Hz (0.5 s) steps - FFT (POWER) / mA
% Groups.Final(1).ProgramSettings.SensingChannel(1).ElectrodeState{1,1}.
% :::: Electrode , mA , fraction????? 
% EACH cell of ElectrodeState is a different contact - LOOK for everything
% BUT CASE 
% Sensing channel is rows for each hemisphere
% 



%% FUNCTIONS with MAIN FUNCTION



function [GROUP_Video_Name , GROUP_Video_MoveIND_Name ,...
    GROUP_Video_CSV, stnL_LFP_St_KinLFP]...
    = getGroupInfo(allGroups, BS_LFP , LFPstreamTimes ,...
    LFPstreamTab , videoFname , videoDirectory)

% switch groupType
    % case 'a'
        % group2use = 'GROUP_A';
        % group2use2 = 'A';
    % case 'b'
        % group2use = 'GROUP_B';
        % group2use2 = 'B';
    % case 'c'
        % group2use = 'GROUP_C';
        % group2use2 = 'C';
% end

% Get Group A rows
groupA_inds = matches(allGroups,'GROUP_A');
stnL_LFP_BS_grA = BS_LFP(groupA_inds,:); % BS LFP 
% stnL_LFP_BStms_grA = BSlfptimes(groupA_inds,:);
stnL_Ttims_grA = LFPstreamTimes(groupA_inds,:);
stnLSt_tab_LH_grA = LFPstreamTab(groupA_inds,:); % Stream LFP

groupB_inds = matches(allGroups,'GROUP_B');
stnL_LFP_BS_grB = BS_LFP(groupB_inds,:);
stnL_Ttims_grB = LFPstreamTimes(groupB_inds,:);
stnLSt_tab_LH_grB = LFPstreamTab(groupB_inds,:);

% Pull in video durations

% Load Video info XslX
stnL_BD_vinfo = readtable(videoFname);
% Pull out Group A and Group B
groupA_Vinfo = stnL_BD_vinfo(contains(stnL_BD_vinfo.ftype,'A'),:);
groupB_Vinfo = stnL_BD_vinfo(contains(stnL_BD_vinfo.ftype,'B'),:);

[~ , stnL_BD_BS_row_grpA] = min(abs(groupA_Vinfo.duration - stnL_Ttims_grA.Duration));
[~ , stnL_BD_BS_row_grpB] = min(abs(groupB_Vinfo.duration - stnL_Ttims_grB.Duration));

% GROUP A - Condition TEST
stnL_BD_BS_testInfo.GroupA = stnL_Ttims_grA(stnL_BD_BS_row_grpA,:);
stnL_LFP_BS_KinLFP.GroupA = stnL_LFP_BS_grA(stnL_BD_BS_row_grpA,:);
stnL_LFP_St_KinLFP.GroupA = stnLSt_tab_LH_grA(stnL_BD_BS_row_grpA,:);

% GROUP A - Condition TEST
stnL_BD_BS_testInfo.GroupB = stnL_Ttims_grB(stnL_BD_BS_row_grpB,:);
stnL_LFP_BS_KinLFP.GroupB = stnL_LFP_BS_grB(stnL_BD_BS_row_grpB,:);
stnL_LFP_St_KinLFP.GroupB = stnLSt_tab_LH_grB(stnL_BD_BS_row_grpB,:);

GROUP_Video_CSV.A = groupA_Vinfo.csvName;
GROUP_Video_CSV.B = groupB_Vinfo.csvName;

cd(videoDirectory)

allVideoNames = dir('*.mp4');
allVideoNames2 = {allVideoNames.name};
allVideoNames3 = cellfun(@(x) [extractBefore(x,'_labeled'),'.mp4'],...
    allVideoNames2, 'UniformOutput',false);

GROUP_Video_Name.A = allVideoNames2(contains(allVideoNames3,...
    extractBefore(GROUP_Video_CSV.A{1},'.')));
GROUP_Video_Name.B = allVideoNames2(contains(allVideoNames3,...
    extractBefore(GROUP_Video_CSV.B{1},'.')));

GROUP_Video_MoveIND_Name.A = [extractBefore(GROUP_Video_Name.A{1},'.'),'_MoveIndex.csv'];
GROUP_Video_MoveIND_Name.B = [extractBefore(GROUP_Video_Name.B{1},'.'),'_MoveIndex.csv'];

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




function [outDLCmatNames] = getMATDLCnames(stnSIDE , matLISTu ,GROUP_Video_Names)

csvNpartsA = split(GROUP_Video_Names.A , {'_','-'});
csvMatSearchA = [csvNpartsA{1},'_',csvNpartsA{3},'_',csvNpartsA{4},'_',stnSIDE,'STN'];
outDLCmatNames.A = matLISTu{contains(matLISTu,csvMatSearchA)};

csvNpartsB = split(GROUP_Video_Names.B , {'_','-'});
csvMatSearchB = [csvNpartsB{1},'_',csvNpartsB{3},'_',csvNpartsB{4},'_',stnSIDE,'STN'];
outDLCmatNames.B = matLISTu{contains(matLISTu,csvMatSearchB)};

end



