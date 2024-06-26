function [] = generateINS_DatasetS_v2(subjectID , stnHemi)
% INPUT ARGUMENT for SUBJECT
% LOAD CSV file with:
% JSON left and right files

% subjectID = 'MDT9';
% stnHemi = 'left';

pcName = getenv("COMPUTERNAME");

switch pcName
    case 'DESKTOP-I5CPDO7'
        mainDir = 'D:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data';
    case 'DESKTOP-EGFQKAI'
        mainDir = 'C:\Users\johna\Dropbox\ConferencePresentationsMeta\INS_2024\MDT2';
    case 'DESKTOP-FAGRV5G'
        sumMetaTLoc = 'E:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data';
        mainDir = 'E:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data';
        saveDIR = 'E:\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data\finalResults';

    otherwise


end
cd(sumMetaTLoc)
metaTAB = readtable('SummaryMETA.xlsx');
subjectTABle = metaTAB(matches(metaTAB.SUBJECT,subjectID),:);
stnhemiTABle = subjectTABle(matches(subjectTABle.STNh,stnHemi),:);

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

switch stnHemi
    case 'left'
        STNLoc.Left = [mainFOLDloc,'\Streaming'];
        STNfileNt = unique(stnhemiTABle.JSON(~matches(stnhemiTABle.JSON,'NaN')));
        STNfile.Left = STNfileNt{1};
        STN_fulFile.Left = [STNLoc.Left , filesep , STNfile.Left];
        body2use = 'right';
        stnLEFTjson  = load(STN_fulFile.Left);
        stnLeftSt_tab = stnLEFTjson.BSTD;
        % Left JSON is Left hemi primary (contralateral always collected)
        stnPrimeSt_tabU.Left = stnLeftSt_tab(contains(stnLeftSt_tab.Channel,'LEFT'),:);
        stnPrimeSt_tabU.Right = stnLeftSt_tab(contains(stnLeftSt_tab.Channel,'RIGHT'),:);

        stnPrimeSt_times.Left = getDAYtimes(stnPrimeSt_tabU.Left.FirstPacketDateTime ,...
            stnPrimeSt_tabU.Left.TimeDomainData);

        stnPrimeSt_times.Right = getDAYtimes(stnPrimeSt_tabU.Right.FirstPacketDateTime ,...
            stnPrimeSt_tabU.Right.TimeDomainData);

        stnBSLFP.Left = stnLEFTjson.BSLFP;
        groupIDchk = stnBSLFP.Left.TherapySnapshot;

        primarySTN = 'Left';
        secondarySTN = 'Right';

    case 'right'
        STNLoc.Right = [mainFOLDloc,'\Streaming'];
        STNfileNt = unique(stnhemiTABle.JSON(~matches(stnhemiTABle.JSON,'NaN')));
        STNfile.Right = STNfileNt{1};
        STN_fulFile.Right = [STNLoc.Right , filesep , STNfile.Right];
        body2use = 'left';
        stnRIGHTjson = load(STN_fulFile.Right);
        stnRightSt_tab = stnRIGHTjson.BSTD;
        % Right JSON is Right hemi primary (contralateral always collected)
        stnPrimeSt_tabU.Left = stnRightSt_tab(contains(stnRightSt_tab.Channel,'LEFT'),:);
        stnPrimeSt_tabU.Right = stnRightSt_tab(contains(stnRightSt_tab.Channel,'RIGHT'),:);

        stnPrimeSt_times.Left = getDAYtimes(stnPrimeSt_tabU.Left.FirstPacketDateTime ,...
            stnPrimeSt_tabU.Left.TimeDomainData);

        stnPrimeSt_times.Right = getDAYtimes(stnPrimeSt_tabU.Right.FirstPacketDateTime ,...
            stnPrimeSt_tabU.Right.TimeDomainData);

        stnBSLFP.Right = stnRIGHTjson.BSLFP;
        groupIDchk = stnBSLFP.Right.TherapySnapshot;

        primarySTN = 'Right';
        secondarySTN = 'Left';
end
% Left STN / Right Body videos
cd(mainFOLDloc)

% SWITCH CASE FOR bodySIDE [input argument]
switch body2use
    case 'right'
        BodyVidDir.Right = [mainFOLDloc , '\RightBody'];
        primaryBody = 'Right';
    case 'left'
        BodyVidDir.Left = [mainFOLDloc , '\LeftBody'];
        primaryBody = 'Left';
end



%% Get Group A (first group)

% EXTRACT ALL Groups from Left [primary] JSON 

groupIDs = cell(height(groupIDchk),1);
for gi = 1:height(groupIDchk)

    if iscell(groupIDchk)
        tmpGrpR = groupIDchk{gi}.ActiveGroup;
    else
        tmpGrpR = groupIDchk(gi).ActiveGroup;
    end

    tmpGrpR1 = extractAfter(tmpGrpR,'.');
    groupIDs{gi} = tmpGrpR1;
end

%%%% LOOP THROUGH GROUPS A AND B --- USING TABLE
outDATAFin = struct;
outDATAFin.SubID = subjectID;
outDATAFin.STNhemi = stnHemi;
groups2id = {'GroupA','GroupB'};
for groupIi = 1:2

    switch groupIi
        case 1 % GroupA
            groupRind = matches(stnhemiTABle.GROUP,'GroupA');
        case 2 % GroupB
            groupRind = matches(stnhemiTABle.GROUP,'GroupB');
    end

    outDATAFin.conditionID.(groups2id{groupIi}) = stnhemiTABle.CONDITION{groupRind};
    outDATAFin.videoID.(groups2id{groupIi}) = stnhemiTABle.VIDEO{groupRind};
    dlcID.(groups2id{groupIi}) = stnhemiTABle.DLC{groupRind};
    moveINDid.(groups2id{groupIi}) = stnhemiTABle.moveIND{groupRind};

    [ stnPRIME_Stream.(primarySTN).(groups2id{groupIi}).BDBS ,...
        stnPRIME_Stream.(primarySTN).(groups2id{groupIi}).BSLFP,...
        stnPRIME_Stream.(primarySTN).(groups2id{groupIi}).TS]...
        = getGroupInfo(groupIDs, groupIi, stnBSLFP.(primarySTN) ,...
        stnPrimeSt_times.(primarySTN) ,...
        stnPrimeSt_tabU.(primarySTN), BodyVidDir.(primaryBody), moveINDid.(groups2id{groupIi}));

    load(dlcID.(groups2id{groupIi}))

    camFields.(groups2id{groupIi}) = fieldnames(outDATA.labelTab);
    d1_labDAT.(groups2id{groupIi}) =...
        outDATA.labelTab.(camFields.(groups2id{groupIi}){1});

    % CLEAN UP with CONFIDENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    likeLIhNs1 = d1_labDAT.(groups2id{groupIi}).Properties.VariableNames;
    likeLIhNs2 = likeLIhNs1(contains(likeLIhNs1,{'likelihood'}));

    for lli = 1:length(likeLIhNs2)

        tmpLLI = d1_labDAT.(groups2id{groupIi}).(likeLIhNs2{lli});
        thresh = 0.75;
        fillwZ = tmpLLI < thresh;

        lliName = extractBefore(likeLIhNs2{lli},'_');

        xLOCation = [lliName,'_x'];
        yLOCation = [lliName,'_y'];
        d1_labDAT.(groups2id{groupIi}).(xLOCation)(fillwZ) = nan;
        d1_labDAT.(groups2id{groupIi}).(yLOCation)(fillwZ) = nan;


        plot(d1_labDAT.(groups2id{groupIi}).(yLOCation))

    end

    % CLEAN UP with CONFIDENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % LFP TABLET FRAMEs
    moveINDsF = readtable(moveINDid.(groups2id{groupIi}));
    lfpTABframes = table2array(moveINDsF(matches(moveINDsF.MoveType,'LFP TABLET'),...
        3:4));
    vidDURATION = round((lfpTABframes(2) - lfpTABframes(1))/60);
    outDATAFin.vidDURATION.(groups2id{groupIi}) = vidDURATION;

    tabletFRAMEstart.(groups2id{groupIi}) = lfpTABframes(1);
    tabletFRAMEend.(groups2id{groupIi}) = lfpTABframes(2);

    dlc_lab2use.(groups2id{groupIi}) = d1_labDAT.(groups2id{groupIi})...
        (tabletFRAMEstart.(groups2id{groupIi}):tabletFRAMEend.(groups2id{groupIi}),:);

    % totalNumSamplesV.(groups2id{groupIi}) = height(dlc_lab2use.(groups2id{groupIi}));
    % totalNumSecs.(groups2id{groupIi}) = totalNumSamplesV.(groups2id{groupIi})/60; % 60 fps

    streamOfInt.(primarySTN).(primarySTN).(groups2id{groupIi}) =...
        stnPRIME_Stream.(primarySTN).(groups2id{groupIi}).TS.TimeDomainData{1};

    % Original signal at 60 Hz
    ts_DLC.(groups2id{groupIi}) = 0:1/60:(height(dlc_lab2use.(groups2id{groupIi}))-1)/60;

    % Target sampling rate at 250 Hz
    ts_LFP.(groups2id{groupIi}) = 0:1/250:(height(streamOfInt.(primarySTN).(primarySTN).(groups2id{groupIi}))-1)/250;

    allColNs = dlc_lab2use.(groups2id{groupIi}).Properties.VariableNames;
    dlc_lab2use2int.(groups2id{groupIi}) = table;
    for coli = 1:width(dlc_lab2use.(groups2id{groupIi}))
        % Tmp col
        tmpCol = allColNs{coli};

        % Interpolation
        x_250A = interp1(ts_DLC.(groups2id{groupIi}),...
            transpose(dlc_lab2use.(groups2id{groupIi}).(tmpCol)),...
            ts_LFP.(groups2id{groupIi}), 'spline');

        dlc_lab2use2int.(groups2id{groupIi}).(tmpCol) = transpose(x_250A);
    end

    %%%% INTERPOLATE FACtor
    interPolfactor = length(ts_LFP.(groups2id{groupIi})) /...
        length(ts_DLC.(groups2id{groupIi}));

    trimFrame_int.(groups2id{groupIi}) =...
        linspace(tabletFRAMEstart.(groups2id{groupIi}),...
        tabletFRAMEend.(groups2id{groupIi}),...
        length(ts_LFP.(groups2id{groupIi})));

    % Load in and loop through movements
    moveINDsFcln = moveINDsF(~moveINDsF.BeginF == 0,:);
    restIItable = moveINDsFcln(matches(moveINDsFcln.MoveType,'REST'),:);
    moveIItable = moveINDsFcln(~matches(moveINDsFcln.MoveType,...
        {'REST','LFP TABLET'}),:);

    % FUNCTION

    % movKIN = [dlc_lab2use2int.B.fTip1_x(startINDi:endINDi) ,...
    %     dlc_lab2use2int.B.fTip1_y(startINDi:endINDi)];
    % movLFP = streamOfInt.Left.Left.GroupB(startINDi:endINDi);

    [lfpKINdata] = computeKIN_LFPcor(restIItable,moveIItable,...
        trimFrame_int.(groups2id{groupIi}),...
        streamOfInt.(primarySTN).(primarySTN).(groups2id{groupIi}),...
        dlc_lab2use2int.(groups2id{groupIi}));


    outDATAFin.tabletFRAMEstart.(groups2id{groupIi}) = tabletFRAMEstart.(groups2id{groupIi});
    outDATAFin.tabletFRAMEend.(groups2id{groupIi}) = tabletFRAMEend.(groups2id{groupIi});
    outDATAFin.dlc_lab2use.(groups2id{groupIi}) = dlc_lab2use.(groups2id{groupIi});
    outDATAFin.streamOfInt.(primarySTN).(primarySTN).(groups2id{groupIi}) =...
        streamOfInt.(primarySTN).(primarySTN).(groups2id{groupIi});
    outDATAFin.restIItable.(groups2id{groupIi}) = restIItable;
    outDATAFin.moveIItable.(groups2id{groupIi}) = moveIItable;
    outDATAFin.trimFrame_int.(groups2id{groupIi}) = trimFrame_int;
    outDATAFin.dlc_lab2use2int.(groups2id{groupIi}) = dlc_lab2use2int.(groups2id{groupIi});
    outDATAFin.lfpKINdata.(groups2id{groupIi}) = lfpKINdata;


end



% SAVE STUFF
cd(saveDIR)
saveNAME = [subjectID , '_' , stnHemi, '.mat'];
save(saveNAME,'outDATAFin');




end


%%%% SOME CURSORY ANALYSIS between BETA and MOVE KIN


%%%% SOME CURSORY ANALYSIS between GROUP A and GROUP B

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


% contact list
% BSLFP - LFPDATA = 2 Hz (0.5 s) steps - FFT (POWER) / mA
% Groups.Final(1).ProgramSettings.SensingChannel(1).ElectrodeState{1,1}.
% :::: Electrode , mA , fraction????? 
% EACH cell of ElectrodeState is a different contact - LOOK for everything
% BUT CASE 
% Sensing channel is rows for each hemisphere
% 




















%% FUNCTIONS with MAIN FUNCTION



function [stnL_BD_BS_testInfo , stnL_LFP_BS_KinLFP ,...
     stnL_LFP_St_KinLFP]...
    = getGroupInfo(allGroups, groupType , BS_LFP , LFPstreamTimes ,...
    LFPstreamTab , videoDirectory , moveIndFname)

switch groupType
    case 1
        group2use = 'GROUP_A';
    case 2
        group2use = 'GROUP_B';
    case 3
        group2use = 'GROUP_C';
end

% Get Group A rows
group_inds = matches(allGroups,group2use);
stnL_LFP_BS_gr = BS_LFP(group_inds,:); % BS LFP 
stnL_Ttims_gr = LFPstreamTimes(group_inds,:);
stnLSt_tab_LH_gr = LFPstreamTab(group_inds,:); % Stream LFP

% Compute length of VIDEO for CONDITION
% USE MOVEINDICES frames LFP TABLET
cd(videoDirectory)
moveINDsF = readtable(moveIndFname);
lfpTABframes = table2array(moveINDsF(matches(moveINDsF.MoveType,'LFP TABLET'),...
    3:4));
vidDURATION = round((lfpTABframes(2) - lfpTABframes(1))/60);

[~ , stnL_BD_BS_row_grp] = min(abs(vidDURATION - stnL_Ttims_gr.Duration));

% GROUP A - Condition TEST
stnL_BD_BS_testInfo = stnL_Ttims_gr(stnL_BD_BS_row_grp,:);
stnL_LFP_BS_KinLFP = stnL_LFP_BS_gr(stnL_BD_BS_row_grp,:);
stnL_LFP_St_KinLFP = stnLSt_tab_LH_gr(stnL_BD_BS_row_grp,:);


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






function  [lfpProcData] = computeKIN_LFPcor(restIItable,moveIItable,...
    trimFrame_int,streamOfInt,dlc_lab2use2int)

allREST = cell(height(restIItable),1);
for rii = 1:height(restIItable)

    tmpRrow = restIItable(rii,:);

    startINDr = tmpRrow.BeginF;
    endINDr = tmpRrow.EndF;


    [~ , startINDir] = min(abs(startINDr - trimFrame_int));
    [~ , endINDir] = min(abs(endINDr - trimFrame_int));

    restLFP = streamOfInt(startINDir:endINDir);

    allREST{rii} = restLFP;

end

allRESTstack = cell2mat(allREST);

[PxxREST,FxxREST] = pwelch(transpose(allRESTstack), hanning(250), 125, 256, 250, 'onesided');
uVp_REST = sqrt(PxxREST).*rms(hanning(250)).*sqrt(2).*2.*250/256;
uPv_RESTsm = smoothdata(uVp_REST,'gaussian',5);
upv50rest = FxxREST < 50;
uPvRESTfreq = FxxREST(upv50rest);
uPvRESTpow = uPv_RESTsm(upv50rest);


groupLFPAmp = []; 
groupKinAmp = [];
betaCWT = {}; 
betaENV = {};
moveXY = {};
moveEUC = {};
timeALL = {};
moveID = {};
obtItemID = {};
betapsdall = {};
betapsdall2 = {};
movepsdall = {};
movepsdall2 = {};

moveCOUNT = 1;

for mii = 1:height(moveIItable)

    close all
    tmpMOve = moveIItable(mii,:);

    if ~matches(tmpMOve.MoveType{1},{'FINGER TAP','HAND OC','HAND PS'})
        continue
    end

    startIND = tmpMOve.BeginF;
    endIND = tmpMOve.EndF;

    [~ , startINDi] = min(abs(startIND - trimFrame_int));
    [~ , endINDi] = min(abs(endIND - trimFrame_int));

    % offset ----- check 
    startINDi = startINDi + 5;
    endINDi = endINDi - 150;


    % FIGURE out OPTIMAL ITEM

    colNames = dlc_lab2use2int.Properties.VariableNames; %
    colNames2 = cellfun(@(x) split(x,'_'), colNames,...
        'UniformOutput',false);
    colNames3 = unique(cellfun(@(x) x{1}, colNames2,...
        'UniformOutput',false));
    colNames4 = colNames3(~matches(colNames3,'frames'));

    optItemInd = colNames4{tmpMOve.OptItem};

    xNAMEot = [optItemInd,'_x'];
    yNAMEot = [optItemInd,'_y'];

    movKIN1 = [dlc_lab2use2int.(xNAMEot)(startINDi:endINDi) ,...
        dlc_lab2use2int.(yNAMEot)(startINDi:endINDi)];

    [movKIN2] = getEucDist(dlc_lab2use2int.(xNAMEot)(startINDi:endINDi) ,...
        dlc_lab2use2int.(yNAMEot)(startINDi:endINDi));

    movLFP = streamOfInt(startINDi:endINDi);

    ts_LFPtmp = 0:1/250:(height(movLFP)-1)/250;

    tiledlayout(7,1,"TileSpacing","tight") %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xTime = ts_LFPtmp;
    pointXsm = smoothdata(movKIN1(:,1),'gaussian',5);
    % pointYsm = smoothdata(movKIN(:,2),'gaussian',5);
    
    nexttile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(xTime,pointXsm,'k');
    % hold on
    % plot(xTime,pointYsm,'r');
    xlim([0 max(ts_LFPtmp)])
    ylabel('X deflection')
    xlabel('Time in seconds')

    nexttile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pointX2sm = smoothdata(movKIN2(:,1),'gaussian',15);
    plot(xTime(1:length(pointX2sm)),pointX2sm,'r')
    xlim([0 max(ts_LFPtmp)])
    ylabel('EDist')
    xlabel('Time in seconds')


    [pxxMM,fMM] = pwelch(pointX2sm,250,125,250,250);
    p2bmm = pow2db(pxxMM);
    p2bMMs = smoothdata(p2bmm,'gaussian',2);
    fMMi = fMM > 1 & fMM < 10;
    fMMf = fMM(fMMi);
    pxxMf = p2bMMs(fMMi);


    nexttile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fs = 250;
    % T = 1/250;
    ecgClean = perceive_ecg(transpose(movLFP),Fs,0);
    cleanLFP = ecgClean.cleandata;

    [cfs,frq] = cwt(cleanLFP,Fs,'FrequencyLimits',[1 50],'VoicesPerOctave',24);
    t = 0:1/Fs:(length(cleanLFP)-1)/Fs;
    absCFS = abs(cfs);
    imagesc(t,frq,absCFS);
    axis xy
    ylabel('Frequency (Hz)')
    xlabel('Time in seconds')

    % test = 1;
    % figure;
    nexttile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pspectrum(cleanLFP,250,'FrequencyLimits',[0 50],'FrequencyResolution',3)
    [bPxx,bFxx] = pspectrum(cleanLFP,250,'FrequencyLimits',[0 50],'FrequencyResolution',3);
    bPxxP = pow2db(bPxx);

    [Pupvxx,upVfxx] = pwelch(transpose(movLFP), hanning(250), 125, 256, 250, 'onesided');
    uVp_t = sqrt(Pupvxx).*rms(hanning(250)).*sqrt(2).*2.*250/256;
    uPv_As = smoothdata(uVp_t,'gaussian',5);
    upv50 = upVfxx < 50;
    upVfxxU = upVfxx(upv50);
    uVp_tU = uPv_As(upv50);

    nexttile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pspectrum(pointX2sm,250,'FrequencyLimits',[1 10],'FrequencyResolution',0.75);
    [mPxx,mFxx] = pspectrum(pointX2sm,250,'FrequencyLimits',[1 10],'FrequencyResolution',0.75);
    mPxxP = pow2db(mPxx);


    nexttile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    betaIND = frq > 13 & frq < 30;
    betaBAND = absCFS(betaIND,:);

    betaBmean = mean(betaBAND,1);
    betaBsm = smoothdata(betaBmean,'gaussian',30);
    plot(xTime,betaBsm)
    xlim([0 max(ts_LFPtmp)])

    ylabel('Ave. beta power')
    xlabel('Time in seconds')

    nexttile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate instantaneous phase using Hilbert transform

    betaPass = bandpass(cleanLFP,[13 30],Fs);
    [betaUP,~] = envelope(betaPass);
    betaPassSm = smoothdata(betaUP,'gaussian',30);

    plot(xTime,betaPassSm);
    xlim([0 max(ts_LFPtmp)])
    ylabel('beta power')
    xlabel('Time in seconds')

    groupLFPAmp(moveCOUNT,1) = mean(betaPassSm);
    groupKinAmp(moveCOUNT,1) = max(abs(pointXsm));

    betaCWT{moveCOUNT,1} = betaBsm;
    betaENV{moveCOUNT,1} = betaPassSm;
    moveXY{moveCOUNT,1} = movKIN1;
    moveEUC{moveCOUNT,1} = movKIN2;
    timeALL{moveCOUNT,1} = xTime;
    moveID{moveCOUNT,1} = tmpMOve.MoveType{1};
    obtItemID{moveCOUNT,1} = optItemInd;

    betapsdall{moveCOUNT,1} = [bPxxP,bFxx];
    betapsdall2{moveCOUNT,1} = [uVp_tU,upVfxxU];
    movepsdall{moveCOUNT,1} = [mPxxP,mFxx];
    movepsdall2{moveCOUNT,1} = [pxxMf,fMMf];

    pause(1.5) 
    close all
    clc

    moveCOUNT = moveCOUNT + 1;

end

lfpProcData.betaCWT = betaCWT;
lfpProcData.betaENV = betaENV;
lfpProcData.moveXY = moveXY;
lfpProcData.moveEUC = moveEUC;
lfpProcData.timeALL = timeALL;
lfpProcData.moveID = moveID;
lfpProcData.obtItemID = obtItemID;
lfpProcData.groupLFPAmp = groupLFPAmp;
lfpProcData.groupKinAmp = groupKinAmp;
lfpProcData.betapsdall = betapsdall;
lfpProcData.movepsdall = movepsdall;
lfpProcData.betapsdall2 = betapsdall2;
lfpProcData.movepsdall2 = movepsdall2;
lfpProcData.REST = [uPvRESTpow,uPvRESTfreq];


end






function [euclidall] = getEucDist(tmpXdata , tmpYdata)

    labelData = [tmpXdata , tmpYdata];
    euclidall = zeros(height(labelData)-1,1);
    for frame_i = 1:height(labelData)
        if frame_i ~= height(labelData)
            point1 = labelData(frame_i,:);
            point2 = labelData(frame_i + 1,:);
            euclidall(frame_i) = pdist2(point1 , point2);
        end
    end

end