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

%%
stnLEFTjson  = jsondecode(fileread(leftFF));
stnRIGHTjson = jsondecode(fileread(rightFF));

%%
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
%%

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

% Get Group A rows
groupA_inds = matches(groupIDs,'GROUP_A');
stnL_LFP_BS_grA = stnLEFTBSlfp(groupA_inds,:);
stnL_LFP_BStms_grA = BSlfptimes(groupA_inds,:);
stnL_Ttims_grA = stnLEFTsTtimes(groupA_inds,:);

% Pull in video durations

% Load Video info XslX
stnL_BD_R_vinfo = readtable(leftSTN_rightBD_Vinfo_fname);
% Pull out Group A and Group B
groupA_Vinfo = stnL_BD_R_vinfo(contains(stnL_BD_R_vinfo.ftype,'A'),:);

[~ , stnL_BD_BS_row] = min(abs(groupA_Vinfo.duration - stnL_Ttims_grA.Duration));

stnL_BD_BSdata = stnL_Ttims_grA(stnL_BD_BS_row,:);

GROUPA_Video_CSV = groupA_Vinfo.csvName;
allVideoNames = dir('*.mp4');
allVideoNames2 = {allVideoNames.name};
allVideoNames3 = cellfun(@(x) [extractBefore(x,'_labeled'),'.mp4'],...
    allVideoNames2, 'UniformOutput',false);
GROUPA_Video_Name = allVideoNames2(contains(allVideoNames3,...
    extractBefore(GROUPA_Video_CSV{1},'.')));
% USE stnL_BD_BSdata to FIND!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 
%%% USE ERIN CODE TO GET CSV data
%%% USE OJEMANN CODE TO GET FRAMES for LFP alignment
GROUPA_LFP = 1;

















% contact list
% BSLFP - LFPDATA = 2 Hz (0.5 s) steps - FFT (POWER) / mA
% Groups.Final(1).ProgramSettings.SensingChannel(1).ElectrodeState{1,1}.
% :::: Electrode , mA , fraction????? 
% EACH cell of ElectrodeState is a different contact - LOOK for everything
% BUT CASE 
% Sensing channel is rows for each hemisphere
% 



%%




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












