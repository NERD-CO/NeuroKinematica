pcName = getenv("COMPUTERNAME");

switch pcName
    case 'DESKTOP-I5CPDO7'
        mainDir = 'D:\Dropbox\ConferencePresentationsMeta\INS_2024\MDT2';
    case 'DESKTOP-EGFQKAI'
        mainDir = 'C:\Users\johna\Dropbox\ConferencePresentationsMeta\INS_2024\MDT2';

    otherwise


end


%%

% loop through and store as structs
leftSTNLoc = [mainDir,'\Streaming\Left STN'];
leftSTNfile = 'Report_Json_Session_Report_20230908T131441.json';

rightSTNLoc = [mainDir,'\Streaming\Right STN'];
rightSTNfile = 'Report_Json_Session_Report_20230911T133035.json';

%%
% cd(leftSTNLoc)
leftFF = [leftSTNLoc , filesep , leftSTNfile];
rightFF = [rightSTNLoc , filesep , rightSTNfile];

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

%%

% first test '11:56:27 AM MDT'
% second test '12:11:18 PM MDT'




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
