pcName = getenv("COMPUTERNAME");

switch pcName
    case 'DESKTOP-I5CPDO7'



%%

% loop through and store as structs
leftSTNLoc = 'C:\Users\johna\Dropbox\ConferencePresentationsMeta\INS_2024\MDT2\Streaming\Left STN';
leftSTNfile = 'Report_Json_Session_Report_20230908T131441.json';

rightSTNLoc = 'C:\Users\johna\Dropbox\ConferencePresentationsMeta\INS_2024\MDT2\Streaming\Right STN';
rightSTNfile = 'Report_Json_Session_Report_20230911T133035.json';

%%

leftFF = [leftSTNLoc , filesep , leftSTNfile];
rightFF = [rightSTNLoc , filesep , rightSTNfile];

%%
stnLEFTjson  = jsondecode(fileread(leftFF));
stnRIGHTjson = jsondecode(fileread(rightFF));

%%
stnLEFTstream = stnLEFTjson.BrainSenseTimeDomain;

% process with a function

%%

% The input date-time string
dateTimeStr = '2023-09-08T17:47:31.000Z';

% Convert the string to a datetime object
dateTimeObj = datetime(dateTimeStr, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z');

% Set the time zone to UTC (Coordinated Universal Time) since the input is in UTC
dateTimeObj.TimeZone = 'UTC';

% Convert to Mountain Time
dateTimeObj_Mountain = datetime(dateTimeObj, 'TimeZone', 'America/Denver');

% Extract the time component in AM/PM format
timeComponent_AMPM = datestr(dateTimeObj_Mountain, 'HH:MM:SS AM');

% Display the time component
disp(['Time in Mountain Time (AM/PM): ', timeComponent_AMPM]);
