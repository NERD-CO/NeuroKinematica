function [] = dlc_processCSV2(inPUTS)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


arguments
    inPUTS.userLOC (1,1) string = "noLOC"
    inPUTS.saveLOC (1,1) string = "noLOC"
    inPUTS.depNUM (1,1) double = 1;
    inPUTS.USERid (1,:) char = 'ER';
end


% DATA DIRECTORY LOCATION

if strcmp(inPUTS.userLOC,"noLOC")
    useFile = uigetfile([],'SELECT LOCATION OF LABEL CSV');
else
    useFile = inPUTS.userLOC;
end


% SAVE DIRECTORY LOCATION
if strcmp(inPUTS.saveLOC,"noLOC")
    saveDIR = uigetdir([],'SELECT LOCATION TO SAVE');
else
    saveDIR = inPUTS.saveLOC;
end

% cd(doLOC)

% [CSV_list] = getFiles(doLOC , 'csv');

% hdRT = getHDRinfo(CSV_list{1});

% outDATA = struct;

% for ci = 1:length(CSV_list)

% Read in frame data
rawCell = getRAWdat(useFile);

[rawTAB] = datCellTab(rawCell , useFile);

% tmpHdr = getHDRinfo(useFile);

% outDATA.hdr.(tmpHdr.camera) = tmpHdr.camera;
outDATA = rawTAB;

% end

% SAVE LOCATION
cd(saveDIR);
% Save name
% Date of surgery, depth, and recording number
switch inPUTS.USERid
    case 'ER'

        % date
        % depth
        % trial ID
        % camera
        firstParse = split(useFile,'\');
        fileNAME = firstParse{length(firstParse)};
        secParse = split(fileNAME,'_');
        date2use = secParse{1};
        depth2use = secParse{3}; 
        trial2use = secParse{2};
        camera2use = secParse{5};
        folder2use = secParse{4};

        saveName = ['dlcDAT_',date2use,'_',depth2use,'_',trial2use,'_',...
            camera2use,'_',folder2use,'.mat'];

    otherwise
        saveName = ['dlcDAT_',hdRT.date,'_',hdRT.depth,'_','R',num2str(inPUTS.depNUM),'.mat'];
end
% Save file
save(saveName,'outDATA');


end





% ----------------------------------------------------------------------- %
% Obtain list of CSV files
% ----------------------------------------------------------------------- %
function [out_list] = getFiles(inLOC , ftype)

cd(inLOC)

fType2f = ['*.',ftype];

tDIR = dir(fType2f);

out_list = {tDIR.name};


end




% ----------------------------------------------------------------------- %
% Extract info of interest from CSV
% ----------------------------------------------------------------------- %
function [outDatCell] = getRAWdat(inCVSlist)

fid = fopen(inCVSlist);
tline = fgetl(fid);
rawTab = cell(1000,1);
tCount = 1;
while ischar(tline)
    rawTab{tCount} = tline;
    tline = fgetl(fid);
    tCount = tCount + 1;
end
fclose(fid);


outDatCell = rawTab(cellfun(@(x) ~isempty(x), rawTab, 'UniformOutput', true));


end


% ----------------------------------------------------------------------- %
% Clean up cell array and convert to Table
% ----------------------------------------------------------------------- %
function [outDataTable] = datCellTab(inCVSlist , csvTname)

varNames1 = strsplit(inCVSlist{2},',');
varNames2 = strsplit(inCVSlist{3},',');

varNames1 = varNames1(2:end);
varNames2 = varNames2(2:end);

% Create column headings
varNamesF = cellfun(@(x,y) [x ,'_',y], varNames1, varNames2,'UniformOutput', false);

rowTemp = inCVSlist(4:end);
rowIDt1 = cellfun(@(x) strsplit(x,',') , rowTemp , 'UniformOutput' , false);
rowIDt2 = cellfun(@(x) x{1}, rowIDt1 , 'UniformOutput', false);
% Create frame index
frameNUM = cellfun(@(x) str2double(x) , rowIDt2 , 'UniformOutput' , true);

% Create data
tmpDattab = readtable(csvTname);
datARRy = table2array(tmpDattab(:,2:end));

% Create output table

outDataTable = array2table(datARRy,'VariableNames',varNamesF);
outDataTable.frames = frameNUM;


end


% ----------------------------------------------------------------------- %
% Get Header data
% ----------------------------------------------------------------------- %
function [hDr] = getHDRinfo(csvTname)

% Header file
hdRt1 = strsplit(csvTname,{'_','-','.'});

hDr.date = hdRt1{1};
hDr.depth = hdRt1{2};
hDr.camera = hdRt1{4};

end
