clear
clc
close all
dataLOC = 'C:\Users\johna\Dropbox\PowerPoint_Meta\2024_INS_Vancouver\Data\videoFRAMEfig';
cd(dataLOC);

% Load in Video file
% videoFile = 't1_20230816_idea08_session002_rightCam-0000DLC_resnet50_INS_2024_MPR8_LSTNMar20shuffle1_100000_labeled.mp4';

videoMat = 'MDT8_T1_videoCdataFT.mat';
load(videoMat)
%%

% dlc_lab_vidObj = VideoReader(videoFile);
% dlc_lab_vid = struct('cdata',zeros(dlc_lab_vidObj.Height,dlc_lab_vidObj.Width,3,'uint8'),'colormap',[]);
% 
% frami = 1;
% while hasFrame(dlc_lab_vidObj)
%     dlc_lab_vid(frami).cdata = readFrame(dlc_lab_vidObj);
%     frami = frami+1;
% end
% dlcLabTab_vid = struct2table(dlc_lab_vid);
% disp('Video1 done!')

%%
% Index CSV
% indexCSV = 't1_20230816_idea08_session002_rightCam-0000DLC_resnet50_INS_2024_MPR8_LSTNMar20shuffle1_100000_labeled_MoveIndex.csv';

% csvIND = readtable(indexCSV);

%%

% dlc_lablab2useVIDE = dlcLabTab_vid(startIND:endIND,:);
% save("MDT8_GroupA_Test1_Finger.mat",'dlc_lablab2useVIDE','-append')

% FT_table = csvIND(find(matches(csvIND.MoveType,'FINGER TAP'),1,'first'),:);
% FT_startIND = FT_table.BeginF;
% FT_endIND = FT_table.EndF;
% 
% dlc_lab2useVIDE = dlcLabTab_vid(FT_startIND:FT_endIND,:);
% save("MDT8_T1_videoCdataFT.mat",'dlc_lab2useVIDE','FT_startIND','FT_endIND','-v7.3')

%%

% Load in dlcKinematic table
dlcMAT = 'dlcDAT_t1_idea08_session002_LSTN_RBODY.mat';
load(dlcMAT)

%%

frameID = FT_startIND:FT_endIND;


% Create a blank RGB matrix with same Y-dimensions, but 5000 columns NANs
blankRGB = ones(size(dlc_lab2useVIDE.cdata{1},1),7000,3,"uint8");

% Extract just the x and y points from the table
dlcTAB = outDATA.labelTab.session002;
colNames = dlcTAB.Properties.VariableNames;
colNames2 = contains(colNames,{'_x','_y'});
dlcTAB2 = dlcTAB(:,colNames2);
allXtab = dlcTAB(:,contains(colNames,'_x'));

frameXmin = zeros(length(frameID),1);
frameXmax = zeros(length(frameID),1);
% For each row find the min and max coordinates
for fi = 1:70%length(frameID)

    tmpF = frameID(fi);
    tmpMIN = min(table2array(allXtab(tmpF,:)));
    tmpMAX = max(table2array(allXtab(tmpF,:)));

    frameXmin(fi) = tmpMIN;
    frameXmax(fi) = tmpMAX;

    tmpSlice = dlc_lab2useVIDE.cdata{fi}(:,tmpMIN:tmpMAX,:);

    % find next NaN
    nanNloc = find(blankRGB(1,:,1)==1,1,'first');
    frameLEN = length(tmpMIN:tmpMAX)-1;
    blankRGB(:,nanNloc:nanNloc+frameLEN,:) = tmpSlice;

    


end


% Add to NaN matrix at first NaN location based on 

% 