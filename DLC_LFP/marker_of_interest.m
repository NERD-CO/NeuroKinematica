%% Directory set-up - Navigate b/t machines
pcname = getenv('COMPUTERNAME');

switch pcname
    case 'DESKTOP-I5CPDO7'   %%% JAT Desktop

        % mainDir = '';

    case 'DSKTP-JTLAB-EMR'   %%% ER Desktop

        mainDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\DLC_LFP';

    case 'NSG-M-FQBPFK3'     %%% ER PC

        mainDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Clinical\DLC_LFP';
end

%% Analyze data isolated by casedate and hemisphere

% Define casedate_hem
casedate_hem = '09_12_2023_LSTN_v2';
% casedate_hem = '09_12_2023_RSTN_v2';

mainDir2 = [mainDir , filesep , casedate_hem];

cd(mainDir2)

%% DLC Movement data folders

DLC_csv_dir = [mainDir2 , filesep , 'csv folder']; % contains dlc label timeseries data as csv files
DLC_mat_dir = [mainDir2 , filesep , 'mat folder']; % contains dlc label timeseries data as mat files
DLC_video_dir = [mainDir2 , filesep , 'video folder', filesep , 'Converted for dual GUI']; % contains labeled videos + MovementIndex csv


%% Isolate dlc outputs of interest

% Generate list of dlc-video-labeled CSV files
% cd(DLC_csv_dir)
% mainCSV = DLC_csv_dir('*.csv');
% mainCSV2 = {mainCSV.name};

% Generate list of dlc-video-labeled MAT files
cd(DLC_mat_dir)
mainMat = dir('*.mat');
mainMAT2 = {mainMat.name};

% Generate list of Motor Index CSVs
cd(DLC_video_dir)
moveCSV = dir('*.csv');
moveCSV2 = {moveCSV.name};


%% Plot dlc label / anatomical marker traces

% Loop through dlc mat files corresponding with move_CSV files
moveCSV_fileLocation = DLC_video_dir;
cd(moveCSV_fileLocation)

for csv_i = 1:length(moveCSV2)

    tmpCSV = moveCSV2{csv_i};

    % Split file names to extract relevant parts (dateID, sessID, and hemID)
    nameParts = split(tmpCSV,'_');
    dateID = nameParts{1};
    sessID = nameParts{3};
    matName_title = [dateID , '-' , sessID];

    % Find and load corresponding dlcDAT MAT file
    cd(DLC_mat_dir)
    matTempfind = [dateID , '_' , sessID];
    matInd = contains(mainMAT2 , matTempfind);
    matName = mainMAT2{matInd};
    load(matName)

    %load(moveCSV_fileName , 'outDATA')

    colNames = outDATA.Properties.VariableNames; %
    colNames2 = cellfun(@(x) split(x,'_'), colNames,...
        'UniformOutput',false);
    colNames3 = unique(cellfun(@(x) x{1}, colNames2,...
        'UniformOutput',false));
    colNames4 = colNames3(~matches(colNames3,'frames'));

    % Initialize 'euclidall' to store Euclidean distances between successive points
    euclidall = zeros(height(outDATA)-1,length(colNames4));

    % Iterate over each label and compute Euclidean distance for each frame
    for label_i = 1:length(colNames4)

        tmpLabel_x = [colNames4{label_i} , '_x'];
        tmpLabel_y = [colNames4{label_i} , '_y'];

        tmpXdata = outDATA.(tmpLabel_x);
        tmpYdata = outDATA.(tmpLabel_y);

        labelData = [tmpXdata , tmpYdata];

        for frame_i = 1:height(labelData)
            if frame_i ~= height(labelData)
                point1 = labelData(frame_i,:);
                point2 = labelData(frame_i + 1,:);
                euclidall(frame_i , label_i) = pdist2(point1 , point2);
            end
        end
    end

    euclidalltmp = euclidall;

    mall = mean(mean(euclidalltmp));
    sall = std(std(euclidalltmp));
    artThresh = mall + (sall*2);

    euclidalltmp2 = euclidalltmp;
    euclidalltmp2(euclidalltmp > artThresh) = 0;

    rangeAll = zeros(1,width(euclidalltmp2));
    for rrii = 1:width(euclidalltmp2)
        rangeAll(rrii) = range(euclidalltmp2(:,rrii));
    end

    % Create y data matrix
    yPLOTDATA = euclidalltmp2;
    pCURRENT = 0;
    yLABELSi = zeros(width(euclidalltmp2),1);
    yLABELSn = cell(width(euclidalltmp2),1);

    for piip = 1:width(euclidalltmp2)
        if piip ~= 1
            pBump = rangeAll(piip-1);
            pCURRENT = pCURRENT + pBump;
        end

        tmpDATA = euclidalltmp2(:,piip) + pCURRENT;
        yPLOTDATA(:,piip) = tmpDATA;

        yLABELSi(piip) = pCURRENT;
        yLABELSn{piip} = ['Marker ', num2str(piip)];
    end

    % Store plot elements
    MarkerPlotData = yPLOTDATA;
    MarkerYlabels = yLABELSn;
    MarkerYindicies = yLABELSi;

    plot(yPLOTDATA,'Color',[0 0 0 0.5])

    xlim([1 height(yPLOTDATA)])

    % fix y lim
    ylim([0 max(yPLOTDATA,[],'all')+5])
    
    % fix labels
    yticks(MarkerYindicies)
    yticklabels(MarkerYlabels)

end