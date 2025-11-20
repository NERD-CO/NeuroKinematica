function MasterZETA_MUA = aggregate_ZETA_MUA_and_plot(FR_Kin_Dir, varargin)

% aggregate_ZETA_MUA_and_plot
%
% MUA analogue of aggregate_ZETA_and_plot.
% Scans case folders for *_ZETA_Summary_MUA.csv in "Zeta Testing MUA",
% attaches vecD_MUA / vecTime_MUA (and optional IFR/PSTH_MUA),
% builds MasterZETA_MUA and basic scatter plots, then saves:
%   MasterZETA_MUA_AllSubjects.csv
%   MasterZETA_MUA_AllData.mat

p = inputParser;
p.addParameter('IO_DataDir','', @(s) ischar(s) || isstring(s));
p.addParameter('SaveDir','', @(s) ischar(s) || isstring(s));
p.addParameter('ZetaCsvNamePattern','*_ZETA_TS_Summary_MUA*.csv', @(s)ischar(s)||isstring(s));
p.addParameter('SigZ', 2, @(x) isscalar(x) && x>0);
p.addParameter('SigP', 0.05, @(x) isscalar(x) && x>0 && x<1);
p.addParameter('DepthMap', containers.Map({'t','c','b'},{'dorsal STN','central STN','ventral STN'}));
p.addParameter('PrettyMoveMap', containers.Map( ...
    {'HAND OC','HAND PS','ARM EF','REST'}, {'Hand OC','Hand PS','Arm EF','Rest'}));
p.addParameter('YMax', 5, @(x) isscalar(x) && x>0);
p.parse(varargin{:});
U = p.Results;

% Output directory
if isempty(U.SaveDir)
    groupOut = fullfile(FR_Kin_Dir, 'Aggregate Zeta MUA Plots');
else
    groupOut = char(U.SaveDir);
end
if ~exist(groupOut,'dir'), mkdir(groupOut); end

%% 1) Find MUA ZETA summary CSVs

fileList = find_ZetaSummary_MUA_files(FR_Kin_Dir, U.ZetaCsvNamePattern);

%% 2) Load & unify ZETA MUA rows

MasterZETA_MUA = table();
canonOrder = {'CaseDate','Hemisphere','MUA_Field','MoveType','Depth', ...
    'nTrials','UseMaxDur_s', ...
    'ZetaP_MUA','ZetaZ_MUA','ZetaD_MUA', ...
    'ZetaTime_MUA','Zeta_Idx', ...
    'MeanStimDur_s','StdStimDur_s', ...
    'MeanZ_MUA','MeanP_MUA'};

for k = 1:numel(fileList)
    csvPath = fileList(k).fullpath;
    Traw = readtable(csvPath, 'TextType','string');
    caseID = string(fileList(k).subject);
    hemID  = string(fileList(k).hemi);

    % Map the core MUA columns
    T = table();
    T.MUA_Field      = string(Traw.MUA_Field);
    T.MoveType       = string(Traw.MoveType);
    T.Depth          = string(Traw.Depth);
    T.nTrials        = double(Traw.nTrials);
    T.UseMaxDur_s    = double(Traw.UseMaxDur_s);
    T.ZetaP_MUA      = double(Traw.ZetaP_MUA);
    T.ZetaZ_MUA      = double(Traw.ZetaZ_MUA);
    T.ZetaD_MUA      = double(Traw.ZetaD_MUA);
    T.ZetaTime_MUA   = double(Traw.ZetaTime_MUA);
    T.Zeta_Idx       = double(Traw.Zeta_Idx);
    T.MeanStimDur_s  = double(Traw.MeanStimDur_s);
    T.StdStimDur_s   = double(Traw.StdStimDur_s);
    T.MeanZ_MUA      = double(Traw.MeanZ_MUA);
    T.MeanP_MUA      = double(Traw.MeanP_MUA);

    % Add CaseDate & Hemisphere
    T.CaseDate   = repmat(caseID, height(T), 1);
    T.Hemisphere = repmat(hemID,  height(T), 1);

    % Reorder columns to canonical order
    T = T(:, canonOrder);

    % -------------------------------------------------------------
    % Attach vecD_MUA / vecTime_MUA from *_ZETA_AllOutputs_MUA.mat
    % -------------------------------------------------------------
    ZETA_vecD_MUA   = cell(height(T),1);
    ZETA_vecT_MUA   = cell(height(T),1);
    ZETA_new_vecT_MUA  = cell(height(T),1);

    zetaFolderMUA = fileparts(csvPath);  % ...\Case\Zeta Testing MUA or ...\Case\LSTN\Zeta Testing MUA
    zetaAll_mat_MUA = fullfile(zetaFolderMUA, sprintf('%s_ZETA_TS_AllOutputs_MUA.mat', caseID));
    hasZetaMat = exist(zetaAll_mat_MUA,'file')==2;

    zetaMap = containers.Map;
    if hasZetaMat
        S = load(zetaAll_mat_MUA, 'ZETA_TS_Summary_MUA','all_sZETA_ts_MUA');
        ZetaTS_sumMUA = S.ZETA_TS_Summary_MUA;
        all_sZETA_MUA_case = S.all_sZETA_ts_MUA;

        for ii = 1:height(ZetaTS_sumMUA)
            key = sprintf('%s|%s|%s', ...
                string(ZetaTS_sumMUA.MUA_Field(ii)), ...
                string(ZetaTS_sumMUA.MoveType(ii)), ...
                string(ZetaTS_sumMUA.Depth(ii)));
            if ~isKey(zetaMap, key)
                zetaMap(key) = ii;
            end
        end
    else
        fprintf('[WARN] ZETA_TS_AllOutputs_MUA.mat not found for %s (hem=%s)\n', caseID, hemID);
    end


    % -------------------------------------------------------------
    % Optional: attach IFR/PSTH_MUA if you saved it
    PSTH_TimeCenters_MUA = cell(height(T),1);
    PSTH_Hz_MUA          = cell(height(T),1);
    IFR_Time_s_MUA       = cell(height(T),1);
    IFR_Hz_MUA           = cell(height(T),1);

    % Assumes (optionally) a file:
    %   ...\Zeta Testing MUA\IFR_PSTH_MUA\<CaseID>_IFR_PSTH_MUA_All.mat
    ifrFolderMUA   = fullfile(zetaFolderMUA, 'IFR_PSTH_MUA');  % adjust if you used a different name
    ifrAll_mat_MUA = fullfile(ifrFolderMUA, sprintf('%s_IFR_PSTH_MUA_All.mat', caseID));
    hasIFRMat      = exist(ifrAll_mat_MUA,'file')==2;

    ifrMap = containers.Map;
    IFRsumMUA   = [];
    all_IFR_MUA = [];

    if hasIFRMat
        S2 = load(ifrAll_mat_MUA);
        if isfield(S2,'IFR_PSTH_MUA_Summary') && isfield(S2,'all_IFR_MUA')
            IFRsumMUA   = S2.IFR_PSTH_MUA_Summary;
            all_IFR_MUA = S2.all_IFR_MUA;

            for ii = 1:height(IFRsumMUA)
                key = sprintf('%s|%s|%s', ...
                    string(IFRsumMUA.SpikeField(ii)), ...
                    string(IFRsumMUA.MoveType(ii)), ...
                    string(IFRsumMUA.Depth(ii)));
                if ~isKey(ifrMap, key)
                    ifrMap(key) = ii;
                end
            end
        else
            % File exists but doesn't have the expected IFR MUA variables
            warning('IFR_PSTH_MUA_Summary not found in %s; skipping IFR/PSTH_MUA attachment.', ifrAll_mat_MUA);
            hasIFRMat = false;
        end
    end
    % -------------------------------------------------------------


    % Per-row attachment
    for r = 1:height(T)
        key = sprintf('%s|%s|%s', ...
            string(T.MUA_Field(r)), string(T.MoveType(r)), string(T.Depth(r)));

        % ZETA vecD / vecTime
        if hasZetaMat && isKey(zetaMap, key)
            idxZ = zetaMap(key);
            sZ   = all_sZETA_MUA_case{idxZ};

            if isfield(sZ,'vecD') && ~isempty(sZ.vecD)
                ZETA_vecD_MUA{r} = sZ.vecD(:)';    % row vector
            else
                ZETA_vecD_MUA{r} = [];
            end

            if isfield(sZ,'vecTime') && ~isempty(sZ.vecTime)
                ZETA_vecT_MUA{r} = sZ.vecTime(:)'; % row vector
            else
                ZETA_vecT_MUA{r} = [];
            end

            if isfield(sZ,'new_vecTime_MUA') && ~isempty(sZ.new_vecTime_MUA)
                ZETA_new_vecT_MUA{r} = sZ.new_vecTime_MUA(:)'; % row vector
            else
                ZETA_new_vecT_MUA{r} = [];
            end

        else
            ZETA_vecD_MUA{r}  = [];
            ZETA_vecT_MUA{r} = [];
            ZETA_new_vecT_MUA{r} = [];
        end

        % IFR/PSTH MUA (optional)
        if hasIFRMat && isKey(ifrMap, key)
            idxI = ifrMap(key);

            % deal with cell vs non-cell gracefully
            cCenters = IFRsumMUA.PSTH_TimeCenters_s(idxI);
            if iscell(cCenters), cCenters = cCenters{1}; end

            cPSTH = IFRsumMUA.PSTH_Hz(idxI);
            if iscell(cPSTH), cPSTH = cPSTH{1}; end

            cIFRt = IFRsumMUA.IFR_Time_s(idxI);
            if iscell(cIFRt), cIFRt = cIFRt{1}; end

            cIFRr = IFRsumMUA.IFR_Hz(idxI);
            if iscell(cIFRr), cIFRr = cIFRr{1}; end

            PSTH_TimeCenters_MUA{r} = cCenters(:)';  % row
            PSTH_Hz_MUA{r}          = cPSTH(:)';
            IFR_Time_s_MUA{r}       = cIFRt(:)';
            IFR_Hz_MUA{r}           = cIFRr(:)';
        else
            PSTH_TimeCenters_MUA{r} = [];
            PSTH_Hz_MUA{r}          = [];
            IFR_Time_s_MUA{r}       = [];
            IFR_Hz_MUA{r}           = [];
        end
    end

    % Add dynamic columns
    T.vecD_MUA          = ZETA_vecD_MUA;
    T.vecTime_MUA       = ZETA_vecT_MUA;
    T.new_vecTime_MUA   = ZETA_new_vecT_MUA;
    T.PSTH_TimeCenters_s_MUA = PSTH_TimeCenters_MUA;
    T.PSTH_Hz_MUA       = PSTH_Hz_MUA;
    T.IFR_Time_s_MUA    = IFR_Time_s_MUA;
    T.IFR_Hz_MUA        = IFR_Hz_MUA;

    % Append this case
    MasterZETA_MUA = [MasterZETA_MUA; T];
end

if isempty(MasterZETA_MUA)
    warning('No ZETA MUA summary rows found. Nothing to plot.');
    return
end


%% 3) Attach Subject + Hemisphere labels from spreadsheet (re-use existing meta logic)
% just replace "MasterZETA" with "MasterZETA_MUA" inside it.

LabelMapCH = [];   % containers.Map for "CaseFolder|HemTag" -> label
NumMapCH   = [];   % containers.Map for "CaseFolder|HemTag" -> SubjectNum (string)
LabelMapC  = [];   % containers.Map for "CaseFolder" -> label (only when unambiguous)
NumMapC    = [];   % containers.Map for "CaseFolder" -> SubjectNum (string)
HemMapC    = [];   % containers.Map for "CaseFolder" -> HemTag (string)
haveMeta   = false;

if ~isempty(U.IO_DataDir)
    metaPath = fullfile(char(U.IO_DataDir),'Subject_Hem_MetaSummary.xlsx');
    if exist(metaPath,'file')
        M = readtable(metaPath,'TextType','string');

        % ---- Normalize / derive HemTag ----
        if ~ismember('HemTag', M.Properties.VariableNames) || any(M.HemTag=="")
            M.HemTag = strings(height(M),1);
            M.HemTag(endsWith(M.MoveCaseFolder,"_LSTN")) = "LSTN";
            M.HemTag(endsWith(M.MoveCaseFolder,"_RSTN")) = "RSTN";
        end
        % Be defensive about random whitespace / case
        M.CaseFolder = strtrim(M.CaseFolder);
        M.HemTag     = upper(strtrim(M.HemTag));
        M.SubjectNum = string(M.SubjectNum);

        % ---- Map 1: Case+Hem -> label / SubjectNum ----
        keyCH   = M.CaseFolder + "|" + M.HemTag;
        lblCH   = "Subject " + M.SubjectNum + " - " + M.HemTag;
        [ukCH, iaCH] = unique(keyCH,'stable');    % keep first occurrence
        LabelMapCH = containers.Map(ukCH, lblCH(iaCH));
        NumMapCH   = containers.Map(ukCH, M.SubjectNum(iaCH));

        % ---- Map 2: Case (only where unambiguous in spreadsheet) ----
        % For unilateral cases (appears once) we can resolve without hem on the left.
        [g, caseCats] = findgroups(M.CaseFolder);
        counts = splitapply(@numel, M.CaseFolder, g);
        unilateralCases = caseCats(counts==1);
        if ~isempty(unilateralCases)
            % Build case-only label/num/hem maps only for those unilateral cases
            umask = ismember(M.CaseFolder, unilateralCases);
            uCase = M.CaseFolder(umask);
            uHem  = M.HemTag(umask);
            uNum  = M.SubjectNum(umask);
            uLbl  = "Subject " + uNum + " - " + uHem;
            [uCuniq, iu] = unique(uCase,'stable');
            LabelMapC = containers.Map(uCuniq, uLbl(iu));
            NumMapC   = containers.Map(uCuniq, uNum(iu));
            HemMapC   = containers.Map(uCuniq, uHem(iu));
        end

        haveMeta = true;
    else
        warning('Meta spreadsheet not found at: %s. Falling back to generic subject labels.', metaPath);
    end
end

% ---- Build per-row labels for MasterZETA ----
caseF = string(strtrim(MasterZETA_MUA.CaseDate));     % case folder in ZETA rows
hemi  = upper(string(strtrim(MasterZETA_MUA.Hemisphere)));  % 'LSTN'/'RSTN' or ""

rowKeyCH = caseF + "|" + hemi;

prettyLabel = strings(height(MasterZETA_MUA),1);
subjNumStr  = strings(height(MasterZETA_MUA),1);
hemiFilled  = hemi;   % will fill from meta when missing

% Fallback generator uses stable per-case index so bilateral days stay split
[~,~,caseStableIdx] = unique(caseF,'stable');

for i = 1:height(MasterZETA_MUA)
    if haveMeta && hemi(i)~="" && isKey(LabelMapCH, rowKeyCH(i))
        % Exact case+hem match (bilateral days)
        prettyLabel(i) = LabelMapCH(rowKeyCH(i));
        subjNumStr(i)  = NumMapCH(rowKeyCH(i));
        % hemiFilled(i) already equals hemi(i)
    elseif haveMeta && isKey(LabelMapC, caseF(i))
        % Unilateral case (hem missing on left; safe to resolve by case only)
        prettyLabel(i) = LabelMapC(caseF(i));
        subjNumStr(i)  = NumMapC(caseF(i));
        hemiFilled(i)  = HemMapC(caseF(i));   % fill hemisphere for completeness
    else
        % Generic fallback (no spreadsheet or unmatched key)
        prettyLabel(i) = "Subject " + string(caseStableIdx(i)) + " - " + hemi(i);
        subjNumStr(i)  = string(caseStableIdx(i));
    end
end

% Store for plotting/export
MasterZETA_MUA.PrettyLabel     = prettyLabel;   % e.g., "Subject 1 - LSTN"
MasterZETA_MUA.SubjectNum      = subjNumStr;    % "1","2",...
MasterZETA_MUA.HemisphereFilled = hemiFilled;   % 'LSTN'/'RSTN' (filled for unilateral)

% sanity summary
unmatched = MasterZETA_MUA(~startsWith(MasterZETA_MUA.PrettyLabel,"Subject "),:);
fprintf('[Map check] unique x-labels: %d | unmatched rows: %d\n', ...
    numel(unique(MasterZETA_MUA.PrettyLabel,'stable')), height(unmatched));


%% 3a) All-categories-in-one scatter (Subjects on x, all MoveType×Depth ZetaZ scores stacked on y with xjitter)

% ---- build friendly labels (global, consistent across whole MasterZETA_MUA) ----
% Subject index in encounter order (Subject 1, Subject 2, ...)
[~, ~, subjIdxAll] = unique(MasterZETA_MUA.CaseDate,'stable');
subjNum = subjIdxAll;  % numeric 1..N

% % ---- build labels from SubjectNum (+ hemisphere, if present) ----
% hasHemi  = ~(MasterZETA_MUA.Hemisphere=="" | ismissing(MasterZETA_MUA.Hemisphere));
% hemiText = strings(height(MasterZETA_MUA),1);
% hemiText(hasHemi) = "–" + string(MasterZETA_MUA.Hemisphere(hasHemi));   % '–LSTN'/'–RSTN'
% % % don't display hemisphere for now
% % hemiText(:) = "";

% Build x-axis groups from meta-driven labels
prettyPerRow = MasterZETA_MUA.PrettyLabel; % per-row label like "Subject 1 - LSTN"

% x positions (with jitter)
[uniqSubsAll, ~, subjIdxX] = unique(prettyPerRow,'stable'); % unique x-tick strings in encounter order
x_all  = subjIdxX; % x position index for each row
jitter_all = (rand(size(x_all)) - 0.5) * 0.30; % ±0.15 jitter
xj_all = x_all + jitter_all;


% Sig vs n.s.
isSig_all = (MasterZETA_MUA.ZetaZ_MUA >= U.SigZ) & (MasterZETA_MUA.ZetaP_MUA <= U.SigP); % ZetaZ > 2 and ZetaP < 0.05
y_all = MasterZETA_MUA.ZetaZ_MUA;


% Plot
hAll = figure('Color','w','Position',[80 80 1200 520]); hold on; grid on;

% n.s. gray dots
scatter(xj_all(~isSig_all), y_all(~isSig_all), 28, [0.6 0.6 0.6], ...
    'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName','n.s.');

% significant red stars
scatter(xj_all(isSig_all), y_all(isSig_all), 60, 'r', '*', ...
    'LineWidth', 1.25, 'DisplayName','p<0.05, Z \geq 2');

% threshold
yline(U.SigZ, '--', 'Color',[0.3 0.3 0.3], 'LineWidth', 1, ...
    'DisplayName','ZETA significance threshold');

% x-ticks & labels
xticks(1:numel(uniqSubsAll));
xticklabels(uniqSubsAll);
xtickangle(60);
xlim([0.5, numel(uniqSubsAll)+0.5]);

% y-axis fixed [0, U.YMax]
ylim([0, U.YMax]);
ylabel('ZETA z-score (ZetaZ)');
title(sprintf('ZETA z-scores | All categories (MoveType × Depth)  (N=%d subjects)', numel(uniqSubsAll)));
legend('Location','northeastoutside');

% Save
fnameAll = 'ZETA_Scatter_AllCategories_AllSubjects.png';
print(hAll, fullfile(groupOut, fnameAll), '-dpng', '-r300');
close(hAll);



%% 3b) All-categories-in-one scatter (by category: color=Depth, shape=MoveType, edge=significance)
% update this to be subplots (tilelayout 3 rows, 1 column) for each depth - use color to differentiate movement type in each subplot


% Build x-axis groups from meta-driven labels
prettyPerRow = MasterZETA_MUA.PrettyLabel; % per-row label like "Subject 1 - LSTN"

% x positions (with jitter)
[uniqSubsAll, ~, subjIdxX] = unique(prettyPerRow,'stable'); % unique x-tick strings in encounter order
x_all  = subjIdxX; % x position index for each row
jitter_all = (rand(size(x_all)) - 0.5) * 0.30; % ±0.15 jitter
xj_all = x_all + jitter_all;


% ----- data -----
y_all    = MasterZETA_MUA.ZetaZ_MUA;
isSig    = (MasterZETA_MUA.ZetaZ_MUA >= U.SigZ) & (MasterZETA_MUA.ZetaP_MUA <= U.SigP);
mv_all   = string(MasterZETA_MUA.MoveType);
dz_all   = string(MasterZETA_MUA.Depth);


% ----- MoveType colors (fill) -----
mtColor = containers.Map( ...
    {'HAND OC','HAND PS','ARM EF','REST'}, ...
    {[0.95 0.60 0.10],  ... % Hand OC = orange
    [0.20 0.65 0.30],  ... % Hand PS = green
    [0.15 0.45 0.85],  ... % Arm EF  = blue
    [0.75 0.75 0.75]});    % REST    = light gray

% Edge colors by significance
edgeGray = [0.6 0.6 0.6];
edgeRed  = [0.85 0.10 0.10];

% Depth order & pretty names (row = depth)
depthKeys  = {'t','c','b'};
depthNames = {'dorsal STN','central STN','ventral STN'};

% Figure with 3 rows (one per depth)
hT = figure('Color','w','Position',[90 90 1200 800]);
tlo = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

for dd = 1:numel(depthKeys)
    dzKey = depthKeys{dd};
    ax = nexttile; hold(ax,'on'); grid(ax,'on');

    % within this depth, plot each MoveType with its color and edge by significance
    movesHere = unique(mv_all(dz_all==dzKey),'stable');

    for im = 1:numel(movesHere)
        mv = movesHere(im);
        % pick fill color
        if isKey(mtColor, char(mv)), fc = mtColor(char(mv));
        else, fc = [0.5 0.5 0.5]; end

        idx = (dz_all==dzKey) & (mv_all==mv) & ~isnan(y_all);
        if ~any(idx), continue; end

        % non-sig
        idx_ns = idx & ~isSig;
        if any(idx_ns)
            scatter(ax, xj_all(idx_ns), y_all(idx_ns), 36, fc, ...
                'filled', 'MarkerEdgeColor', edgeGray, 'MarkerFaceAlpha', 0.95, 'MarkerEdgeAlpha', 0.95);
        end

        % sig
        idx_sig = idx & isSig;
        if any(idx_sig)
            scatter(ax, xj_all(idx_sig), y_all(idx_sig), 60, fc, ...
                'filled', 'MarkerEdgeColor', edgeRed, 'LineWidth', 1.2, ...
                'MarkerFaceAlpha', 0.95, 'MarkerEdgeAlpha', 0.95);
        end
    end

    % threshold and axes for this depth
    yline(ax, U.SigZ, '--', 'Color',[0.3 0.3 0.3], 'LineWidth', 1);
    xticks(ax, 1:numel(uniqSubsAll));
    xticklabels(ax, uniqSubsAll);
    xtickangle(ax, 60);
    xlim(ax, [0.5, numel(uniqSubsAll)+0.5]);
    ylim(ax, [0, U.YMax]);
    ylabel(ax, 'ZETA z-score (ZetaZ)');
    title(ax, sprintf('All MoveTypes | %s', depthNames{dd}));
end

% --- Legend (robust across MATLAB versions) ---
% Dummy handles for MoveType colors + edge significance
hOC = scatter(nan,nan,50,mtColor('HAND OC'),'filled','MarkerEdgeColor',edgeGray,'DisplayName','Hand OC');
hPS = scatter(nan,nan,50,mtColor('HAND PS'),'filled','MarkerEdgeColor',edgeGray,'DisplayName','Hand PS');
hEF = scatter(nan,nan,50,mtColor('ARM EF'), 'filled','MarkerEdgeColor',edgeGray,'DisplayName','Arm EF');
hRE = scatter(nan,nan,50,mtColor('REST'),   'filled','MarkerEdgeColor',edgeGray,'DisplayName','Rest');
hNS = scatter(nan,nan,36,[1 1 1],'filled','MarkerEdgeColor',edgeGray,'DisplayName','n.s. edge');
hSG = scatter(nan,nan,36,[1 1 1],'filled','MarkerEdgeColor',edgeRed, 'DisplayName','sig edge');

% Try to dock in the tiledlayout header (newer MATLAB), otherwise overlay legend on an invisible axes
try
    lgd = legend([hOC hPS hEF hRE hNS hSG], 'Orientation','horizontal');
    if isprop(lgd,'Layout') && isprop(lgd.Layout,'Tile')
        lgd.Layout.Tile = 'north';     % docks legend above the tiles
    else
        % Fallback: overlay legend in figure using an invisible, full-figure axes
        delete(lgd); % remove the axes-attached legend first
        axLeg = axes('Parent',hT,'Units','normalized','Position',[0 0 1 1], 'Visible','off');
        lgd = legend(axLeg,[hOC hPS hEF hRE hNS hSG], 'Orientation','horizontal','Box','off');
        % Tweak position near the top center (x y w h in normalized fig coords)
        lgd.Units = 'normalized';
        lgd.Position = [0.32 0.965 0.36 0.03];   % adjust if needed
    end
catch
    % Absolute fallback: put a standard legend below the bottom tile
    lgd = legend([hOC hPS hEF hRE hNS hSG], 'Orientation','horizontal', 'Location','southoutside');
end

% Figure title for the whole tiledlayout
title(tlo, sprintf('ZETA z-scores | All categories by depth '));

fnameAll2 = 'ZETA_Scatter_AllCategories_ByDepth_Tiles.png';
print(hT, fullfile(groupOut, fnameAll2), '-dpng', '-r300');
close(hT);


%% 4) Scatter per category (MoveType × Depth) on x, ZetaZ on y with jitter

cats = unique(MasterZETA_MUA(:, {'MoveType','Depth'}), 'rows', 'stable');
for i = 1:height(cats)
    mv = cats.MoveType(i);
    dz = cats.Depth(i);

    sel = MasterZETA_MUA(MasterZETA_MUA.MoveType==mv & MasterZETA_MUA.Depth==dz, :);
    if isempty(sel), continue; end

    % ---- build friendly labels (safe string construction) ----
    % Subject index in encounter order (Subject 1, Subject 2, ...)
    [~, ~, subjIdxAll] = unique(sel.CaseDate,'stable');
    subjNum = subjIdxAll;  % numeric 1..N

    % % ---- labels by SubjectNum (+ hem) for the selected subset ----
    % hasHemi  = ~(sel.Hemisphere=="" | ismissing(sel.Hemisphere));
    % hemiText = strings(height(sel),1); % initialize all empty
    % hemiText(hasHemi) = "–" + string(sel.Hemisphere(hasHemi));
    % % % don't display hemisphere for now
    % % hemiText(:) = "";

    prettyPerRow = sel.PrettyLabel;
    [uniqSubs, ~, xIdx] = unique(prettyPerRow,'stable');
    x = xIdx;  % (no jitter for the per-category plots)

    % sig vs nonsig
    isSig = (sel.ZetaZ_MUA >= U.SigZ) & (sel.ZetaP_MUA <= U.SigP); % ZetaZ > 2, ZetaP < 0.05
    y = sel.ZetaZ_MUA;

    % Pretty labels
    depthLbl = mapOrDefault(U.DepthMap, dz, dz);
    mvPretty = mapOrDefault(U.PrettyMoveMap, mv, mv);

    % Plot
    h = figure('Color','w','Position',[100 100 1100 500]);
    hold on; grid on;

    % nonsig: gray dots
    scatter(x(~isSig), y(~isSig), 28, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName','n.s.');
    % sig: red stars
    scatter(x(isSig),  y(isSig),  60, 'r', '*', 'LineWidth', 1.25, 'DisplayName','p<0.05, Z \geq 2');

    % dashed threshold line at Z = SigZ
    yline(U.SigZ, '--', 'Color',[0.3 0.3 0.3], 'LineWidth', 1, 'DisplayName','ZETA significance threshold');

    % x-axis ticks per subject, rotate labels
    xticks(1:numel(uniqSubs));
    xticklabels(uniqSubs);
    xtickangle(60);
    xlim([0.5, numel(uniqSubs)+0.5]);

    % y-axis limits (start at 0 for clarity)
    ylim([0, U.YMax]);              % fixed y-axis

    ylabel('ZETA z-score (ZetaZ)');
    title(sprintf('ZETA z-scores | %s × %s  (N=%d subjects)', mvPretty, depthLbl, numel(uniqSubs)));

    legend('Location','northeastoutside');

    % Save
    fname = sprintf('ZETA_Scatter_%s_%s.png', sanitize_filename(mv), sanitize_filename(dz));
    print(h, fullfile(groupOut, fname), '-dpng', '-r300');
    close(h);
end



%% 5) Save MasterZETA_MUA

writetable(MasterZETA_MUA, fullfile(groupOut, 'MasterZETA_MUA_AllSubjects.csv'));

save(fullfile(groupOut, 'MasterZETA_MUA_AllData.mat'), 'MasterZETA_MUA', '-v7.3');

fprintf('[OK] Aggregated %d MUA ZETA rows across %d subjects. Outputs saved to:\n  %s\n', ...
    height(MasterZETA_MUA), numel(unique(MasterZETA_MUA.CaseDate)), groupOut);

end


%% ------- Helper: find MUA ZETA CSVs -------

function list = find_ZetaSummary_MUA_files(baseFRKinDir, zetaCsvPattern)
% Same logic as find_ZetaSummary_files, but looks in "Zeta Testing MUA"

list = struct('fullpath',{},'subject',{},'hemi',{});

caseDirs = dir(fullfile(baseFRKinDir,'*'));
caseDirs = caseDirs([caseDirs.isdir]);
caseDirs = caseDirs(~ismember({caseDirs.name},{'.','..'}));

for i = 1:numel(caseDirs)
    caseName = caseDirs(i).name;
    casePath = fullfile(baseFRKinDir, caseName);

    % Unilateral: Case/Zeta Testing MUA/*.csv
    zetaDir = fullfile(casePath,'Zeta Testing MUA');
    files = dir(fullfile(zetaDir, zetaCsvPattern));

    if isempty(files)
        % Bilateral: Case/LSTN/Zeta Testing MUA/*.csv etc.
        hemiDirs = dir(fullfile(casePath,'*STN'));
        hemiDirs = hemiDirs([hemiDirs.isdir]);
        for h = 1:numel(hemiDirs)
            hemiName = hemiDirs(h).name;
            zetaDirH = fullfile(casePath, hemiName, 'Zeta Testing MUA');
            filesH = dir(fullfile(zetaDirH, zetaCsvPattern));
            for k = 1:numel(filesH)
                list(end+1) = struct( ...
                    'fullpath', fullfile(filesH(k).folder, filesH(k).name), ...
                    'subject', caseName, ...
                    'hemi', hemiName);
            end
        end
    else
        for k = 1:numel(files)
            list(end+1) = struct( ...
                'fullpath', fullfile(files(k).folder, files(k).name), ...
                'subject', caseName, ...
                'hemi', '');
        end
    end
end
end


function val = mapOrDefault(mp, key, defaultVal)
% mapOrDefault: safe lookup for containers.Map with a fallback
%
% mp         : containers.Map
% key        : string/char key to look up
% defaultVal : value to return if key not found or mp not a Map

val = defaultVal;

if ~isa(mp,'containers.Map')
    return;
end

% Normalize key to char (your maps use char keys like 't','HAND OC', etc.)
if isstring(key)
    k = char(key);
else
    k = key;
end

try
    if isKey(mp, k)
        val = mp(k);
    end
catch
    % In case of weird key types, just fall back to defaultVal
    val = defaultVal;
end
end

function s = sanitize_filename(s)
% sanitize_filename: safe filename generator
% Converts MoveType/Depth strings into valid filesystem names.
%
% Examples:
%   "HAND OC" -> "HAND_OC"
%   "ARM EF"  -> "ARM_EF"
%   "b"       -> "b"
%
% Removes characters not allowed in filenames and replaces spaces.

if isstring(s)
    s = char(s);
end

% replace spaces with underscores
s = strrep(s, ' ', '_');

% remove problematic characters
badChars = ['\/:*?"<>|']; % Windows illegal characters
for c = badChars
    s = strrep(s, c, '');
end

% Optionally trim leading/trailing underscores
s = regexprep(s, '^_+', '');
s = regexprep(s, '_+$', '');
end
