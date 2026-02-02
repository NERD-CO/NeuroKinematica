function plot_ZETA_PCA_3D_byDepth_MUA(MasterZETA_MUA, depthCode, varargin)

% plot_ZETA_PCA_3D_byDepth_MUA
%
% 3-D PCA scatterplots (PC1, PC2, PC3) of MUA ZETA temporal deviation
% vector scores (vecD_MUA) for each STN depth, with colored points per
% MoveType and optional trajectory lines.
%
% Modes:
%   depthCode = 't','c','b'      -> single-depth figure
%                 + optional 3-depth tiled figure if DoTiled = true
%   depthCode = 'all','all depths','all_depths'
%                               -> 3x1 tiled figure for all depths
%
% Each point = one MUA unit (subject × channel) at that depth.
% Coordinates = PC1, PC2, PC3 scores.
% Row = STN depth ('t' = dorsal, 'c' = central, 'b' = ventral).
% Color = MoveType.
%
% INPUTS
%   MasterZETA_MUA : table with MUA ZETA data, containing at least:
%       - MoveType
%       - Depth
%       - vecD_MUA
%       - new_vecTime_MUA
%   depthCode      : 't','c','b' OR 'all','all depths','all_depths'
%
% Optional Name–Value:
%   'SavePath'  – directory to save PNGs (default = no save)
%   'DoLines'   – true/false to draw trajectory lines (default = true)
%   'DoTiled'   – true/false to make 3-depth tiled figure (default = false)

%% ---- Parse inputs ----
p = inputParser;
p.addParameter('SavePath', '', @(x) ischar(x) || isstring(x));
p.addParameter('DoLines', true, @islogical);
p.addParameter('DoTiled', false, @islogical);
p.parse(varargin{:});
U = p.Results;

%% --- Sanity check required columns ---
reqVars = {'MoveType','Depth','vecD_MUA','new_vecTime_MUA'};
missing = setdiff(reqVars, MasterZETA_MUA.Properties.VariableNames);
if ~isempty(missing)
    error('MasterZETA_MUA is missing required variables: %s', strjoin(missing,', '));
end

%% --- MoveType + color mapping ---

moveTypes = {'HAND OC','HAND PS','ARM EF','REST'};  % consistent order

% mtColors = containers.Map( ...
%     {'HAND OC','HAND PS','ARM EF','REST'}, ...
%     {[0.95 0.60 0.10], ... % Hand OC = orange
%      [0.20 0.65 0.30], ... % Hand PS = green
%      [0.15 0.45 0.85], ... % Arm EF  = blue
%      [0.60 0.60 0.60]});   % Rest    = gray

% -----------------------
% JNE color scheme (canonical)
% -----------------------
purpleShades = ([ ...
    118,42,131;      % dorsal STN (t)
    175,141,195;     % central STN (c)
    231-15, 212-15, 232-15] ... % ventral STN (b)
    ./ 255);

% greenShades = ([ ...
%     128,128,128;     % REST
%     217-15,240-15,211-15;     % HAND OC
%     127,191,123;     % HAND PS
%     27,120,55] ...   % ARM EF
%     ./ 255);
JNE_move = ([ ...
    128,128,128;    % REST  (grey)
    38,116,183;     % HAND OC  (blue)
    53,183,121;     % HAND PS  (green/teal)
    243,120,98] ...   % ARM EF   (coral)
    ./ 255);  % /255 = standard

depthColorMap = containers.Map({'t','c','b'}, {purpleShades(1,:), purpleShades(2,:), purpleShades(3,:)});
moveColorMap  = containers.Map({'REST','HAND OC','HAND PS','ARM EF'}, {JNE_move(1,:), JNE_move(2,:), JNE_move(3,:), JNE_move(4,:)});
fallbackCol   = [0.5 0.5 0.5];

% Define Actve movement list/color map for part 3 and 4
Active_moveTypes = {'HAND OC','HAND PS','ARM EF'};
Active_moveColorMap = containers.Map( ...
    Active_moveTypes, ...
    {moveColorMap('HAND OC'), moveColorMap('HAND PS'), moveColorMap('ARM EF')} );

% Depth labels
depthNames = containers.Map({'t','c','b'}, ...
    {'dorsal STN','central STN','ventral STN'});


%% --- Interpret depthCode / mode ---
depthCodeStr = string(depthCode);
isAllDepths  = any(strcmpi(depthCodeStr, ["all","all depths","all_depths"]));

if ~isAllDepths
    depthCodeSingle = char(depthCodeStr(1));
    if isKey(depthNames, depthCodeSingle)
        depthLabel = depthNames(depthCodeSingle);
    else
        depthLabel = sprintf('depth %s', depthCodeSingle);
    end
else
    depthCodeSingle = '';
    depthLabel      = '';
end

% common masks
hasZ = ~cellfun(@isempty, MasterZETA_MUA.vecD_MUA);
hasT = ~cellfun(@isempty, MasterZETA_MUA.new_vecTime_MUA);

%% =========================================================
%  PART 1: Single-depth 3D PCA (PC1,PC2,PC3) for depthCode
% ==========================================================
if ~isAllDepths
    isDepth = MasterZETA_MUA.Depth == string(depthCodeSingle);
    M = MasterZETA_MUA(isDepth & hasZ & hasT, :);
    if isempty(M)
        error('No MUA rows with vecD_MUA + new_vecTime_MUA at depth "%s".', depthCodeSingle);
    end

    fprintf('[3D PCA MUA] Depth %s | usable rows: %d\n', depthCodeSingle, height(M));

    % Build time-aligned matrix using helper
    [X, t_use, mtPerUnit, labelsUnit] = build_MUA_ZETA_timeAlignedMatrix(M, ...
        'DepthCode', depthCodeSingle); 

    if size(X,2) < 3
        error('Need at least 3 MUA units at this depth to plot PC1–PC3 (found %d).', size(X,2));
    end

    % --- PCA: observations = units, variables = timepoints
    [coeff, score, latent] = pca(X', 'Centered', true);

    if size(score,2) < 3
        warning('Only %d PCs available; padding PC3 with zeros.', size(score,2));
        score(:,3) = 0;
    end

    % Optional flip of PC1 sign (post-onset > 0)
    pc1 = coeff(:,1);
    postMask = t_use >= 0 & t_use <= (max(t_use)*0.5);
    if any(postMask) && mean(pc1(postMask)) < 0
        coeff(:,1) = -pc1;
        score(:,1) = -score(:,1);
    end

    PC1 = score(:,1);
    PC2 = score(:,2);
    PC3 = score(:,3);

    %% --- Single-depth figure ---
    hFig = figure('Color','w', ...
        'Units','pixels', ...
        'Position',[100 100 900 700]);

    ax = axes('Parent', hFig);
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    ax.Position = [0.13 0.15 0.70 0.72];

    mtOrder  = moveTypes;
    hLegend  = gobjects(1,numel(mtOrder));

    for mIdx = 1:numel(mtOrder)
        mv = mtOrder{mIdx};
        idx = (mtPerUnit == mv);
        if ~any(idx), continue; end

        % if isKey(mtColors,mv), c = mtColors(mv); else, c = [0.5 0.5 0.5]; end
        if isKey(moveColorMap, mv), c = moveColorMap(mv); else, c = fallbackCol; end


        hLegend(mIdx) = scatter3(ax, PC1(idx), PC2(idx), PC3(idx), ...
            50, c, 'filled', 'MarkerFaceAlpha',0.9);

        if U.DoLines && ~strcmp(mv,'REST')
            idxMT = find(idx);
            [~,ord] = sort(PC1(idxMT));
            idxOrd  = idxMT(ord);
            plot3(ax, PC1(idxOrd), PC2(idxOrd), PC3(idxOrd), ...
                'Color', c, 'LineWidth', 1.5, 'HandleVisibility','off');
        end
    end

    xlabel(ax,'PC1 score','FontSize',14,'FontWeight','bold');
    ylabel(ax,'PC2 score','FontSize',14,'FontWeight','bold');
    zlabel(ax,'PC3 score','FontSize',14,'FontWeight','bold');

    title(ax, sprintf('3-D PCA of MUA ZETA deviation | %s', depthLabel), ...
        'FontSize',16);

    legend(ax, hLegend, mtOrder, 'Location','northeastoutside', 'FontSize',11);

    view(ax, [-40 20]);
    axis(ax,'vis3d');

    margin = 0.1;
    xlim(ax, [min(PC1)-margin, max(PC1)+margin]);
    ylim(ax, [min(PC2)-margin, max(PC2)+margin]);
    zlim(ax, [min(PC3)-margin, max(PC3)+margin]);
    ax.LooseInset = max(ax.TightInset, 0.05);

    fprintf('[3D PCA MUA] Depth %s: %d units, %d time points | PCs var %%: %.1f / %.1f / %.1f\n', ...
        depthCodeSingle, size(X,2), size(X,1), ...
        100*latent(1)/sum(latent), ...
        100*latent(2)/sum(latent), ...
        100*latent(3)/sum(latent));

    if ~isempty(U.SavePath)
        if ~exist(U.SavePath,'dir'), mkdir(U.SavePath); end
        outSingle = fullfile(U.SavePath, sprintf('ZETA_PCA_3D_MUA_depth_%s.png', char(depthCodeSingle)));
        exportgraphics(hFig, outSingle, 'Resolution',300);
        fprintf('Saved single-depth MUA 3D PCA figure to:\n  %s\n', outSingle);
    end
end

%% =========================================================
%  PART 2: Tiled layout (1×3) across depths (t,c,b)
% ==========================================================
if U.DoTiled || isAllDepths
    depthList = {'t','c','b'};

    hTile = figure('Color','w','Units','pixels','Position',[150 40 950 1100]);
    tlo   = tiledlayout(hTile, 3, 1, 'TileSpacing','compact','Padding','tight');

    mtOrder   = moveTypes;
    axHandles = gobjects(numel(depthList),1);

    % global limits
    globalXmin = Inf;  globalXmax = -Inf;
    globalYmin = Inf;  globalYmax = -Inf;
    globalZmin = Inf;  globalZmax = -Inf;

    for dIdx = 1:numel(depthList)
        dz = depthList{dIdx};
        ax = nexttile(tlo);
        axHandles(dIdx) = ax;
        hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        ax.FontSize = 11;

        isDepth_d = MasterZETA_MUA.Depth == string(dz);
        Md        = MasterZETA_MUA(isDepth_d & hasZ & hasT, :);
        if isempty(Md)
            title(ax, sprintf('No data | depth %s', dz));
            continue;
        end

        [Xd, t_used, mtPerUnit_d, ~] = build_MUA_ZETA_timeAlignedMatrix(Md, ...
            'DepthCode', dz);

        if size(Xd,2) < 3
            title(ax, sprintf('Too few units for 3D PCA | depth %s', dz));
            continue;
        end

        [coeff_d, score_d, latent_d] = pca(Xd', 'Centered', true);
        if size(score_d,2) < 3
            score_d(:,3) = 0;
        end

        pc1_d      = coeff_d(:,1);
        postMask_d = t_used >= 0 & t_used <= (max(t_used)*0.5);
        if any(postMask_d) && mean(pc1_d(postMask_d)) < 0
            coeff_d(:,1) = -pc1_d;
            score_d(:,1) = -score_d(:,1);
        end

        PC1d = score_d(:,1);
        PC2d = score_d(:,2);
        PC3d = score_d(:,3);

        globalXmin = min(globalXmin, min(PC1d));
        globalXmax = max(globalXmax, max(PC1d));
        globalYmin = min(globalYmin, min(PC2d));
        globalYmax = max(globalYmax, max(PC2d));
        globalZmin = min(globalZmin, min(PC3d));
        globalZmax = max(globalZmax, max(PC3d));

        for mIdx = 1:numel(mtOrder)
            mv   = mtOrder{mIdx};
            mask = (mtPerUnit_d == mv);
            if ~any(mask), continue; end

            % c = mtColors(mv);
            if isKey(moveColorMap, mv), c = moveColorMap(mv); else, c = fallbackCol; end
            scatter3(ax, PC1d(mask), PC2d(mask), PC3d(mask), ...
                30, c, 'filled', 'MarkerFaceAlpha', 0.85);

            if U.DoLines && ~strcmp(mv,'REST')
                idx    = find(mask);
                [~,ord] = sort(PC1d(idx));
                idxOrd  = idx(ord);
                plot3(ax, PC1d(idxOrd), PC2d(idxOrd), PC3d(idxOrd), ...
                    'Color', c, 'LineWidth', 1.3);
            end
        end

        if isKey(depthNames, dz)
            ttl = depthNames(dz);
        else
            ttl = sprintf('depth %s', dz);
        end
        title(ax, ttl, 'FontSize', 14,'FontWeight','bold');

        xlabel(ax,'PC1','FontSize',12);
        ylabel(ax,'PC2','FontSize',12);
        zlabel(ax,'PC3','FontSize',12);
        view(ax, [-40 20]);
        axis(ax,'vis3d');
    end

    % uniform limits
    validAxes = axHandles(isgraphics(axHandles));
    if ~isempty(validAxes) && isfinite(globalXmin)
        mx = 0.1 * max(1, globalXmax - globalXmin);
        my = 0.1 * max(1, globalYmax - globalYmin);
        mz = 0.1 * max(1, globalZmax - globalZmin);

        for ax = reshape(validAxes,1,[])
            xlim(ax, [globalXmin-mx, globalXmax+mx]);
            ylim(ax, [globalYmin-my, globalYmax+my]);
            zlim(ax, [globalZmin-mz, globalZmax+mz]);
        end
    end

    % global legend with dummy handles
    firstAx = validAxes(1);
    hold(firstAx,'on');
    hLeg = gobjects(numel(mtOrder),1);
    for mIdx = 1:numel(mtOrder)
        mv = mtOrder{mIdx};
        % c  = mtColors(mv);
        if isKey(moveColorMap, mv), c = moveColorMap(mv); else, c = fallbackCol; end
        hLeg(mIdx) = scatter3(firstAx, NaN,NaN,NaN, 40, c, ...
            'filled', 'MarkerEdgeColor','none');
    end
    lgd = legend(firstAx, hLeg, mtOrder, 'Location','eastoutside');
    lgd.Title.String = 'MoveType';

    % global title as annotation (keeps space above dorsal tile)
    annotation(hTile,'textbox',[0.05 0.96 0.9 0.04], ...
        'String','3-D PCA of MUA ZETA deviation across STN depths', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'LineStyle','none', ...
        'FontSize',16, ...
        'FontWeight','bold');

    % Save Figure
    if ~isempty(U.SavePath)
        outTile = fullfile(U.SavePath, 'ZETA_PCA_3D_MUA_allDepths.png');
        exportgraphics(hTile, outTile, 'Resolution',300);
        fprintf('Saved tiled MUA 3-depth 3D PCA figure to:\n  %s\n', outTile);
    end
end


%% =========================================================
%  PART 3: Exclude REST - Single-depth 3D PCA (PC1,PC2,PC3)
% ==========================================================
if ~isAllDepths
    isDepth = MasterZETA_MUA.Depth == string(depthCodeSingle);
    M = MasterZETA_MUA(isDepth & hasZ & hasT, :);

    % exclude REST
    M = M(~strcmpi(string(M.MoveType), "REST"), :);

    if isempty(M)
        warning('No non-REST MUA rows at depth "%s". Skipping Part 3.', depthCodeSingle);
    else
        fprintf('[3D PCA MUA NO-REST] Depth %s | usable rows: %d\n', depthCodeSingle, height(M));

        [X, t_use, mtPerUnit, ~] = build_MUA_ZETA_timeAlignedMatrix(M, 'DepthCode', depthCodeSingle);

        if size(X,2) < 3
            warning('Too few MUA units for 3D PCA (No-REST) at depth %s.', depthCodeSingle);
        else
            [coeff, score, latent] = pca(X', 'Centered', true);
            if size(score,2) < 3, score(:,3) = 0; end

            pc1 = coeff(:,1);
            postMask = t_use >= 0 & t_use <= (max(t_use)*0.5);
            if any(postMask) && mean(pc1(postMask)) < 0
                score(:,1) = -score(:,1);
            end

            PC1 = score(:,1); PC2 = score(:,2); PC3 = score(:,3);

            hFig = figure('Color','w','Units','pixels','Position',[100 100 900 700]);
            ax = axes('Parent', hFig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
            ax.Position = [0.13 0.15 0.70 0.72];

            mtOrder = Active_moveTypes;
            hLegend = gobjects(1,numel(mtOrder));

            for mIdx = 1:numel(mtOrder)
                mv = mtOrder{mIdx};
                idx = strcmpi(string(mtPerUnit), string(mv));
                if ~any(idx), continue; end
                c = Active_moveColorMap(mv);

                hLegend(mIdx) = scatter3(ax, PC1(idx), PC2(idx), PC3(idx), ...
                    50, c, 'filled', 'MarkerFaceAlpha',0.9);

                if U.DoLines
                    idxMT = find(idx);
                    [~,ord] = sort(PC1(idxMT));
                    idxOrd  = idxMT(ord);
                    plot3(ax, PC1(idxOrd), PC2(idxOrd), PC3(idxOrd), ...
                        'Color', c, 'LineWidth', 1.5, 'HandleVisibility','off');
                end
            end

            xlabel(ax,'PC1 score','FontSize',14,'FontWeight','bold');
            ylabel(ax,'PC2 score','FontSize',14,'FontWeight','bold');
            zlabel(ax,'PC3 score','FontSize',14,'FontWeight','bold');
            title(ax, sprintf('3-D PCA of MUA ZETA deviation (No REST) | %s', depthLabel), 'FontSize',16);

            legend(ax, hLegend, mtOrder, 'Location','northeastoutside', 'FontSize',11);

            view(ax, [-40 20]);
            axis(ax,'vis3d');

            % Match tiled geometry style
            set(ax, 'PlotBoxAspectRatio', [1 1 1]);
            ax.PositionConstraint = 'innerposition';
            daspect(ax, [1 1 1]);

            margin = 0.1;
            xlim(ax, [min(PC1)-margin, max(PC1)+margin]);
            ylim(ax, [min(PC2)-margin, max(PC2)+margin]);
            zlim(ax, [min(PC3)-margin, max(PC3)+margin]);

            fprintf('[3D PCA MUA NO-REST] Depth %s | PCs var %%: %.1f / %.1f / %.1f\n', ...
                depthCodeSingle, 100*latent(1)/sum(latent), 100*latent(2)/sum(latent), 100*latent(3)/sum(latent));

            if ~isempty(U.SavePath)
                if ~exist(U.SavePath,'dir'), mkdir(U.SavePath); end
                outSingle = fullfile(U.SavePath, sprintf('ZETA_PCA_3D_MUA_depth_NoRest_%s.png', char(depthCodeSingle)));
                exportgraphics(hFig, outSingle, 'Resolution',300);
                fprintf('Saved single-depth MUA 3D PCA (No REST) to:\n  %s\n', outSingle);
            end
        end
    end
end


%% =========================================================
%  PART 4: Exclude REST - Tiled layout (1×3) across depths
% ==========================================================
if U.DoTiled || isAllDepths
    depthList = {'t','c','b'};
    hTile = figure('Color','w','Units','pixels','Position',[150 40 950 1100]);
    tlo   = tiledlayout(hTile, 3, 1, 'TileSpacing','compact','Padding','tight');

    mtOrder   = Active_moveTypes;
    axHandles = gobjects(numel(depthList),1);

    globalXmin = Inf; globalXmax = -Inf;
    globalYmin = Inf; globalYmax = -Inf;
    globalZmin = Inf; globalZmax = -Inf;

    for dIdx = 1:numel(depthList)
        dz = depthList{dIdx};
        ax = nexttile(tlo);
        axHandles(dIdx) = ax;
        hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        ax.FontSize = 11;

        isDepth_d = MasterZETA_MUA.Depth == string(dz);
        Md = MasterZETA_MUA(isDepth_d & hasZ & hasT, :);

        % exclude REST
        Md = Md(~strcmpi(string(Md.MoveType), "REST"), :);

        if isempty(Md)
            title(ax, sprintf('No non-REST data | depth %s', dz));
            continue;
        end

        [Xd, t_used, mtPerUnit_d, ~] = build_MUA_ZETA_timeAlignedMatrix(Md, 'DepthCode', dz);
        if size(Xd,2) < 3
            title(ax, sprintf('Too few units (No REST) | depth %s', dz));
            continue;
        end

        [coeff_d, score_d] = pca(Xd', 'Centered', true);
        if size(score_d,2) < 3, score_d(:,3) = 0; end

        pc1_d      = coeff_d(:,1);
        postMask_d = t_used >= 0 & t_used <= (max(t_used)*0.5);
        if any(postMask_d) && mean(pc1_d(postMask_d)) < 0
            score_d(:,1) = -score_d(:,1);
        end

        PC1d = score_d(:,1); PC2d = score_d(:,2); PC3d = score_d(:,3);

        globalXmin = min(globalXmin, min(PC1d));  globalXmax = max(globalXmax, max(PC1d));
        globalYmin = min(globalYmin, min(PC2d));  globalYmax = max(globalYmax, max(PC2d));
        globalZmin = min(globalZmin, min(PC3d));  globalZmax = max(globalZmax, max(PC3d));

        for mIdx = 1:numel(mtOrder)
            mv = mtOrder{mIdx};
            mask = strcmpi(string(mtPerUnit_d), string(mv));
            if ~any(mask), continue; end
            c = Active_moveColorMap(mv);

            scatter3(ax, PC1d(mask), PC2d(mask), PC3d(mask), ...
                30, c, 'filled', 'MarkerFaceAlpha', 0.85);

            if U.DoLines
                idx = find(mask);
                [~,ord] = sort(PC1d(idx));
                idxOrd = idx(ord);
                plot3(ax, PC1d(idxOrd), PC2d(idxOrd), PC3d(idxOrd), ...
                    'Color', c, 'LineWidth', 1.3);
            end
        end

        title(ax, depthNames(dz), 'FontSize', 14,'FontWeight','bold');
        xlabel(ax,'PC1','FontSize',12);
        ylabel(ax,'PC2','FontSize',12);
        zlabel(ax,'PC3','FontSize',12);
        view(ax, [-40 20]);
        axis(ax,'vis3d');

        % match SU/MUA geometry
        set(ax, 'PlotBoxAspectRatio', [1 1 1]);
        ax.PositionConstraint = 'innerposition';
        daspect(ax, [1 1 1]);
    end

    validAxes = axHandles(isgraphics(axHandles));
    if ~isempty(validAxes) && isfinite(globalXmin)
        mx = 0.1 * max(1, globalXmax - globalXmin);
        my = 0.1 * max(1, globalYmax - globalYmin);
        mz = 0.1 * max(1, globalZmax - globalZmin);

        for ax = reshape(validAxes,1,[])
            xlim(ax, [globalXmin-mx, globalXmax+mx]);
            ylim(ax, [globalYmin-my, globalYmax+my]);
            zlim(ax, [globalZmin-mz, globalZmax+mz]);

            % re-assert after limits
            set(ax, 'PlotBoxAspectRatio', [1 1 1]);
            ax.PositionConstraint = 'innerposition';
        end
    end

    % legend
    firstAx = validAxes(1); hold(firstAx,'on');
    hLeg = gobjects(numel(mtOrder),1);
    for mIdx = 1:numel(mtOrder)
        mv = mtOrder{mIdx};
        c  = Active_moveColorMap(mv);
        hLeg(mIdx) = scatter3(firstAx, NaN,NaN,NaN, 40, c, 'filled', 'MarkerEdgeColor','none');
    end
    lgd = legend(firstAx, hLeg, mtOrder, 'Location','eastoutside');
    lgd.Title.String = 'MoveType';

    annotation(hTile,'textbox',[0.05 0.96 0.9 0.04], ...
        'String','3-D PCA of MUA ZETA deviation across STN depths (No REST)', ...
        'HorizontalAlignment','center', 'VerticalAlignment','top', ...
        'LineStyle','none', 'FontSize',16, 'FontWeight','bold');

    if ~isempty(U.SavePath)
        outTile = fullfile(U.SavePath, 'ZETA_PCA_3D_MUA_allDepths_NoRest.png');
        exportgraphics(hTile, outTile, 'Resolution',300);
        fprintf('Saved tiled MUA 3-depth 3D PCA (No REST) to:\n  %s\n', outTile);
    end
end


end
