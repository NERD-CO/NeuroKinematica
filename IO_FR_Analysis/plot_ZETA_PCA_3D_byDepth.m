function plot_ZETA_PCA_3D_byDepth(MasterZETA, depthCode, varargin)

% plot_ZETA_PCA3D_byDepth
%
% Creates 3-D PCA scatterplots (PC1, PC2, PC3) of single-unit ZETA
% temporal deviation vector scores for each STN depth, with colored
% points per MoveType and optional trajectory lines.
%
% Modes:
% - if depthCode = 't','c','b', then single-depth 3D PCA figure
%     + optional 3-depth tiled figure if DoTiled = true
% - if depthCode = 'all','all depths','all_depths', then tiled layout (3×1)
%                                                   figure for all depths
%
% Each point = one unit (subject × cluster) at that depth.
% Coordinates = PC1, PC2, PC3 scores.
% Row = STN depth ('t' = dorsal, 'c' = central, 'b' = ventral).
% Color = MoveType.
%
% INPUTS
%   MasterZETA : table from aggregate_ZETA_and_plot (single-unit)
%   depthCode  : 't','c','b' OR 'all','all depths','all_depths'
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

reqVars = {'MoveType','Depth','ZETA_vecD','ZETA_vecT', ...
    'PSTH_TimeCenters_s','PrettyLabel'};
missing = setdiff(reqVars, MasterZETA.Properties.VariableNames);
if ~isempty(missing)
    error('MasterZETA is missing required variables: %s', strjoin(missing,', '));
end


%% --- MoveType + color mapping ---

moveTypes = {'HAND OC','HAND PS','ARM EF','REST'};  % consistent plotting order

mtColors = containers.Map( ...
    {'HAND OC','HAND PS','ARM EF','REST'}, ...
    {[0.95 0.60 0.10], ... % Hand OC = orange
    [0.20 0.65 0.30], ... % Hand PS = green
    [0.15 0.45 0.85], ... % Arm EF  = blue
    [0.60 0.60 0.60]});   % Rest    = gray

% Depth labels
depthNames = containers.Map({'t','c','b'}, ...
    {'dorsal STN','central STN','ventral STN'});

%% --- Interpret depthCode / mode ---

depthCodeStr = string(depthCode);
isAllDepths  = any(strcmpi(depthCodeStr, ["all","all depths","all_depths"]));

% For single-depth mode, just use first entry
if ~isAllDepths
    depthCodeSingle = char(depthCodeStr(1));
    if isKey(depthNames, depthCodeSingle)
        depthLabel = depthNames(depthCodeSingle);
    else
        depthLabel = sprintf('depth %s', depthCodeSingle);
    end
else
    depthCodeSingle = '';   % not used
    depthLabel      = '';
end

% Common masks for valid ZETA rows
hasZ = ~cellfun(@isempty, MasterZETA.ZETA_vecD);
hasT = ~cellfun(@isempty, MasterZETA.PSTH_TimeCenters_s);


%% =========================================================
%  PART 1: Single-depth 3D PCA (PC1,PC2,PC3) for depthCode
% ==========================================================

isDepth = MasterZETA.Depth == string(depthCode);
hasZ    = ~cellfun(@isempty, MasterZETA.ZETA_vecD);
hasT    = ~cellfun(@isempty, MasterZETA.PSTH_TimeCenters_s);

if ~isAllDepths
    % Subset rows for requested depth
    isDepth = MasterZETA.Depth == string(depthCodeSingle);
    M = MasterZETA(isDepth & hasZ & hasT, :);
    if isempty(M)
        error('No rows with ZETA_vecD + PSTH_TimeCenters_s at depth "%s".', depthCodeSingle);
    end

    fprintf('[3D PCA] Depth %s | usable SU rows: %d\n', depthCodeSingle, height(M));

    % Build time-aligned matrix using helper
    [X, t_use, mtPerUnit, labelsUnit] = build_SU_ZETA_timeAlignedMatrix(M, ...
        'DepthCode', depthCodeSingle);

    if size(X,2) < 3
        error('Need at least 3 units at this depth to plot PC1–PC3 (found %d).', size(X,2));
    end

    %% --- PCA (observations = units, variables = timepoints) ---

    [coeff, score, latent] = pca(X', 'Centered', true); % score: [nUnits × nPC]

    if size(score,2) < 3
        warning('Only %d PCs available; padding PC3 with zeros.', size(score,2));
        % pad to allow plotting without error
        score(:,3) = 0;
    end

    % Optional: flip PC1 sign so mean post-onset deflection is positive
    pc1 = coeff(:,1);
    postMask = t_use >= 0 & t_use <= (max(t_use)*0.5);
    if any(postMask) && mean(pc1(postMask)) < 0
        coeff(:,1)   = -pc1;
        score(:,1)   = -score(:,1);
    end

    PC1 = score(:,1);
    PC2 = score(:,2);
    PC3 = score(:,3);


    %% --- FIGURE: clean, not-enormous 3D scatter ---

    hFig = figure('Color','w', ...
        'Units','pixels', ...
        'Position',[100 100 900 700]);   % controlled size

    ax = axes('Parent', hFig);
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    % Make axes slightly smaller than full figure to avoid clipping labels
    ax.Position = [0.13 0.15 0.70 0.72];   % [left bottom width height]

    % Plot each MoveType with its own color
    mtOrder = moveTypes; % {'HAND OC','HAND PS','ARM EF','REST'};
   
    % Collect handles for legend
    hLegend = gobjects(1,numel(mtOrder));

    for mIdx = 1:numel(mtOrder)
        mv = mtOrder{mIdx};
        idx = mtPerUnit == mv;
        if ~any(idx), continue; end

        if isKey(mtColors, mv)
            c = mtColors(mv);
        else
            c = [0.5 0.5 0.5];
        end

        % % Scatter points
        % scatter3(ax, PC1(idx), PC2(idx), PC3(idx), ...
        %     50, c, 'filled', 'MarkerFaceAlpha',0.9, 'DisplayName',mv);

        % Scatter and capture handle
        hLegend(mIdx) = scatter3(ax, PC1(idx), PC2(idx), PC3(idx), ...
             50, c, 'filled', 'MarkerFaceAlpha',0.9);

        % Optional trajectory lines (exclude REST)
        if U.DoLines && ~strcmp(mv,'REST')
            % indices for this MoveType
            idx = find(mtPerUnit == mv);
            [~,ord] = sort(PC1(idx));
            idxOrd  = idx(ord);
            plot3(ax, PC1(idxOrd), PC2(idxOrd), PC3(idxOrd), ...
                'Color', c, 'LineWidth', 1.5, 'HandleVisibility','off');
        end
    end

    % Axis labels
    xlabel(ax,'PC1 score','FontSize',14,'FontWeight','bold');
    ylabel(ax,'PC2 score','FontSize',14,'FontWeight','bold');
    zlabel(ax,'PC3 score','FontSize',14,'FontWeight','bold');

    % Title
    title(ax, sprintf('3-D PCA of ZETA deviation | %s', depthLabel), ...
        'FontSize',16);

    % Legend
    legend(ax, hLegend, mtOrder, 'Location','northeastoutside', 'FontSize', 11);

    % View / scaling
    view(ax, [-40 20]);   % azimuth, elevation
    axis(ax,'vis3d');     % preserves aspect ratio under rotation

    % Optional: center axes around data with small margin
    margin = 0.1;
    xlim(ax, [min(PC1)-margin, max(PC1)+margin]);
    ylim(ax, [min(PC2)-margin, max(PC2)+margin]);
    zlim(ax, [min(PC3)-margin, max(PC3)+margin]);

    % Small padding so labels never touch figure edges
    ax.LooseInset = max(ax.TightInset, 0.05);

    fprintf('[3D PCA] Depth %s: %d units, %d time points | PCs var %%: %.1f / %.1f / %.1f\n', ...
        depthCodeSingle, size(X,2), size(X,1), ...
        100*latent(1)/sum(latent), ...
        100*latent(2)/sum(latent), ...
        100*latent(3)/sum(latent));

    % Save
    if ~isempty(U.SavePath)
        if ~exist(U.SavePath,'dir'), mkdir(U.SavePath); end
        outSingle = fullfile(U.SavePath, sprintf('ZETA_PCA_3D_depth_%s.png', char(depthCodeSingle)));
        exportgraphics(hFig, outSingle, 'Resolution',300);
        fprintf('Saved single-depth 3D PCA figure to:\n  %s\n', outSingle);
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

    % For global, uniform axes
    globalXmin = Inf;  globalXmax = -Inf;
    globalYmin = Inf;  globalYmax = -Inf;
    globalZmin = Inf;  globalZmax = -Inf;

    for dIdx = 1:numel(depthList)
        dz = depthList{dIdx};
        ax = nexttile(tlo);
        axHandles(dIdx) = ax;
        hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        ax.FontSize = 11;

        % Subset rows for this depth
        isDepth_d = MasterZETA.Depth == string(dz);
        Md        = MasterZETA(isDepth_d & hasZ & hasT, :);
        if isempty(Md)
            title(ax, sprintf('No data | depth %s', dz));
            continue;
        end

        % Build aligned matrix and PCA for this depth
        [Xd, t_used, mtPerUnit_d, ~] = build_SU_ZETA_timeAlignedMatrix(Md, ...
            'DepthCode', dz);

        if size(Xd,2) < 3
            title(ax, sprintf('Too few units for 3D PCA | depth %s', dz));
            continue;
        end

        [coeff_d, score_d, latent_d] = pca(Xd', 'Centered', true);
        if size(score_d,2) < 3
            score_d(:,3) = 0;
        end

        % Flip PC1 sign to positive post-onset
        pc1_d      = coeff_d(:,1);
        postMask_d = t_used >= 0 & t_used <= (max(t_used)*0.5);
        if any(postMask_d) && mean(pc1_d(postMask_d)) < 0
            coeff_d(:,1) = -pc1_d;
            score_d(:,1) = -score_d(:,1);
        end

        PC1d = score_d(:,1);
        PC2d = score_d(:,2);
        PC3d = score_d(:,3);

        % Update global limits
        globalXmin = min(globalXmin, min(PC1d));
        globalXmax = max(globalXmax, max(PC1d));
        globalYmin = min(globalYmin, min(PC2d));
        globalYmax = max(globalYmax, max(PC2d));
        globalZmin = min(globalZmin, min(PC3d));
        globalZmax = max(globalZmax, max(PC3d));

        % Scatter + optional trajectories per MoveType
        for mIdx = 1:numel(mtOrder)
            mv   = mtOrder{mIdx};
            mask = (mtPerUnit_d == mv);
            if ~any(mask), continue; end

            c = mtColors(mv);

            % scatter
            scatter3(ax, PC1d(mask), PC2d(mask), PC3d(mask), ...
                30, c, 'filled', 'MarkerFaceAlpha', 0.85, 'DisplayName', mv);

            % trajectory lines (exclude REST)
            if U.DoLines && ~strcmp(mv,'REST')
                idx    = find(mask);
                [~,ord] = sort(PC1d(idx));
                idxOrd  = idx(ord);
                plot3(ax, PC1d(idxOrd), PC2d(idxOrd), PC3d(idxOrd), ...
                    'Color', c, 'LineWidth', 1.3);
            end
        end

        % Per-depth title: use mapping directly (no "STN STN")
        if isKey(depthNames, dz)
            ttl = depthNames(dz);     % e.g. 'dorsal STN'
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

    % Apply uniform limits across all valid axes
    validAxes = axHandles(isgraphics(axHandles));
    if ~isempty(validAxes) && isfinite(globalXmin)
        % add a small margin around the global range
        mx = 0.1 * max(1, globalXmax - globalXmin);
        my = 0.1 * max(1, globalYmax - globalYmin);
        mz = 0.1 * max(1, globalZmax - globalZmin);

        for ax = reshape(validAxes,1,[])
            xlim(ax, [globalXmin-mx, globalXmax+mx]);
            ylim(ax, [globalYmin-my, globalYmax+my]);
            zlim(ax, [globalZmin-mz, globalZmax+mz]);
        end
    end

    % --- Global legend (dummy handles so all MoveTypes appear) ---
    firstAx = validAxes(1);
    hold(firstAx,'on');

    hLeg = gobjects(numel(mtOrder),1);
    for mIdx = 1:numel(mtOrder)
        mv = mtOrder{mIdx};
        c  = mtColors(mv);
        % dummy invisible scatter for legend only
        hLeg(mIdx) = scatter3(firstAx, NaN, NaN, NaN, ...
            40, c, 'filled', 'MarkerEdgeColor','none');
    end

    lgd = legend(firstAx, hLeg, mtOrder, 'Location','eastoutside');
    lgd.Title.String = 'MoveType';


    % --- Global title as annotation (avoids overlap with dorsal title) ---
    annotation(hTile,'textbox',[0.05 0.96 0.9 0.04], ...
        'String','3-D PCA of ZETA deviation across STN depths', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'LineStyle','none', ...
        'FontSize',16, ...
        'FontWeight','bold');


    % Save tiled figure
    if ~isempty(U.SavePath)
        outTile = fullfile(U.SavePath, 'ZETA_PCA_3D_allDepths.png');
        exportgraphics(hTile, outTile, 'Resolution',300);
        fprintf('Saved tiled 3-depth 3D PCA figure to:\n  %s\n', outTile);
    end
end

end
