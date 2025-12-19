function hFig = plot_IFR_summary_grid_v2(IFR_PSTH_Summary, all_IFR, spikeFieldKey, depthOrder, moveOrder, varargin)
% plot_IFR_summary_grid
% Rows = depths (each depth has 2 rows: raster then IFR/PSTH)
% Cols = MoveTypes
% Saves tiled summary figure optionally.

% -----------------------
% Defaults / ordering
% -----------------------
spikeFieldKey = char(string(spikeFieldKey));

if nargin < 4 || isempty(depthOrder)
    depthOrder = {'t','c','b'};
end
if nargin < 5 || isempty(moveOrder)
    moveOrder  = {'HAND OC','HAND PS','ARM EF'};
end

% -----------------------
% Optional inputs (fix: assign varargin first)
% -----------------------
p = inputParser;
p.addParameter('SaveDir', '', @(s) ischar(s) || isstring(s));
p.addParameter('CaseDate','', @(s) ischar(s) || isstring(s));
p.addParameter('XLim', [-0.05 1.60], @(x) isnumeric(x) && numel(x)==2);
p.parse(varargin{:});
U = p.Results;

% -----------------------
% JNE color scheme
% -----------------------

% Depth colors (purple shades)
purpleShades = ([ ...
    118,42,131;      % dorsal STN (t)
    175,141,195;     % central STN (c)
    231-15,212-15,232-15] ... % ventral STN (b)
    ./ 255); % /255 = standard

% Movement-context colors (greens)
greenShades = ([ ...
    217,240,211;     % HAND OC  (light green)
    127,191,123;     % HAND PS  (green)
    27,120,55] ...   % ARM EF   (dark green)
    ./ 255);  % /255 = standard

depthColor = containers.Map({'t','c','b'}, {purpleShades(1,:), purpleShades(2,:), purpleShades(3,:)});
moveColor  = containers.Map({'HAND OC','HAND PS','ARM EF'}, {greenShades(1,:), greenShades(2,:), greenShades(3,:)});

depthName = containers.Map({'t','c','b'}, {'dorsal','central','ventral'});

nDepths = numel(depthOrder);
nMoveTs = numel(moveOrder);

% -----------------------
% Build key -> index map from all_IFR
% -----------------------
key2idx = containers.Map;
for i = 1:numel(all_IFR)
    s = all_IFR{i};
    sf = char(string(s.SpikeField));
    mv = char(string(s.MoveType));
    dz = char(string(s.Depth));
    key = sprintf('%s|%s|%s', sf, mv, dz);
    if ~isKey(key2idx, key)
        key2idx(key) = i;
    end
end

% -----------------------
% Figure layout
% -----------------------
hFig = figure('Color','w','Position',[50 50 1500 900]);
tlo  = tiledlayout(2*nDepths, nMoveTs, 'TileSpacing','compact','Padding','compact');

for di = 1:nDepths
    dz = char(string(depthOrder{di}));
    cDepth = [0 0 0];
    if isKey(depthColor, dz), cDepth = depthColor(dz); end

    dzLbl = dz;
    if isKey(depthName, dz), dzLbl = depthName(dz); end

    for mi = 1:nMoveTs
        mv = char(string(moveOrder{mi}));
        cMove = [0 0 0];
        if isKey(moveColor, mv), cMove = moveColor(mv); end

        key = sprintf('%s|%s|%s', spikeFieldKey, mv, dz);

        tileRaster = (2*(di-1))*nMoveTs + mi;       % rows 1,3,5...
        tileCurve  = (2*(di-1)+1)*nMoveTs + mi;     % rows 2,4,6...

        % ----- Raster tile -----
        axR = nexttile(tlo, tileRaster);
        cla(axR); hold(axR,'on'); grid(axR,'on');

        if isKey(key2idx, key)
            sI = all_IFR{key2idx(key)};

            if ~isempty(sI.raster_t)
                scatter(axR, sI.raster_t, sI.raster_trial, 6, cDepth, 'filled');
            end
            xline(axR, 0, '-', 'Color',[0.3 0.3 0.3]);

            % FIXED x-axis
            xlim(axR, U.XLim);

            % y-axis = repetition count (total trials/events)
            if isfield(sI,'nReps_behavioral') && ~isempty(sI.nReps_behavioral)
                nRep = double(sI.nReps_behavioral);   % rep-centric y-axis
            else
                nRep = double(sI.nTrials);            % fallback
            end


            ylim(axR, [0.5 nRep+0.5]);

            % tick spacing
            if nRep >= 25 % or 30
                yTickStep = 5;
            elseif nRep >= 10 % or 12
                yTickStep = 2;
            else
                yTickStep = 1;
            end
            yticks(axR, 1:yTickStep:nRep);
            % ylabel(axR, 'Repetition');

            % Decide raster y-axis label (SU vs MUA)
            isMUA = contains(upper(string(spikeFieldKey)), "MUA");
            if isMUA
                yLbl = 'Trial index';     % or: 'Event #', 'Window #', 'Trial index'
            else
                yLbl = 'Repetition';      % SU
            end

            % Only show the raster y-label on the first column
            if mi == 1
                ylabel(axR, yLbl);
            else
                ylabel(axR, '');
            end

        else
            axis(axR,'off');
        end

        if di == 1
            title(axR, mv, 'Interpreter','none');
        end
        if mi == 1
            % depth label once per depth block
            text(axR, U.XLim(1)-0.02*(U.XLim(2)-U.XLim(1)), 0.5, dzLbl, ...
                'Rotation', 90, 'HorizontalAlignment','right', 'VerticalAlignment','bottom', ...
                'FontWeight','bold');
        end

        % ----- IFR + PSTH tile -----
        axC = nexttile(tlo, tileCurve);
        cla(axC); hold(axC,'on'); grid(axC,'on');

        if isKey(key2idx, key)
            sI = all_IFR{key2idx(key)};

            % PSTH dashed (move color)
            if isfield(sI,'centers') && ~isempty(sI.centers) && isfield(sI,'psth_Hz') && ~isempty(sI.psth_Hz)
                plot(axC, sI.centers, sI.psth_Hz, '--', 'LineWidth', 1.6, 'Color', cMove);
            end

            % IFR solid (depth color)
            if ~isempty(sI.vecTime) && ~isempty(sI.vecRate)
                t_ifr = sI.vecTime(:);
                r_ifr = sI.vecRate(:);

                % If vecTime is stitched, your per-category plots already used modulo.
                % For the grid, keep it simple: plot what falls in the fixed peri window.
                mask  = (t_ifr >= U.XLim(1)) & (t_ifr <= U.XLim(2));
                plot(axC, t_ifr(mask), r_ifr(mask), 'LineWidth', 1.8, 'Color', cDepth);
            end

            xline(axC, 0, '-', 'Color',[0.3 0.3 0.3]);
            xlim(axC, U.XLim);
            yticks(axC, []);

        else
            axis(axC,'off');
        end

        if di == nDepths
            xlabel(axC,'Time from onset (s)');
        end
    end
end

title(tlo, sprintf('SpikeField %s: Raster + IFR/PSTH by Depth × MoveType', spikeFieldKey), ...
    'Interpreter','none');

% -----------------------
% Save alongside per-category IFR plots
% -----------------------
if ~isempty(U.SaveDir)
    if ~exist(U.SaveDir,'dir'), mkdir(U.SaveDir); end

    if ~isempty(U.CaseDate)
        base = sprintf('%s_IFRGrid_%s', char(string(U.CaseDate)), spikeFieldKey);
    else
        base = sprintf('IFRGrid_%s', spikeFieldKey);
    end
    base = strrep(base,' ','_');

    print(hFig, fullfile(U.SaveDir, [base,'.png']), '-dpng','-r300');
    savefig(hFig, fullfile(U.SaveDir, [base,'.fig']));
end

end
