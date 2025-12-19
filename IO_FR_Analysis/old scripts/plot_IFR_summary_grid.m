function hFig = plot_IFR_summary_grid(IFR_PSTH_Summary, all_IFR, spikeFieldKey, depthOrder, moveOrder)

% plot_IFR_summary_grid
% Rows = depths (each depth has 2 rows: raster then IFR/PSTH)
% Cols = MoveTypes
% Saves tiled summary figure optionally.


% -----------------------
% Defaults / ordering
% -----------------------
spikeFieldKey = char(string(spikeFieldKey));

if nargin < 4 || isempty(depthOrder)
    depthOrder = {'t','c','b'};                  % dorsal, central, ventral
end
if nargin < 5 || isempty(moveOrder)
    moveOrder  = {'HAND OC','HAND PS','ARM EF'}; % columns
end


% -----------------------
% Optional inputs (fix: assign varargin first)
% -----------------------
args = varargin;   % <- this avoids brace-indexing issues in some patterns
p = inputParser;
p.addParameter('SaveDir', '', @(s) ischar(s) || isstring(s));
p.addParameter('CaseDate','', @(s) ischar(s) || isstring(s));
p.addParameter('XLim', [-0.05 1.60], @(x) isnumeric(x) && numel(x)==2);
p.parse(args{:});
U = p.Results;

% -----------------------
% JNE color scheme (yours)
% -----------------------
purpleShades = [118,42,131;
    175,141,195;
    231,212,232] ./ 256;  % t,c,b

greenShades  = [217,240,211;   % HAND OC
    127,191,123;   % HAND PS
    27,120,55] ./ 256;    % ARM EF

depthColor = containers.Map({'t','c','b'}, {purpleShades(1,:), purpleShades(2,:), purpleShades(3,:)});
moveColor  = containers.Map({'HAND OC','HAND PS','ARM EF'}, {greenShades(1,:), greenShades(2,:), greenShades(3,:)});

depthName = containers.Map({'t','c','b'}, {'dorsal','central','ventral'});
nD = numel(depthOrder);
nM = numel(moveOrder);

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
%   rows: 2 per depth (raster row, IFR/PSTH row)
%   cols: MoveTypes in moveOrder
% -----------------------
hFig = figure('Color','w','Position',[50 50 1500 900]);
tlo  = tiledlayout(2*nD, nM, 'TileSpacing','compact','Padding','compact');

% -----------------------
% Plot each cell using explicit tile indices
% -----------------------
for di = 1:nD
    dz = char(string(depthOrder{di}));
    if isKey(depthColor, dz), cDepth = depthColor(dz); else, cDepth = [0 0 0]; end

    for mi = 1:nM
        mv = char(string(moveOrder{mi}));
        if isKey(moveColor, mv), cMove = moveColor(mv); else, cMove = [0 0 0]; end

        key = sprintf('%s|%s|%s', spikeFieldKey, mv, dz);

        % Compute tile indices so raster is ALWAYS above its IFR/PSTH
        tileRaster = (2*(di-1))*nM + mi;       % row 1,3,5...
        tileCurve  = (2*(di-1)+1)*nM + mi;     % row 2,4,6...

        % ----- Raster tile -----
        axR = nexttile(tlo, tileRaster);
        cla(axR); hold(axR,'on'); grid(axR,'on');

        if isKey(key2idx, key)
            sI = all_IFR{key2idx(key)};

            if ~isempty(sI.raster_t)
                scatter(axR, sI.raster_t, sI.raster_trial, 6, [0 0 0], 'filled');
            end
            xline(axR, 0, '-', 'Color',[0.3 0.3 0.3]);

            xlim(axR, [-0.05, sI.useMaxDur_s]); % matches your prewindow default; adjust if needed
            xlim(axR, U.XLim); %% new

            % ylim(axR, [0.5, max(1, max(sI.raster_trial))+0.5]);
            % yticks(axR, []);

            % Raster y-axis: Show repetition count      %% new
            nReps = max(1, max(sI.raster_trial));
            ylim(axR, [0.5 nReps+0.5]);

            % % If lots of reps, don't label every tick
            % if nReps <= 20
            %     yticks(axR, 1:nReps);
            % else
            %     yticks(axR, unique(round(linspace(1,nReps,6))));
            % end
            %
            % if mi == 1
            %     ylabel(axR, sprintf('%s\nRep #', dzPretty), 'Interpreter','none');
            % else
            %     ylabel(axR, ''); % keep clean
            % end

            % Choose tick spacing adaptively
            if nRep >= 15
                yTickStep = 5;
            elseif nRep >= 6
                yTickStep = 2;
            else
                yTickStep = 1;
            end

            yticks(axR, 1:yTickStep:nRep);
            ylabel(axR, 'Repetition');

        else
            axis(axR,'off');
        end

        if di == 1
            title(axR, mv, 'Interpreter','none');
        end

        depthPretty = containers.Map({'t','c','b'}, {'dorsal','central','ventral'});
        dzPretty = dz;
        if isKey(depthPretty, dz), dzPretty = depthPretty(dz); end

        if mi == 1
            ylabel(axR, dzPretty, 'Interpreter','none');
        end

        % ----- IFR + PSTH tile -----
        axC = nexttile(tlo, tileCurve);
        cla(axC); hold(axC,'on'); grid(axC,'on');

        if isKey(key2idx, key)
            sI = all_IFR{key2idx(key)};

            % PSTH (MoveType color)
            if isfield(sI,'centers') && ~isempty(sI.centers) && isfield(sI,'psth_Hz') && ~isempty(sI.psth_Hz)
                plot(axC, sI.centers, sI.psth_Hz, '--', 'LineWidth', 1.6, 'Color', cMove);
            end

            % IFR (Depth color) — plot peri-stim portion only
            if ~isempty(sI.vecTime) && ~isempty(sI.vecRate)
                t_ifr = sI.vecTime(:);
                r_ifr = sI.vecRate(:);

                % if vecTime is already peri-stim, this is fine; if it's stitched, you may want
                % your modulo logic here instead.
                mask  = (t_ifr >= -0.05) & (t_ifr <= sI.useMaxDur_s);
                plot(axC, t_ifr(mask), r_ifr(mask), 'LineWidth', 1.8, 'Color', cDepth);
            end

            xline(axC, 0, '-', 'Color',[0.3 0.3 0.3]);
            xlim(axC, [-0.05, sI.useMaxDur_s]);
            xlim(axC, U.XLim); %% new
            yticks(axC, []);
        else
            axis(axC,'off');
        end

        if di == nD
            xlabel(axC,'Time from onset (s)');
        end
    end
end

title(tlo, sprintf('SpikeField %s: Raster + IFR/PSTH by Depth × MoveType', spikeFieldKey), ...
    'Interpreter','none');

%% Save tiled layout summary figure

if ~isempty(U.SaveDir)
    if ~exist(U.SaveDir,'dir'); mkdir(U.SaveDir); end

    sfTag = sanitize_filename(spikeFieldKey); % you already have this helper elsewhere; or just replace spaces
    if ~isempty(U.CaseDate)
        base = sprintf('%s_IFR-SUMMARYGRID_%s', char(U.CaseDate), sfTag);
    else
        base = sprintf('IFR-SUMMARYGRID_%s', sfTag);
    end

    print(hFig, fullfile(char(U.SaveDir), [base '.png']), '-dpng','-r300');
    savefig(hFig, fullfile(char(U.SaveDir), [base '.fig']));
end


%%
% Optional: add a legend once (dummy lines):
%{
axL = nexttile(tlo, 1);
hold(axL,'on');
hIFR = plot(axL, nan, nan, '-',  'LineWidth',2, 'Color', purpleShades(1,:));
hPSTH= plot(axL, nan, nan, '--', 'LineWidth',2, 'Color', greenShades(1,:));
legend(axL, [hIFR hPSTH], {'IFR (Depth color)','PSTH (Move color)'}, 'Location','best');
%}

end
