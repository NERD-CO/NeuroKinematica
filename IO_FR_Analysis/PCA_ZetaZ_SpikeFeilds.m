%% PCA of ZETA z-scores per MoveType × STN Depth
clear; clc;

FR_Kin_Dir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\Ephys_Kinematics\FR_Kinematic_Analyses';
Aggr_ZETA_dir = fullfile(FR_Kin_Dir, 'Aggregate Zeta Plots');
cd(Aggr_ZETA_dir);

% Load master ZETA file (produced by aggregate_ZETA_and_plot)
load('MasterZeta_AllCases.mat', 'MasterZETA');  % adjust name as needed
% Expect fields: MasterZETA.ZetaZ (time × units), MasterZETA.MoveType, MasterZETA.Depth, MasterZETA.UnitID, MasterZETA.Time_ms

%% Extract metadata
ZetaZ     = MasterZETA.ZetaZ;       % [timeBins × units]
MoveTypes = MasterZETA.MoveType;    % 1×units cell array
Depths    = MasterZETA.Depth;       % 1×units cell array
time_ms   = MasterZETA.Time_ms(:);  % ensure column

uniqMoves = unique(MoveTypes, 'stable');
uniqDepths = unique(Depths, 'stable');

%% Parameters
nTime = numel(time_ms);
fprintf('Data shape: %d time bins × %d units\n', size(ZetaZ));

%% Initialize output
PCA_results = struct();

%% Loop by depth and movement
for d = 1:numel(uniqDepths)
    depth = uniqDepths{d};
    figure('Name', ['PCA_', depth], 'Color', 'w');
    t = tiledlayout(ceil(numel(uniqMoves)/3), 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(t, sprintf('ZETA PCA – %s STN', depth));

    for m = 1:numel(uniqMoves)
        mv = uniqMoves{m};
        idx = strcmp(MoveTypes, mv) & strcmp(Depths, depth);
        if ~any(idx)
            continue
        end

        % Build data matrix: rows=time bins, cols=units
        X = ZetaZ(:, idx);

        % Check for valid variance
        if size(X,2) < 3 || all(isnan(X(:)))
            continue
        end

        % Mean-center across time (for each unit)
        X = fillmissing(X, 'linear', 1, 'EndValues', 'nearest'); % handle NaNs
        Xc = X - mean(X,1);    % column-wise centering

        % PCA across units (columns)
        [coeff, score, latent, ~, explained] = pca(Xc','Centered',false);

        % Store
        PCA_results.(depth).(mv).coeff = coeff;
        PCA_results.(depth).(mv).score = score;
        PCA_results.(depth).(mv).latent = latent;
        PCA_results.(depth).(mv).explained = explained;
        PCA_results.(depth).(mv).Xc = Xc;
        PCA_results.(depth).(mv).time_ms = time_ms;

        % Plot first PC time-series (unit projection on PC1)
        nexttile;
        plot(time_ms, coeff(:,1), 'LineWidth', 1.5);
        xlabel('Time (ms)');
        ylabel('PC1 loading');
        title(sprintf('%s', mv), 'Interpreter','none');
        grid on;
    end
end

%% Save results
save(fullfile(Aggr_ZETA_dir,'PCA_ZETA_byMoveDepth.mat'),'PCA_results');
