function row_indices = plot_LFP_byDepthTrial(LFPs_tbl, STN_depth, trialNum, repIdx, fs, lfpCol)

% plotting funtion for LFP segments by STN depth and trial

% LFPs_tbl : All_LFPsPerMove_Tbl_filt
% STN_depth: 'dorsal'|'central'|'ventral'
% trialNum : e.g., 3 or '3'
% repIdx   : numeric index (e.g., 1) OR 'all'
% fs       : sampling rate of processed LFP (e.g., 500)
% lfpCol   : optional (e.g., 'LFP_E1_proc500'); if omitted, picks first *_proc500

% Returns: rows_used (indices in tbl used for plotting)

    if nargin < 6 || isempty(lfpCol)
        var_names = string(LFPs_tbl.Properties.VariableNames);
        cols = var_names(startsWith(var_names,"LFP_E") & endsWith(var_names,"_proc500"));
        if isempty(cols), error('No processed LFP columns (LFP_E*_proc500) found.'); end
        lfpCol = char(cols(1));
    end
    if nargin < 5 || isempty(fs), fs = 500; end

    % Map depth word -> letter
    switch lower(STN_depth)
        case 'dorsal',  depth_prefix = 't';
        case 'central', depth_prefix = 'c';
        case 'ventral', depth_prefix = 'b';
        otherwise, error('Depth must be dorsal|central|ventral');
    end

    % Trial string + simple safe pattern: start with depth+trial (e.g., 't3') and then boundary/underscore/end
    trial_str = string(trialNum);
    pat  = sprintf('^%s%s(\\b|_|$)', depth_prefix, trial_str);
    ids  = string(LFPs_tbl.move_trial_ID);
    rows = find(~cellfun('isempty', regexp(ids, pat, 'once')));

    if isempty(rows)
        error('No rows found for depth "%s" (letter %s) and trial "%s".', STN_depth, depth_prefix, trial_str);
    end
    row_indices = rows;

    % Which reps to plot?
    if ischar(repIdx) || isstring(repIdx)
        if ~strcmpi(string(repIdx), 'all'), error('repIdx must be numeric or ''all''.'); end
        repList = 1:numel(rows);
    else
        if repIdx < 1 || repIdx > numel(rows)
            warning('repIdx=%d is out of range (1..%d). Using 1.', repIdx, numel(rows));
            repIdx = 1;
        end
        repList = repIdx;
    end

    % Plot each requested rep (time-series + PSD)
    for rep_i = repList
        LFP_seg = LFPs_tbl.(lfpCol){rows(rep_i)};
        t_fs   = (0:numel(LFP_seg)-1)/fs;

        plot_title = sprintf('%s STN | Trial %s | Rep %d', ...
            STN_depth, trial_str, rep_i);

        % % Time series
        % figure('Name', plot_title, 'Color','w');
        % plot(t_fs, LFP_seg, 'LineWidth', 1);
        % xlabel('Time (s)'); ylabel('LFP amplitude (\muV)'); grid on;
        % title(plot_title);

        % PSD 0–50 Hz (power in dB)
        [Pxx, Fxx] = pspectrum(LFP_seg, fs, 'FrequencyLimits',[0 50], 'FrequencyResolution', 3);
        Pxx_db = pow2db(Pxx);
        figure('Name', [sprintf('%s STN', STN_depth) ' PSD']);
        plot(Fxx, Pxx_db, 'LineWidth', 1);
        xlim([0 50]); xlabel('Frequency (Hz)'); 
        ylim([-20 40]); ylabel('Power (dB)'); 
        grid on;
        title([plot_title ' | PSD (0–50 Hz)']);
    end
end
