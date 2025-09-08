function LFP_PostSpecInterp_Tbl = apply_LFPfilters_HP_LP_resample(LFP_PostSpecInterp_Tbl, lfpCols, Fs_in)
    if ischar(lfpCols), lfpCols = string(lfpCols); end
    if iscell(lfpCols), lfpCols = string(lfpCols); end
    have = string(LFP_PostSpecInterp_Tbl.Properties.VariableNames);
    lfpCols = lfpCols(ismember(lfpCols, have));
    if isempty(lfpCols), warning('No LFP cols found.'); return; end

    for c = 1:numel(lfpCols)
        col = lfpCols(c);
        out = col + "_preproc500";
        LFP_PostSpecInterp_Tbl.(out) = cell(height(LFP_PostSpecInterp_Tbl),1);
        for r = 1:height(LFP_PostSpecInterp_Tbl)
            x = LFP_PostSpecInterp_Tbl.(col){r};
            if isempty(x) || ~isnumeric(x) || ~isvector(x)
                LFP_PostSpecInterp_Tbl.(out){r} = x; continue
            end
            try
                LFP_PostSpecInterp_Tbl.(out){r} = preprocess_lfp_vector(x, Fs_in);
            catch ME
                warning('Row %d, %s: %s', r, col, ME.message);
                LFP_PostSpecInterp_Tbl.(out){r} = x(:);
            end
        end
    end
end
