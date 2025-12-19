function LMM_tbl = makeLMM_ContrastsTable(C, varargin)
% Make a publication-ready contrasts table from runIFR_PlannedContrasts_SU output

p = inputParser;
p.addParameter('DropReverse', true, @islogical);
p.addParameter('Alpha', 0.05, @(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:});
U = p.Results;

LMM_tbl = C(:, {'DV','DepthLabel','Contrast','Est_CI','p', 'p_adj'});
LMM_tbl.Properties.VariableNames = {'DV','Depth','Contrast','Estimate_95CI', 'p_raw', 'p_adj'};

% --- Drop reverse-direction duplicates (recommended for manuscript) ---
if U.DropReverse
    % Drop "REST vs ..." lines
    drop1 = startsWith(string(LMM_tbl.Contrast), "REST vs ");

    % Drop explicit reverse active-pair duplicates by keeping only one orientation
    % Rule: keep contrasts where the left side is alphabetically <= right side
    % e.g., keep "ARM EF vs HAND OC" OR keep "HAND OC vs HAND PS" consistently
    isPair = contains(string(LMM_tbl.Contrast), " vs ") & ~contains(string(LMM_tbl.Contrast), "Mean(Active)");
    parts  = split(string(LMM_tbl.Contrast), " vs ");
    left   = parts(:,1); right = parts(:,2);

    % Only apply to active-pair style rows (not REST comparisons)
    drop2 = false(height(LMM_tbl),1);
    for i = 1:height(LMM_tbl)
        if isPair(i) && ~startsWith(left(i),"REST") && ~startsWith(right(i),"REST")
            drop2(i) = left(i) > right(i); % drop reverse ordering
        end
    end

    % Drop the redundant "REST vs Mean(Active)" line
    drop3 = contains(string(LMM_tbl.Contrast), "REST vs Mean(Active)");

    LMM_tbl = LMM_tbl(~(drop1 | drop2 | drop3), :);
end

% --- p formatting + stars ---
p_use = LMM_tbl.p_adj;
p_fmt = strings(height(LMM_tbl),1);
sig   = strings(height(LMM_tbl),1);

for i = 1:height(LMM_tbl)
    if isnan(p_use(i))
        p_fmt(i) = "";
        sig(i)   = "";
    elseif p_use(i) < 0.001
        p_fmt(i) = "<0.001";
        sig(i)   = "***";
    elseif p_use(i) < 0.01
        p_fmt(i) = sprintf('%.3f', p_use(i));
        sig(i)   = "**";
    elseif p_use(i) < U.Alpha % 0.05
        p_fmt(i) = sprintf('%.3f', p_use(i));
        sig(i)   = "*";
    else
        p_fmt(i) = sprintf('%.3f', p_use(i));
        sig(i)   = "";
    end
end

LMM_tbl.p_adj_fmt = p_fmt;
LMM_tbl.Sig = sig;

% Order rows nicely
LMM_tbl = sortrows(LMM_tbl, {'DV','Depth','Contrast'});


end
