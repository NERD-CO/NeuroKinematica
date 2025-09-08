function LFPsPerMove_T = wrapper_applySpectrumInterpolation_LFPs(LFPsPerMove_T, lfpCols, Fs, Fl, neighborsToSample, neighborsToReplace, makeNewCols) 
% Runs spectrumInterpolation() on each 1D vector in given LFP channel columns. 

% Normalize lfpCols to MATLAB string type, whether you pass a char scalar 
% or a cell array of char vectors. 
if ischar(lfpCols), lfpCols = string(lfpCols); end 
if iscell(lfpCols), lfpCols = string(lfpCols); end 
% Default: if you didn't pass makeNewCols, create new columns (keeps raw data intact) 
if nargin < 7 || isempty(makeNewCols) 
    makeNewCols = true; 
end 

have = string(LFPsPerMove_T.Properties.VariableNames); 
lfpCols = lfpCols(ismember(lfpCols, have)); 
if isempty(lfpCols) 
    warning('No matching LFP columns found to filter. Skipping.'); 
    return 
end 

% Loop over each LFP_E column 
for col_i = 1:numel(lfpCols) 
    colName = lfpCols(col_i);  % current LFP column name (e.g., "LFP_E1"). 
    if makeNewCols 
        outName = colName + "_filt"; % e.g., "LFP_E1_filt" 
    else 
        outName = colName; % overwrite original 
    end 

    % Initialize the output variable for filtered LFPs 
    LFPsPerMove_T.(outName) = cell(height(LFPsPerMove_T),1); % cell column with the same height as T. 
    
    % Iterate over each row (movement trial seg) fill each row with a filtered vector 
    for row_i = 1:height(LFPsPerMove_T) 
        raw_vec = LFPsPerMove_T.(colName){row_i}; % Extract raw vector from the cell 
        
        % Only process numeric vectors; otherwise pass through silently 
        if isempty(raw_vec) || ~(isnumeric(raw_vec) && isvector(raw_vec)) 
            LFPsPerMove_T.(outName){row_i} = raw_vec; 
            continue 
        end 

        raw_vec = double(raw_vec(:)); % ensures a column vector as a double 
        
        % Apply spectrum interpolation to the raw vector (once per LFP col vec) 
        try 
            vec_filt = spectrumInterpolation(raw_vec, Fs, Fl, neighborsToSample, neighborsToReplace); 
        catch ME 
            warning('Row %d, %s: interpolation failed (%s). Keeping original.', row_i, colName, ME.message); 
            vec_filt = raw_vec; 
        end 
        LFPsPerMove_T.(outName){row_i} = vec_filt; 
    end 
end 
end 
