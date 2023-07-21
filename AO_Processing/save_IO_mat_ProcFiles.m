function mat_ProcFiles = save_IO_mat_ProcFiles(mat_filelist, Case_DataDir)

% hardcode directories
% IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';  % directory where all IO data is located
RawDataDir = [Case_DataDir, filesep, 'Raw Electrophysiology MATLAB'];                    % directory where raw MATLAB data files are located (case-specific)
ProcDataDir = [Case_DataDir, filesep, 'Processed Electrophysiology'];                    % directory where processed MATLAB data should be saved (case-specific)

% navigate to directory where raw MATLAB data files are located
cd(RawDataDir)

% extract relevant info from relevant .mat files in mat_filelist
for i = 1:height(mat_filelist)

    tmpFilename = mat_filelist{i};          % retrieve i-th .mat filename from mat_filelist
    matFileInfo = matfile(tmpFilename);     % create matfile object representing the .mat file specified by tmpFilename
    matFileVars1 = whos(matFileInfo);       % use 'whos' function to get info about variables stored in .mat file represented by matFileInfo
    matFileVars2 = {matFileVars1.name};     % extract names of all variables in the .mat file and store them in cell array matFileVars2.

    % isolate fields of interest - look at code in GitHub repo: save_DLCprocFiles_er

    % list conditions
    ftypes = {'CSPK', 'CLFP', 'CMacro_LFP', 'CDIG'};

    for f = 1:4
        switch ftypes{f}
            case 'CSPK'
                % 1. Find all spike files
                % find relevant fields
                % save to outStruct

            case 'CLFP'
                % 2. Find all LFP
                % find relevant fields
                % save to outStruct

            case 'CMacro_LFP'
                % 3. Find all mLFP
                % find relevant fields
                % save to outStruct

            case 'CDIG'
                % 4. Find all TTL
                % find relevant fields
                % save to outStruct
        end

    end

    % 5. Find EMG when used
    % 6. Find ACC when used

    % create new struct containing fields of interest

    % save into one struct
    ProcEphys.Spike = outStructSPK;
    ProcEphys.LFP = outStructLFP;
    ProcEphys.MLFP = outStructMLFP;
    ProcEphys.TTL = outStructTTL;

    % save into new directory with new name
    saveName = ['EphysIO_',tmpFilename];
    cd(ProcDataDir)
    save(saveName,'ProcEphys');

end

