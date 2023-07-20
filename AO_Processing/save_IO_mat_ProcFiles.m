function mat_ProcFiles = save_IO_mat_ProcFiles(mat_filelist, RawDataDir, ProcDataDir)

   


for i = 1:height(mat_files)

     cd(RawDataDir)

    tmpFilename = mat_filelist.FullFile{mi};
    matFileInfo = matfile(tmpFilename);
    matFileVars1 = whos(matFileInfo);
    matFileVars2 = {matFileVars1.name};
    
    % use fuction(s) to isolate fields of interest 
        % look at code in GitHub repo: save_DLCprocFiles_er
        
        % 2. Find all LFP 'CLFP'
        % 3. Find all mLFP 'CMacro_LFP'
        % 4. Find all TTL 'CDIG'
        % Optional Find EMG when used
    % 1. Find all spike files 'CSPK'
    [outStructSPK] = getFILEinfo('CSPK',matFileVars2,tmpFilename);
    % 2. Find all LFP
    [outStructLFP] = getFILEinfo('CLFP',matFileVars2,tmpFilename);
    % 3. Find all mLFP
    [outStructMLFP] = getFILEinfo('CMacro_LFP',matFileVars2,tmpFilename);
    % 4. Find all TTL
    [outStructTTL] = getFILEinfo('CDIG',matFileVars2,tmpFilename);
    % Optional Find EMG when used
    
    % create new struct containing fields of interest

    % save into one struct

    % save into new directory with new name

    % save(filename, new struct var) % filename ~ ProcDataDir combined with mat fileame

    % 1. Find all spike files
    % Save into one Struct
    dlcDepths.Spike = outStructSPK;
    dlcDepths.LFP = outStructLFP;
    dlcDepths.MLFP = outStructMLFP;
    dlcDepths.TTL = outStructTTL;

    % Save into new directory with new name
    sAVEname = ['DLCao_',tmpFilename];
    cd(dlcDATAdir)
    save(sAVEname,'dlcDepths');
end