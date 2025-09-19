function [varExtract] = extractVarName(inMfname,inMvar)

tmpLoadF = load(inMfname,inMvar);
tmpLoadFns = fieldnames(tmpLoadF);
varExtract = tmpLoadF.(tmpLoadFns{1});

end