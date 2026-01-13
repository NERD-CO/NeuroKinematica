function OUT = run_IFR_LMM_Pipeline(MasterTbl, varargin)
% run_IFR_LMM_Pipeline
% Unified IFR-LMM pipeline for SU and MUA tables (with/without REST).
%
% Inputs:
%   MasterTbl : MasterZETA (SU) or MasterZETA_MUA (MUA)
%
% Key options:
%   'DataType'     : "SU" (default) or "MUA"
%   'IncludeREST'  : true (default) or false (NoREST pipeline)
%   'DVList'       : string array of DVs to run
%   'SavePath'     : folder to write CSV outputs ("" = no saving)
%   'PAdjust'      : "bonferroni"/"holm"/"fdr"/"none" (passed to planned contrasts)
%   'AdjustScope'  : "all" or "family"
%   'DropReverse'  : true/false for makeLMM_ContrastsTable
%
% Returns:
%   OUT struct with:
%     .T           analysis table
%     .LME         output from runIFR_LMMs
%     .C           struct of planned-contrast tables by DV
%     .JNE         struct of manuscript tables by DV

p = inputParser;
p.addParameter('DataType', "SU", @(x) isstring(x) || ischar(x));
p.addParameter('IncludeREST', true, @islogical);

p.addParameter('DVList', string.empty, @(x) isstring(x) || ischar(x));
p.addParameter('SavePath', "", @(x) isstring(x) || ischar(x));

p.addParameter('PAdjust', "holm", @(x) isstring(x) || ischar(x));
p.addParameter('AdjustScope', "family", @(x) isstring(x) || ischar(x));
p.addParameter('DropReverse', true, @islogical);

p.parse(varargin{:});
U = p.Results;

dataType = upper(string(U.DataType));
includeREST = U.IncludeREST;

% -----------------------------
% Default DV list
% -----------------------------
if isempty(U.DVList)
    if dataType == "SU"
        dvList = ["IFR_mean_baselineNorm", "IFR_max_baselineNorm" , ...
            "IFR_peakLatency","IFR_peakOnset_Latency"];
        % exclude: "IFR_mean_Znorm", "IFR_mean_Hz", "IFR_max_Hz", ...
    else
        % MUA default: 
        dvList = ["IFR_mean_baselineNorm", "IFR_max_baselineNorm"];
        % exclude: "IFR_mean_Znorm", "IFR_mean_Hz", "IFR_max_Hz"
    end
else
    dvList = string(U.DVList);
end

% -----------------------------
% Build the correct analysis table
% -----------------------------
switch dataType
    case "SU"
        if includeREST
            T = build_IFR_LMM_Table(MasterTbl); % MasterZETA
            % Ensure explicit category order (REST ref; Depth b ref)
            T.MoveType = categorical(string(T.MoveType), {'REST','HAND OC','HAND PS','ARM EF'}); % REST ref
        else % No REST
            T = build_IFR_LMM_Table_NoREST(MasterTbl);
            % Ensure explicit category order (HAND OC ref; Depth b ref)
            T.MoveType = categorical(string(T.MoveType), {'HAND OC','HAND PS','ARM EF'}); % HAND OC ref
        end

    case "MUA"
        if includeREST
            T = build_IFR_LMM_Table_MUA(MasterTbl); % MasterZETA_MUA
            T.MoveType = categorical(string(T.MoveType), {'REST','HAND OC','HAND PS','ARM EF'}); % REST ref
        else % No REST
            T = build_IFR_LMM_Table_MUA_NoREST(MasterTbl);
            T.MoveType = categorical(string(T.MoveType), {'HAND OC','HAND PS','ARM EF'}); % HAND OC ref
        end

    otherwise
        error('Unknown DataType: %s (use "SU" or "MUA")', dataType);
end

% Shared coding
T.Depth   = categorical(string(T.Depth), {'b','t','c'}); % b ref
T.Subject = categorical(string(T.Subject));

% Safety cleanup (prevents ghost categories)
if ~includeREST
    if any(string(categories(T.MoveType))=="REST")
        T = T(T.MoveType ~= categorical("REST"), :);
    end
    T.MoveType = removecats(T.MoveType);
    T.Depth    = removecats(T.Depth);
end


% -----------------------------
% Fit LMEs (DV-agnostic)
% -----------------------------
LME = run_LMMs(T, dvList);

% -----------------------------
% Planned contrasts per DV
% -----------------------------
C = struct();
JNE = struct();
CoeffNames = struct();

for i = 1:numel(dvList)
    dv = dvList(i);

    if ~isfield(LME, dv) || ~isfield(LME.(dv), 'lme') || isempty(LME.(dv).lme)
        fprintf('[SKIP contrasts] DV=%s (no model fit)\n', dv);
        continue;
    end

    lme = LME.(dv).lme;
    CoeffNames.(dv) = string(lme.CoefficientNames(:));   

    if includeREST
        % REST-inclusive planned contrasts
        Cdv = runLMM_PlannedContrasts(lme, ...
            'DVLabel', dv, ...
            'PAdjust', string(U.PAdjust), ...
            'AdjustScope', string(U.AdjustScope));
    else
        % NoREST planned contrasts
        Cdv = runLMM_PlannedContrasts_NoREST(lme, ...
            'DVLabel', dv, ...
            'PAdjust', string(U.PAdjust), ...
            'AdjustScope', string(U.AdjustScope));
    end

    C.(dv) = Cdv;
    JNE.(dv) = makeLMM_ContrastsTable(Cdv, 'DropReverse', U.DropReverse);
    

    % -----------------------------
    % Optional saving
    % -----------------------------
    if strlength(string(U.SavePath)) > 0
        SavePath = char(U.SavePath);
        if ~exist(SavePath,'dir'), mkdir(SavePath); end

        tagData = char(dataType);
        tagRest = "WithREST";
        if ~includeREST, tagRest = "NoREST"; end
        % new
        dvTag   = regexprep(char(dv), '[^A-Za-z0-9_]+', '_');

        fnameFull = sprintf('%s_%s_%s_PlannedContrasts_FULL.csv', tagData, tagRest, dvTag);
        fnameJNE  = sprintf('%s_%s_%s_PlannedContrasts_JNE.csv',  tagData, tagRest, dvTag);

        writetable(Cdv, fullfile(SavePath, fnameFull));
        writetable(JNE.(dv), fullfile(SavePath, fnameJNE));
    end
end

% -----------------------------
% Combine ALL DVs tables for easy review + optional saving
% -----------------------------
All_FULL = table();
All_JNE  = table();

dvNames = string(fieldnames(C));
for i = 1:numel(dvNames)
    dv = dvNames(i);
    if ~isempty(C.(dv)),   All_FULL = [All_FULL; C.(dv)]; end 
    if ~isempty(JNE.(dv)), All_JNE  = [All_JNE;  JNE.(dv)]; end 
end

if strlength(string(U.SavePath)) > 0 && height(All_FULL) > 0
    tagData = char(dataType);
    tagRest = "WithREST"; if ~includeREST, tagRest = "NoREST"; end
    writetable(All_FULL, fullfile(char(U.SavePath), sprintf('%s_%s_ALLDVs_PlannedContrasts_FULL.csv', tagData, tagRest)));
    writetable(All_JNE,  fullfile(char(U.SavePath), sprintf('%s_%s_ALLDVs_PlannedContrasts_JNE.csv',  tagData, tagRest)));
end

% -----------------------------
% Assemble outputs
% -----------------------------
OUT = struct();
OUT.DataType     = dataType;
OUT.IncludeREST  = includeREST;
OUT.DVList       = dvList;
OUT.T            = T;
OUT.LME          = LME;
OUT.C            = C;
OUT.JNE          = JNE;
OUT.CoeffNames   = CoeffNames;
OUT.All_FULL     = All_FULL;
OUT.All_JNE      = All_JNE;

fprintf('IFR LMM pipeline finished | %s | IncludeREST=%d | DVs=%d | N=%d\n', ...
    dataType, includeREST, numel(dvList), height(T));

end
