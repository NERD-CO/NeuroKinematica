function T = build_ZETA_LMM_Table_MUA_NoREST(MasterZETA_MUA)

% build_ZETA_LMM_Table_MUA_NoREST
% Build LMM-ready table for MUA ZETA outcomes, excluding REST.
%
% Reference levels:
%   MoveType ref = HAND OC
%   Depth ref    = b

T = build_ZETA_LMM_Table_MUA(MasterZETA_MUA);

% Drop REST and re-level MoveType
T = T(T.MoveType ~= categorical("REST"), :);
T.MoveType = removecats(T.MoveType);
T.Depth    = removecats(T.Depth);
T.Subject  = categorical(T.Subject);

T.MoveType = categorical(string(T.MoveType), {'HAND OC','HAND PS','ARM EF'});  % HAND OC: ref
T.Depth    = categorical(string(T.Depth),    {'b','c','t'});                   % b (ventral): ref

end
