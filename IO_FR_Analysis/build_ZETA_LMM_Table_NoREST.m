function T = build_ZETA_LMM_Table_NoREST(MasterZETA)

% Build LMM-ready table for ZETA outcomes, excluding REST.
%
% Reference levels:
%   MoveType ref = HAND OC
%   Depth ref    = b

T = build_ZETA_LMM_Table(MasterZETA);

% Drop REST and re-level MoveType
T = T(T.MoveType ~= categorical("REST"), :);
T.MoveType = removecats(T.MoveType);
T.Depth    = removecats(T.Depth);
T.Subject  = categorical(T.Subject);

T.MoveType = categorical(string(T.MoveType), {'HAND OC','HAND PS','ARM EF'}); % ref: HAND OC
T.Depth    = categorical(string(T.Depth),    {'b','t','c'});                  % ref: b (ventral)

end
