function C = getJNEColorMaps()
% getJNEColorMaps
% Canonical JNE palette + maps for STN depth + MoveType.
%
% Returns struct C with fields:
%   C.purpleShades (3x3)  [t;c;b]
%   C.greenShades  (4x3)  [REST;HAND OC;HAND PS;ARM EF]
%   C.depthColorMap, C.moveColorMap, C.activeMoveColorMap
%   C.depthOrder, C.moveOrder, C.activeMoveOrder
%   C.depthNames, C.fallbackCol

C.purpleShades = ([ ...
    118,42,131; ...
    175,141,195; ...
    231-15, 212-15, 232-15] ./ 255);

C.greenShades = ([ ...
    128,128,128; ...
    217,240,211; ...
    127,191,123; ...
    27,120,55] ./ 255);

C.depthOrder = {'t','c','b'};
C.moveOrder  = {'HAND OC','HAND PS','ARM EF','REST'};
C.activeMoveOrder = {'HAND OC','HAND PS','ARM EF'};

C.depthNames = containers.Map({'t','c','b'}, ...
    {'dorsal STN','central STN','ventral STN'});

C.depthColorMap = containers.Map({'t','c','b'}, ...
    {C.purpleShades(1,:), C.purpleShades(2,:), C.purpleShades(3,:)});

C.moveColorMap = containers.Map({'REST','HAND OC','HAND PS','ARM EF'}, ...
    {C.greenShades(1,:), C.greenShades(2,:), C.greenShades(3,:), C.greenShades(4,:)});

C.activeMoveColorMap = containers.Map({'HAND OC','HAND PS','ARM EF'}, ...
    {C.greenShades(2,:), C.greenShades(3,:), C.greenShades(4,:)});

C.fallbackCol = [0.5 0.5 0.5];
end
