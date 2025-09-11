function [start_indices, end_indices] = findConsecutiveIndices_fun(logical_array)
% Find start and end indices of consecutive true regions in a logical array
d = diff([0; logical_array(:); 0]);
start_indices = find(d > 0);
end_indices = find(d < 0) - 1;
end