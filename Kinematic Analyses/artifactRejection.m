%% functions

% artifact rejection function
% input table:
%   1. o.g. data: outDATA
% plot
% method
% loop through each marker
% per dlc label, plot raw data
% output fraction of frames that fall below CT
% confidence threshhold ~0.75
% if frame loss% < 5-10%,
% elseif ~~~~, increase/decrease CT
% replace values in rejected frames with interp value
% finder border accepted frame valuse --> insert interp value b/t
% output tables
% 2. decision per frame - binary status (accept/reject): outData Index
% 3. interpolated / cleaned data: outData interp
% plot
% save(current file, index, interp, '-append')
