function [offset_LFP_samples, meta_Offset] = useOffset_LFP(TTL_fs, AO_LFP_fs, offset_ms, useOffset)

% useOffset_LFP
% Returns # of samples to pre-pad in LFP domain based on pre-offset time. 
% If useOffset == false or offset_ms <= 0, Returns 0.

% OUT:
%   offset_LFP_samp : integer samples in LFP domain to subtract from start
%   meta_Offset     : struct with .ms, .seconds, .ttl_samples, .lfp_samples

    if nargin < 4 || isempty(useOffset), useOffset = true; end

    meta_Offset.ms = offset_ms;
    meta_Offset.seconds = offset_ms/1000;
    if ~useOffset || offset_ms <= 0
        meta_Offset.ttl_samples = 0;
        meta_Offset.lfp_samples = 0;
        offset_LFP_samples  = 0;
        return
    end

    % For completeness, keep both TTL and LFP counts available in meta
    meta_Offset.ttl_samples = round(TTL_fs * meta_Offset.seconds);
    meta_Offset.lfp_samples = round(AO_LFP_fs * meta_Offset.seconds);

    % What we actually need when indexing LFP data:
    offset_LFP_samples = meta_Offset.lfp_samples;
end
