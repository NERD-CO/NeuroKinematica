function [offset_spike_samples, meta_Offset_spk] = useOffset_spikes(TTL_fs, AO_spike_fs, pre_offset_ms, useOffset)

% useOffset_LFP
% Returns # of samples to pre-pad in spike sample domain based on pre-offset time. 
% If useOffset == false or pre_offset_ms <= 0, Returns 0.

% OUT:
%   offset_spike_samp : integer samples in spike domain to subtract from start
%   meta_Offset     : struct with .ms, .seconds, .ttl_samples, .spike_samples

    if nargin < 4 || isempty(useOffset), useOffset = true; end

    meta_Offset_spk.ms = pre_offset_ms;
    meta_Offset_spk.seconds = pre_offset_ms/1000;
    if ~useOffset || pre_offset_ms <= 0
        meta_Offset_spk.ttl_samples = 0;
        meta_Offset_spk.spk_samples = 0;
        offset_spike_samples  = 0;
        return
    end

    % For completeness, keep both TTL and LFP counts available in meta
    meta_Offset_spk.ttl_samples = round(TTL_fs * meta_Offset_spk.seconds);
    meta_Offset_spk.spk_samples = round(AO_spike_fs * meta_Offset_spk.seconds);

    % What we actually need when indexing LFP data:
    offset_spike_samples = meta_Offset_spk.spk_samples;
end