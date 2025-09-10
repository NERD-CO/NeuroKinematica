function [epochDur_LFP_samples, meta_epochDur] = UniformEpochs_LFP(TTL_fs, AO_LFP_fs, epochDur_ms, UniformEpochs)

% UniformEpochs_LFP
% Returns # of samples to pad in LFP domain (after trial start)
% based on uniform epoch duration time.
% If UniformEpochs_LFP == false or epochDur_ms <= 0, Returns 0.

% OUT:
%   epochDur_LFP_samples : integer samples in LFP domain to add from start
%   meta_epochDur        : struct with .ms, .seconds, .ttl_samples, .lfp_samples

if nargin < 4 || isempty(UniformEpochs), UniformEpochs = true; end

meta_epochDur.ms = epochDur_ms;
meta_epochDur.seconds = epochDur_ms/1000;

if ~UniformEpochs || epochDur_ms <= 0
    meta_epochDur.ttl_samples = 0;
    meta_epochDur.lfp_samples = 0;
    epochDur_LFP_samples  = 0;
    return
end

% For completeness, keep both TTL and LFP counts available in meta
meta_epochDur.ttl_samples = round(TTL_fs * meta_epochDur.seconds);
meta_epochDur.lfp_samples = round(AO_LFP_fs * meta_epochDur.seconds);

% What we actually need when indexing LFP data:
epochDur_LFP_samples = meta_epochDur.lfp_samples;

end
