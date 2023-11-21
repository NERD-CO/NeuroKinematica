% Let's do Gaussian smoothing with arbitrary underlying sample rate.

% xin = K.Time;
% yin = fillpos;
% xout = K.Time; to just do smoothing at original timepoints, without resampling or interpolation, set xout the same as xin
% gkern = 0.05; standard deviation of gaussian (aka width gaussian). standard value to use is 50 ms (0.05s). 
% sigmax = 4; sigmax is how far to extend tails of the gaussian (essentially how many samples away from the current sample to consider)

function yout = GaussSmooth_Arbitrary(xin, yin, xout, gkern, sigmax)

% xin is the timestamps of the original samples (in seconds)
% yin is the unfiltered raw signal sampled at xin
% xout is the timestamps of the desired output (in seconds, same clock as xin)
% gkern is the size of the Gaussian kernel (in seconds)
% sigmax is the number of sigmas to look forward and backwards to

% Right now if there are nans, the function tries to extend the good signal
% as far as possible into the nan regions, essentially weighting any
% non-nan values more strongly to place a smooth value at times when the
% original signal may have been nan.

nout = length(xout);
ny = size(yin,2);
lookout = gkern*sigmax;

yout = nan(nout,ny);

for i = 1:nout
    xrel = xin-xout(i);
    kmask = abs(xrel) < lookout;
    kpoints = normpdf(xrel(kmask), 0, gkern);
    yinslice = yin(kmask,:);
    notnan = ~isnan(yinslice(:,1));
    kpoints = kpoints(notnan);
    kpoints = kpoints/sum(kpoints);
    if ~isempty(kpoints)
        yout(i,:) = yinslice(notnan,:)'*kpoints;
    end
end