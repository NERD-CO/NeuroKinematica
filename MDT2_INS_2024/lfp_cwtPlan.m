[cfs,frq] = cwt(cleanLFP,250,'FrequencyLimits',[2 100],'VoicesPerOctave',48);
t = 0:1/250:(1796-1)/250;

absCFS = abs(cfs);

imagesc(t,frq,abs(cfs));
axis xy
% set(gca, 'YScale', 'log');
% imagesc(cfs)

betaIND = frq > 13 & frq < 30;
betaBAND = absCFS(betaIND,:)

betaBmean = mean(betaBAND,1)
betaBsm = smoothdata(betaBmean,'gaussian',30)
plot(betaBmean)
hold on
plot(betaBsm)