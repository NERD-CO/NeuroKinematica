fs = 250;

[pxxMM,fMM] = pwelch(pointX2sm,250,125,250,250);
p2bmm = pow2db(pxxMM);
p2bMMs = smoothdata(p2bmm,'gaussian',2);
fMMi = fMM > 1 & fMM < 10;
fMMf = fMM(fMMi);
pxxMf = p2bMMs(fMMi);

% figure;
% plot(fMMf,pxxMf)




% plot(f,10*log10(pxx))
% 
% xlabel('Frequency (Hz)')
% ylabel('PSD (dB/Hz)')
% 
% xlim([1 10])


%%

[Pupvxx,upVfxx] = pwelch(transpose(movLFP), hanning(250), 125, 256, 250, 'onesided');
uVp_t = sqrt(Pupvxx).*rms(hanning(250)).*sqrt(2).*2.*250/256;
% PxxP = 10*log10(Pxx);
% uPv_A = uVp_t;

uPv_As = smoothdata(uVp_t,'gaussian',5);
% pxxps = smoothdata(PxxP,'gaussian',5);

upv50 = upVfxx < 50;
upVfxxU = upVfxx(upv50);
uVp_tU = uPv_As(upv50);


figure;
plot(upVfxxU,uVp_tU)
% xlim([0 50])

% figure;
% plot(pxxps)
% xlim([0 50])
hold on
plot(uPvRESTfreq , uPvRESTpow)
legend('epoch','rest')


%%
% PCA plot of movement kin

% 3D plot of movement kin

% cool animation plot

% wavelet cluster analysis ---- see pubicaiton 
