figure;

subplot(2,1,1)
plot(tmpXwave,tmpYwaveNan,...
    'LineWidth',0.5,"Color",'k')

spkDat2plot = app.spikeDataRaw(app.EpochInds(app.CurEpoch,1):...
    app.EpochInds(app.CurEpoch,2));
colORvolt = [0.9412 0.9412 0.9412 0.8];

testPLOTmain.spkDataall = spkDat2plot;
testPLOTmain.color = colORvolt;

save('testFIG.mat','testPLOTmain')

testPLOT.tmpXwave = tmpXwave;
testPLOT.tmpYwave = tmpYwaveNan;

save('testFIG.mat','testPLOT','-append')


spkDat2plot = app.spikeDataRaw(app.EpochInds(1,1):...
    app.EpochInds(1,2));
colORvolt = [0.9412 0.9412 0.9412 0.8];

subplot(2,1,2)
plot(spkDat2plot,'Color',colORvolt);

%%

figure;


testALL = lowpass(double(testPLOTmain.spkDataall),5000,44000);
smALL = smoothdata(testALL,'gaussian',22);
subplot(2,1,1)
plot(testPLOT.tmpXwave,testPLOT.tmpYwave,...
    'LineWidth',0.5,"Color",'k')
subplot(2,1,2)
plot(smALL,'Color','k');
hold on
plot(testPLOT.tmpXwave,testPLOT.tmpYwave,...
    'LineWidth',0.5,"Color",'r')

%% 

tmpRAW = ProcEphys.Spike.E1.rawData;

lpRAW = lowpass(double(tmpRAW),5000,44000);
smRAW = smoothdata(lpRAW,'gaussian',22);

ProcEphys.Spike.E1.rawData = int16(smRAW);


%%

matfiles1 = dir('*.mat');
matfiles2 = {matfiles1.name};

for mi = 1:length(matfiles2)

    tmpMfile = matfiles2{mi};
    load(tmpMfile,'ProcEphys')

    tmpEls = fieldnames(ProcEphys.Spike);

    for eti = 1:length(tmpEls)

        tmpDATA = ProcEphys.Spike.(tmpEls{eti}).rawData;

        lpRAW = lowpass(double(tmpDATA),5000,44000);
        smRAW = smoothdata(lpRAW,'gaussian',22);

        ProcEphys.Spike.(tmpEls{eti}).rawData = int16(smRAW);
    end
    save(tmpMfile,"ProcEphys")
end


%% LFP

cd('E:\Dropbox\MUD_STUDY\IO_MUD1\MUD_1\KC_NAc_Right\Processed Electrophysiology')

tmpFile = 'Processed_RT1D2.572F0001.mat';

load(tmpFile)

testLFP = double(ProcEphys.LFP.E1.rawData);
testFS = round(ProcEphys.LFP.E1.Hz*1000);

%%
% [pxx , fxx] = pspectrum(testLFP,testFS,"power","FrequencyLimits",[1 100],...
    % 'FrequencyResolution',1);

[pxx , fxx] = pspectrum(testLFP,testFS,"power","FrequencyLimits",[1 55],...
    'FrequencyResolution',1.25);

save('testFreq.mat','pxx','fxx')

pxx2 = pow2db(pxx);

%%

matLIST = dir('*.mat');
matLIST2 = {matLIST.name};

allE1 = zeros(4401,length(matLIST2));
allE2 = zeros(4401,length(matLIST2));

for mi = 1:length(matLIST2)

    load(matLIST2{mi},'ProcEphys');

    eFields = fieldnames(ProcEphys.LFP);

    for eii = 1:length(eFields)

        tMPrawE1 = double(ProcEphys.LFP.E1.rawData);
        tMPrawE2 = double(ProcEphys.LFP.E2.rawData);

        [pxxE1 , ~] = pspectrum(tMPrawE1,testFS,"power","FrequencyLimits",[1 55],...
            'FrequencyResolution',1.25);

        [pxxE2 , ~] = pspectrum(tMPrawE2,testFS,"power","FrequencyLimits",[1 55],...
            'FrequencyResolution',1.25);

        allE1(:,mi) = pow2db(pxxE1);

        allE2(:,mi) = pow2db(pxxE2);

    end
end

%%

% find outliers

e1mean = mean(allE1);
e1std = std(allE1);

kurtVal = zeros(width(allE1),1);
for ki = 1:width(allE1)
    kurtVal(ki) = kurtosis(allE1(:,ki));
end

% kurtVal(channel,1) = kurtosis(data,:)
kurtMAD = mad(kurtVal,1);
kurtMED = median(kurtVal);
kurtZ = ((kurtVal - kurtMED)./kurtMAD); % Z-scored kurtosis
I = kurtZ > 3 | kurtZ < -3;
plot(allE1(:,~I),'k')
hold on
plot(allE1(:,I),'r')
% xlim([1 55])

%%

e1mean = mean(allE1,2);
e1std = std(allE1,[],2);
thresh = e1mean + (e1std*1.25);

fracAbove = zeros(width(allE1),1);
for ki = 1:width(allE1)
    fracAbove(ki) = sum(allE1(:,ki) > thresh) / height(allE1);
end
close all
I = fracAbove > 0.25;

plot(allE1(:,~I),'k')
hold on
plot(allE1(:,I),'r')

%%

keepE1 = allE1(:,~I);

%%
e2mean = mean(allE2,2);
e2std = std(allE2,[],2);
thresh = e1mean + (e2std*1.25);

fracAbove = zeros(width(allE2),1);
for ki = 1:width(allE2)
    fracAbove(ki) = sum(allE2(:,ki) > thresh) / height(allE2);
end
close all
I = fracAbove > 0.25;

plot(allE2(:,~I),'k')
hold on
plot(allE2(:,I),'r')

%%

keepE2 = allE2(:,~I);
%%
plot(keepE2)

