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
