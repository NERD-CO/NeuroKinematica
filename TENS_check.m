% TENS check

% CD to location of JSON file
dir = '';
cd(dir);

jsonSelect = uigetfile('*.json');
jsrt = jsondecode(fileread(jsonSelect));
bstd = jsrt.BrainSenseTimeDomain;
bstdTab = struct2table(bstd);
tmpLFP = bstdTab.TimeDomainData{1} % left
plot(tmpLFP)
hold on
tmpLFP = bstdTab.TimeDomainData{2} % right
plot(tmpLFP)