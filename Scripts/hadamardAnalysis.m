close all;
%% Load and Present Hadamard Results
load('D:\Results\Results-14-Jul-2019 16-35-15-Consistency-wHad-2msQuant-100mW-2sPiez\Results.mat')
load('D:\Results\Results-14-Jul-2019 16-35-15-Consistency-wHad-2msQuant-100mW-2sPiez\Vars.mat')

hadRes = res;
hadVars = curVars;

figure();
errorbar(curVars.scan.timeFrames, res.phiFrame(7,:), res.phiFrameStd(7,:))

figure();
errorbar(curVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, res.phiFrame(:,1), res.phiFrameStd(:,1))

%% Load and Present Naive Results
close all;
load('D:\Results\Results-14-Jul-2019 11-36-57-Consistency-noHad-2msQuant-100mW-2sPiez\Results.mat')
load('D:\Results\Results-14-Jul-2019 11-36-57-Consistency-noHad-2msQuant-100mW-2sPiez\Vars.mat')

naiveRes = res;
naiveVars = curVars;

figure();
errorbar(curVars.scan.timeFrames, res.phiFrame(7,:), res.phiFrameStd(7,:))

figure();
errorbar(curVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, res.phiFrame(:,1), res.phiFrameStd(:,1))

%% Compare
% 
% figure();
% subplot(1,2,1)
% errorbar(hadVars.scan.timeFrames, hadRes.phiFrame(1,:), hadRes.phiFrameStd(1,:)); hold on
% errorbar(naiveVars.scan.timeFrames, naiveRes.phiFrame(1,:), naiveRes.phiFrameStd(1,:));
% subplot(1,2,2)
% errorbar(hadVars.scan.timeFrames, hadRes.phiFrame(1,:), hadRes.phiFrameStd(1,:)); hold on
% errorbar(naiveVars.scan.timeFrames, naiveRes.phiFrame(1,:), naiveRes.phiFrameStd(1,:));
% legend('with Hadamard', 'w/o Hadamard');

figure();
errorbar(hadVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, hadRes.phiFrame(:,1), hadRes.phiFrameStd(:,1)); hold on
errorbar(naiveVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, naiveRes.phiFrame(:,1), naiveRes.phiFrameStd(:,1)); hold off
legend('with Hadamard', 'w/o Hadamard');
title('Frame Duration: 2[s]');


figure();
subplot(2,2,1)
errorbar(hadVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, hadRes.phiFrame(:,1), hadRes.phiFrameStd(:,1)); hold on
errorbar(naiveVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, naiveRes.phiFrame(:,1), naiveRes.phiFrameStd(:,1)); hold off
legend('with Hadamard', 'w/o Hadamard');
title('Frame Duration: 2[s]');

subplot(2,2,2)
errorbar(hadVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, hadRes.phiFrame(:,2), hadRes.phiFrameStd(:,2)); hold on
errorbar(naiveVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, naiveRes.phiFrame(:,2), naiveRes.phiFrameStd(:,2)); hold off
legend('with Hadamard', 'w/o Hadamard');
title('Frame Duration: 4[s]');

subplot(2,2,3)
errorbar(hadVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, hadRes.phiFrame(:,3), hadRes.phiFrameStd(:,3)); hold on
errorbar(naiveVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, naiveRes.phiFrame(:,3), naiveRes.phiFrameStd(:,3)); hold off
legend('with Hadamard', 'w/o Hadamard');
title('Frame Duration: 6[s]');

subplot(2,2,4)
errorbar(hadVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, hadRes.phiFrame(:,4), hadRes.phiFrameStd(:,4)); hold on
errorbar(naiveVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, naiveRes.phiFrame(:,4), naiveRes.phiFrameStd(:,4)); hold off
legend('with Hadamard', 'w/o Hadamard');
title('Frame Duration: 8[s]');

SNRHad   = hadRes.phiFrame.^2 ./ hadRes.phiFrameStd.^2;
SNRNaive = naiveRes.phiFrame.^2 ./ naiveRes.phiFrameStd.^2;

figure();
subplot(2,2,1)
plot(hadVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, SNRHad(:,1)); hold on;
plot(naiveVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, SNRNaive(:,1));hold off
legend('with Hadamard', 'w/o Hadamard');
title('Frame Duration: 2[s]');

subplot(2,2,2)
plot(hadVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, SNRHad(:,2)); hold on;
plot(naiveVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, SNRNaive(:,2));hold off
legend('with Hadamard', 'w/o Hadamard');
title('Frame Duration: 4[s]');

subplot(2,2,3)
plot(hadVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, SNRHad(:,3)); hold on;
plot(naiveVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, SNRNaive(:,3));hold off
legend('with Hadamard', 'w/o Hadamard');
title('Frame Duration: 6[s]');

subplot(2,2,4)
plot(hadVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, SNRHad(:,4)); hold on;
plot(naiveVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, SNRNaive(:,4));hold off
legend('with Hadamard', 'w/o Hadamard');
title('Frame Duration: 8[s]');

%% Search for best Delay
hadData = load('D:\Results\Results-14-Jul-2019 16-35-15-Consistency-wHad-2msQuant-100mW-2sPiez\rawData\F2.00S1-rawData.mat');
hadVars = load('D:\Results\Results-14-Jul-2019 16-35-15-Consistency-wHad-2msQuant-100mW-2sPiez\rawData\F2.00S1-vars.mat');

hadData = hadData.res;
hadVars = hadVars.curVars;

ao = acoustoOptics();
uVars = hadVars.acoustoOptics.uVars;
uVars.exportRawData = true;

gReq = algoGraphics.createGraphicsRunVars();
gReq.ch    = 1;
gReq.zIdx  = 8;
gReq.quant = 1;
gReq.intExt = 'int';

names = fieldnames(gReq.validStruct);
for i=1:length(names)
    gReq.validStruct.(names{i}) = false;
end

uVars.gReq = gReq;

figure()
ax = axes();

inPhantSamp = hadVars.acoustoOptics.algoVars.samples.inPhantomPropSamples;

peakMean = NaN*ones(1,inPhantSamp);
peakStd = NaN*ones(1,inPhantSamp);

for i = 1:hadVars.acoustoOptics.algoVars.samples.inPhantomPropSamples
    ao.setMeasVars(uVars);
    algoVars = ao.getAlgoVars;
    ao.algo.samples.inPhantomPropSamples = i;
    ao.algo.data = hadData;
    disp(i)
    res = ao.algo.analyse();
    results(i).phi = res.phi;
    results(i).phiStd = res.phiStd;
    
    peakMean(i) = results(i).phi(7);
    peakStd(i)  = results(i).phiStd(7);
    errorbar(ax, 1:i, peakMean(1:i), peakStd(1:i))
    xlim([0,i+1])
    drawnow()
end









