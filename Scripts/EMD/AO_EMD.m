close all
clear all
clc

%% Load Data
res(1) = load('D:\WavelengthScan\SNR-Vs.-Quant-PMT-2s-0Hz\19-Jul-2021 17-31-27-PMT-2s-0Hz-q.0.0001.mat');
res(2) = load('D:\WavelengthScan\SNR-Vs.-Quant-PMT-2s-5Hz\19-Jul-2021 16-13-52-PMT-2s-q.0.0020.mat');
res(3) = load('D:\WavelengthScan\SNR-Vs.-Quant-PMT-2s-10Hz\19-Jul-2021 18-56-11-PMT-2s-10Hz-q.0.0001.mat');
res(4) = load('D:\WavelengthScan\SNR-Vs.-Quant-PMT-2s-25Hz\19-Jul-2021 19-24-18-PMT-2s-25Hz-q.0.0001.mat');
res(5) = load('D:\WavelengthScan\SNR-Vs.-Quant-PMT-2s-50Hz-2\19-Jul-2021 20-04-13-PMT-2s-50Hz-q.0.0001.mat');

for i=1:5
    resMeas(i) = res(i).resMeas;
    resVars(i) = res(i).aoVars;
end

% depthAx= resVars(1).measVars.algo.len.depthVecUS*1e3;
% fBar = resVars(1).measVars.algo.freq.fBar*1e-6;

%% Configure Algo Object
close all
% clear resRecalc fittedChAvgFFT fittedChAvgFFTAvg2
legStr = {'0Hz', '5Hz', '10Hz', '25Hz', '50Hz'};
figure()
algo = Algo();
frameTime = [0.5, 0.5, 0.5, 0.5, 0.5];
env = 1000;
posIdx = 102;
reps = 30;
% ax(1) = subplot(1,3,1);
% ax(2) = subplot(1,3,2);
ax(3) = subplot(1,3,3);
ax(3) = subplot(1,1,1);
for i=1:5
    fprintf("i=%d\n", i)
%     for j = 1:reps
%         fprintf("j=%d ", j)
%         algoVars = convertAlgoParams(resVars(i).extVars.algo);
%         algoVars.frameTime = frameTime(i);
%         algo.setVars(algoVars);
%         algoVars = algo.getVars();
%         algo.setRawData(resMeas(i).res(1).rawData);
%         resRecalc(i,j) = algo.analyse();
%         fittedChAvgFFT(i,j,:) = resRecalc(i,j).fittedChAvgFFT(:,posIdx);
%     end
%     fittedChAvgFFTAvg2(i,:) = sqrt(squeeze(mean(fittedChAvgFFT(i,:,:),2)))-1;
%     
    depthAx = algoVars.len.depthVec*1e3;
    fBar    = algoVars.freq.fBar*1e-6;
    fIdx    = algoVars.freq.fUSIdx;
%     envUp{i} =  envelope(resRecalc(i).fittedFFT(fIdx-env:fIdx+env, 102), 20);
    envUp{i} =  envelope(fittedChAvgFFTAvg2(i,fIdx-env:fIdx+env), 100, 'peaks');
    envBar   =  fBar(fIdx-env:fIdx+env);
    
    fwhmVec(i) = fwhm(envBar*1e6, envUp{i});
    
%     plot(ax(1), depthAx, resRecalc(i,1).phi); hold(ax(1), 'on')
%     plot(ax(2), fBar,    resRecalc(i,1).fittedFFT(:, posIdx)); hold(ax(2), 'on')
    plot(ax(3), envBar,  envUp{i}); hold(ax(3), 'on')
end

% legend(ax(1), legStr)
% legend(ax(2), legStr)
legend(ax(3), legStr)
xlabel(ax(3), "Frequency [MHz]")
ylabel(ax(3), "Amplitude [AU]")

%% Process
algo.setRawData(resMeas.res(1).rawData)
resOrig = algo.analyse();




%%
figure();
subplot(1,2,1)
plot(resVars.measVars.algo.len.depthVec, resOrig.phiNorm)
subplot(1,2,2)
plot(algoVars.freq.fBar*1e-6, resOrig.unFittedAvgFFTStd(:,102))

resNew = algo.analyseEMD('Pos');

figure();
plot(resNew(1).phiNorm)

figure();
for i = 1:10
    ax(i) = subplot(2,5,i);
    plot(resVars.measVars.algo.len.depthVec, resNew(i).phiNorm)
end
linkaxes(ax)

figure();
for i = 1:10
    ax(i) = subplot(2,5,i);
    plot(algoVars.freq.fBar*1e-6, resNew(i).unFittedAvgFFT(:,102)); hold on
end
linkaxes(ax)

figure();
plot(resVars.measVars.algo.len.depthVec, resOrig.phiNorm, '-o'); hold on
plot(resVars.measVars.algo.len.depthVec,resNew(1).phiNorm, '-+')

