% In this file:
% first banana shape
% consistency fit analysis

A = load('D:\ResultsToKeep\Scan-Quants-0.05-t-2\Results.mat');
B = load('D:\ResultsToKeep\Scan-Quants-0.05-t-2\Vars.mat');

figure();
subplot(2,1,1)
imagesc('XData', B.curVars.stages.yVec,...
        'YData', B.curVars.stages.xVec,...
        'CData', A.res.phiRep(:,:,9));
colorbar;
title('Point Value (Quants Mean)');
axis equal
axis tight

subplot(2,1,2)
imagesc('XData', B.curVars.stages.yVec,...
        'YData', B.curVars.stages.xVec,...
        'CData', A.res.phiRepStd(:,:,9));
colorbar;
title('Point Std (Quants Std)');
axis equal
axis tight


%% show naive results
A = load('D:\ResultsToKeep\Consistency-Laser-2s-rawData\Results.mat');
B = load('D:\ResultsToKeep\Consistency-Laser-2s-rawData\Vars.mat');

figure()
% yyaxis left
errorbar(B.curVars.scan.timeFrames, A.res.phiFrame(8,:), A.res.phiFrameStd(9,:));hold on
% yyaxis right
errorbar(B.curVars.scan.tmeFrames, A.res.phiFrame(1,:), A.res.phiFrameStd(9,:));
title('Consistency Results');
legend('Strongest Point', 'Weakest Point');
xlabel('Frame Duration [s]');
ylabel('Mean Value');
xlim([0, 22])

SNR = (A.res.phiFrame.^2)./(A.res.phiFrameStd.^2);

figure()
subplot(1,2,1)
plot(B.curVars.scan.timeFrames, SNR(8,:)); hold on
plot(B.curVars.scan.timeFrames, SNR(1,:))
title('SNR vs Time Frame')
legend('Strongest Point' , 'Weakest Point');
xlabel ('Frame Duration [s]')
ylabel ('SNR')

subplot(1,2,2)
plot(B.curVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, SNR(:,1), '-+'); hold on;
plot(B.curVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, SNR(:,4), '-+'); hold on;
plot(B.curVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, SNR(:,8), '-+'); hold on;
plot(B.curVars.acoustoOptics.algoVars.len.zVecUSRes*1e3, SNR(:,10), '-+'); hold on;
legend('2[s]', '8[s]', '16[s]', '20[s]'); 
title('SNR vs Depth (Several time Frames)')
xlabel ('Z [mm]')
ylabel ('SNR')

figure();
subplot(2,2,1)
stem(squeeze(A.res.phiQuant(8,1,1,:)))
title('Time Frame = 2[s]')
subplot(2,2,2)
stem(squeeze(A.res.phiQuant(8,4,1,:)))
title('Time Frame = 8[s]')
subplot(2,2,3)
stem(squeeze(A.res.phiQuant(8,8,1,:)))
title('Time Frame = 16[s]')
subplot(2,2,4)
stem(squeeze(A.res.phiQuant(8,10,1,:)))
title('Time Frame = 20[s]')

figure()
subplot(2,2,1)
stem(squeeze(A.res.phiQuant(8,10,1,:)))
title('Set 1')
subplot(2,2,2)
stem(squeeze(A.res.phiQuant(8,10,5,:)))
title('Set 5')
subplot(2,2,3)
stem(squeeze(A.res.phiQuant(8,10,10,:)))
title('Set 10')
subplot(2,2,4)
stem(squeeze(A.res.phiQuant(8,10,20,:)))
title('Set 20')

naiveData = A;
naiveVars = B;
%%
%Find speckle decorrelation time: Mean of fft and fit to lorentzian
A = load('D:\ResultsToKeep\Consistency-Laser-2s-rawData\rawData\F20.00S2-rawData.mat');
B = load('D:\ResultsToKeep\Consistency-Laser-2s-rawData\rawData\F20.00S1-Vars.mat');

ao = acoustoOptics();
uVars = ao.uVarsCreate();
uVars = B.curVars.acoustoOptics.uVars;
uVars.quantTime = 0.5;
uVars.useQuant = false;
uVars.exportRawData = [];
uVars.exportRawData.rawData       = false;
uVars.exportRawData.netSignal     = false;
uVars.exportRawData.deMultiplexed = false;
uVars.exportRawData.reshape       = false;
uVars.exportRawData.FFT           = true;

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
uVars.useHadamard = false;

uVars.fs.saveFullData = false;
uVars.fs.saveFigs     = false;
uVars.fs.saveResults  = false;
uVars.fs.resDirPath   = false;
uVars.fs.scanName     = false;
uVars.fs.saveAny      = false;

ao.setMeasVars(uVars);
algoVars = ao.getAlgoVars;

ao.algo.data = A.res;
ao.algo.analyse();

C = ao.algo.res.fftRes;

C = mean(abs(C), 1);

fIdx = algoVars.freq.fSinIdx;
fBar = algoVars.freq.frequencyBar*1e-6;

figure()
plot(fBar, abs(squeeze(mean(C(1,:,:,8),2))))
xlim(fBar(fIdx+[-400 +400]));

env = [5, 10, 15, 20, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 50, 60, 70, 85, 100];
figure()
ax(1) = subplot(2,2,1);
ax(2) = subplot(2,2,2);
ax(3) = subplot(2,2,3);
ax(4) = subplot(2,2,4);

figure()
axS = axes();

for i = 1:length(env)
%     figure()
    idxs = fIdx+(-env(i):env(i));
    residual = 0;
    
    data = squeeze(C(1,1,idxs,8));
    data = data - min(data);
    [YPRIME, PARAMS, RESNORM, RESIDUAL,] = lorentzfit(fBar(idxs)*1e6, double(data'));
    Gamma(i,1) = PARAMS(1)*2;
    residual(i,1) = sum(RESIDUAL);
    plot(ax(1), fBar(idxs), data, 'x'); hold(ax(1), 'on')
    plot(ax(1), fBar(idxs), YPRIME, 'r-');hold(ax(1), 'off')
    title(ax(1), num2str(env(i)));

    data = squeeze(C(1,2,idxs,8));
    data = data - min(data);
    [YPRIME, PARAMS, RESNORM, RESIDUAL,] = lorentzfit(fBar(idxs)*1e6, double(data'));
    Gamma(i,2) = PARAMS(1)*2;
    residual(i,2) = mean(RESIDUAL);
    plot(ax(2), fBar(idxs), data, 'x');hold(ax(2), 'on')
    plot(ax(2), fBar(idxs), YPRIME, 'r-'); hold(ax(2), 'off')

    data = squeeze(C(1,3,idxs,8));
    data = data - min(data);
    [YPRIME, PARAMS, RESNORM, RESIDUAL,] = lorentzfit(fBar(idxs)*1e6, double(data'));
    Gamma(i,3) = PARAMS(1)*2;
    residual(i,3) =  mean(RESIDUAL);
    plot(ax(3), fBar(idxs), data, 'x'); hold(ax(3), 'on')
    plot(ax(3), fBar(idxs), YPRIME, 'r-'); hold(ax(3), 'off')

    data = squeeze(C(1,4,idxs,8));
    data = data - min(data);
    [YPRIME, PARAMS, RESNORM, RESIDUAL,] = lorentzfit(fBar(idxs)*1e6, double(data'));
    Gamma(i,4) = PARAMS(1)*2;
    residual(i,4) = mean(RESIDUAL);
    plot(ax(4), fBar(idxs), data, 'x');hold(ax(4), 'on')
    plot(ax(4), fBar(idxs), YPRIME, 'r-'); hold(ax(4), 'off')

%     figure()
    stem(axS, residual(i,:))
    xlim(axS,[0,5])
    title(num2str(env(i)));

    meanRes(i) = mean(residual(i,:));
    stdRes(i) = std(residual(i,:));

    IdealQuantTimeVec(i) = 1/mean(Gamma(i,:))
    
    pause(1)
end

figure();
subplot(1,2,1)
errorbar(env, log(abs(meanRes)), log(stdRes))
title('Mean Residue Value (log)')
xlim([20, 50])
subplot(1,2,2)
stem(env, IdealQuantTimeVec)
title('Quant time vs. Fit Enviroment Size')

[~, idealEnvIdx] = min( stdRes);

idealQuantTime = IdealQuantTimeVec(idealEnvIdx);

idxs = fIdx+(-env(idealEnvIdx):env(idealEnvIdx));

figure();
subplot(2,2,1)
data = squeeze(C(1,1,idxs,8));
data = data - min(data);
[YPRIME, PARAMS, RESNORM, RESIDUAL,] = lorentzfit(fBar(idxs)*1e6, double(data'));
plot(fBar(idxs), data, 'x'); hold on
plot(fBar(idxs), YPRIME, 'r-');hold off
title('ch1')

subplot(2,2,2)
data = squeeze(C(1,2,idxs,8));
data = data - min(data);
[YPRIME, PARAMS, RESNORM, RESIDUAL,] = lorentzfit(fBar(idxs)*1e6, double(data'));
plot(fBar(idxs), data, 'x'); hold on
plot(fBar(idxs), YPRIME, 'r-');hold off
title('ch2')

subplot(2,2,3)
data = squeeze(C(1,3,idxs,8));
data = data - min(data);
[YPRIME, PARAMS, RESNORM, RESIDUAL,] = lorentzfit(fBar(idxs)*1e6, double(data'));
plot(fBar(idxs), data, 'x'); hold on
plot(fBar(idxs), YPRIME, 'r-'); hold off
title('ch3')

subplot(2,2,4)
data = squeeze(C(1,4,idxs,8));
data = data - min(data);
[YPRIME, PARAMS, RESNORM, RESIDUAL,] = lorentzfit(fBar(idxs)*1e6, double(data'));
plot(fBar(idxs), data, 'x'); hold on
plot(fBar(idxs), YPRIME, 'r-');  hold off
title('ch4')

figure(); 
subplot(1,2,1)
plot(fBar(idxs), data, 'x'); hold on
plot(fBar(idxs), YPRIME, 'r-')
title("Fit Average FFT to Lorentzian (Ch 4)")
xlabel('f [MHz]')
ylabel('Amp[v]')

subplot(1,2,2)
stem(Gamma(idealEnvIdx, :))
title('Gamma Vs. Channels')
xlabel('Channel')
ylabel('Gamma')
xlim([0, 5])

%% Analyse results with Ideal quant time
framesVec = 2:2:20;
dirName = 'D:\Results\Results-11-Jul-2019 01-17-59-Consistency-Laser-2s-rawDara\rawData\';
for i = 1:10
    for j =1:20
        disp([i, j])
        dataName = sprintf("F%.2fS%d-rawData.mat", framesVec(i), j );
        varsName = sprintf("F%.2fS%d-Vars.mat", framesVec(i), j );
        
        A = load(sprintf("%s%s", dirName, dataName));
        B = load(sprintf("%s%s", dirName, varsName));
        
        uVars = B.curVars.acoustoOptics.uVars;
        uVars.quantTime = idealQuantTime;
        uVars.exportRawData = true;
        uVars.useHadamard = false;
        uVars.gReq = gReq;
        
        if (j ==1)
            ao.setMeasVars(uVars);
            algoVars = ao.getAlgoVars;
        end
        
        k = algoVars.samples.numOfQuant;
        quantLenVec(i) = k;
        ao.algo.data = A.res;
        res = ao.algo.analyse();
        
        results.phiCh(:,i,j,1:k,:)  = permute(gather(res.phiCh), [3,4,5,6,1,2]);
        results.phiQuant(:,i,j,1:k) = permute(gather(res.phiQuant), [2,3,4,5,1]);
        results.phiSets(:,i,j)      = gather(res.phi);
        results.phiSetsStd(:,i,j)   = gather(res.phiStd);
        
    end
    results.phiFrame(:,i)    = mean(results.phiSets(:,i,:), 3);
    results.phiFrameStd(:,i) = std(results.phiSets(:,i,:), 0, 3);
end


figure(); 
errorbar(naiveVars.curVars.scan.timeFrames, naiveData.res.phiFrame(8,:), naiveData.res.phiFrameStd(8,:));hold on;
errorbar(2:2:20, results.phiFrame(8, :), results.phiFrameStd(8, :));
legend('before Quant opt', 'after Quant opt')
xlim([0, 22])
xlabel('Frame Duration[s]')
ylabel ('Amp[v]');
title('Consistency Results - Comparison')

SNRIdeal =  results.phiFrame.^2 ./ results.phiFrameStd.^2;

figure();
plot(2:2:20, log(SNR(8,:)), '-+'); hold on ;
plot(2:2:20, log(SNRIdeal(8, :)), '-+');
xlabel ('Frame time [s]')
ylabel ('SNR')
legend('before quant opti', 'after quant opt')
title('SNR Comparison (log Scale)')
xlim([0, 22])

figure();
stem(2:2:20, results.phiFrameStd(8, :))
title('Std Trendline')
xlabel('Frame Duration')
ylabel('STD')
xlim([0, 22])