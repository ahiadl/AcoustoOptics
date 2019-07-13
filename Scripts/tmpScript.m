A = load('D:\Results\Results-09-Jul-2019 15-46-52-2DScan-Quants-0.05-t-2\Results.mat');
B = load('D:\Results\Results-09-Jul-2019 15-46-52-2DScan-Quants-0.05-t-2\Vars.mat');

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
A = load('D:\Results\Results-11-Jul-2019 01-17-59-Consistency-Laser-2s-rawDara\Results.mat');
B = load('D:\Results\Results-11-Jul-2019 01-17-59-Consistency-Laser-2s-rawDara\Vars.mat');

figure()
yyaxis left
errorbar(B.curVars.scan.timeFrames, A.res.phiFrame(8,:), A.res.phiFrameStd(9,:));hold on
yyaxis right
errorbar(B.curVars.scan.timeFrames, A.res.phiFrame(1,:), A.res.phiFrameStd(9,:));
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
%%
%Find speckle decorrelation time: Mean of fft and fit to lorentzian
A = load('D:\Results\Results-11-Jul-2019 01-17-59-Consistency-Laser-2s-rawDara\rawData\F20.00S1-rawData.mat');
B = load('D:\Results\Results-11-Jul-2019 01-17-59-Consistency-Laser-2s-rawDara\rawData\F20.00S1-Vars.mat');

ao = acoustoOptics();
uVars = ao.uVarsCreate();
uVars = B.curVars.acoustoOptics.uVars;
uVars.quantTime = 0.5;
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
uVars.useHadamard = false;

ao.setMeasVars(uVars);
algoVars = ao.getAlgoVars;

ao.algo.data = A.res;
ao.algo.analyse();

C = ao.algo.res.fftRes;

C = mean(abs(C), 1);

fIdx = algoVars.freq.fSinIdx;
fBar = algoVars.freq.frequencyBar*1e-6;

figure()
plot(abs(squeeze(C(1,1,:,9))))
xlim([fIdx+[-100 +100]]);

idxs = fIdx+(-25:25); 
data = squeeze(C(1,1,idxs,9));
data = data - min(data);
[YPRIME, PARAMS, RESNORM, RESIDUAL,] = lorentzfit(fBar(idxs)*1e6, double(data'));
Gamma(1) = PARAMS(1)*2;
data = squeeze(C(1,2,idxs,9));
data = data - min(data);
[YPRIME, PARAMS, RESNORM, RESIDUAL,] = lorentzfit(fBar(idxs)*1e6, double(data'));
Gamma(2) = PARAMS(1)*2;
data = squeeze(C(1,3,idxs,9));
data = data - min(data);
[YPRIME, PARAMS, RESNORM, RESIDUAL,] = lorentzfit(fBar(idxs)*1e6, double(data'));
Gamma(3) = PARAMS(1)*2;
data = squeeze(C(1,4,idxs,9));
data = data - min(data);
[YPRIME, PARAMS, RESNORM, RESIDUAL,] = lorentzfit(fBar(idxs)*1e6, double(data'));
Gamma(4) = PARAMS(1)*2;

IdealQuantTime = 1/mean(Gamma);

figure(); 
subplot(1,2,1)
plot(fBar(idxs), data, 'x'); hold on
plot(fBar(idxs), YPRIME, 'r-')
title("Fit Average FFT to Lorentzian (Ch 4)")
xlabel('f [MHz]')
ylabel('Amp[v]')

subplot(1,2,2)
stem(Gamma)
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
        uVars.quantTime = IdealQuantTime;
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
errorbar(2:2:20, results.phiFrame(9, 2:2:20), results.phiFrameStd(9, 2:2:20)); hold on;
errorbar(2:2:20, results.phiFrame(9, 2:2:20), results.phiFrameStd(9, 2:2:20));
legend('before Quant opt', 'after Quant opt')

SNRIdeal =  results.phiFrame(9, 2:2:20).^2 ./ results.phiFrameStd(9, 2:2:20).^2;

figure();
plot(2:2:20, SNRIdeal, '-+')
plot(2:2:20, SNRIdeal, '-+');
xlabel ('Frame time [s]')
ylabel ('SNR')
