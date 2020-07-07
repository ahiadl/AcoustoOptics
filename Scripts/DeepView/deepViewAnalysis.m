depthIdx = 29;
timeVec = dv.time.timeVec;
logPhi = log10(sqrt(res.phi(:,:,depthIdx)));

xVec = dv.s2DVars.grid.firstZero;
dx   = xVec(2)-xVec(1);

phi_sq   = res.phi(:, :, depthIdx);
phi      = sqrt(phi_sq);
grad_phi = gradient(phi);
muEff    = grad_phi./(dx*phi);

avg_muEff = mean(grad_phi./(dx*phi));
std_mueff = std(grad_phi./(dx*phi));

figure();
for i=1:10
    plot(xVec, res.phi(:,i, depthIdx));
    hold on
    legStr(i) = sprintf("T =%d", dv.time.timeVec(i));
end
legend(legStr);
xlabel("X [mm");
ylabel("\phi^2")

figure();
for i=1:10
    plot(xVec, phi(:,i));
    hold on
    legStr(i) = sprintf("T =%d", dv.time.timeVec(i));
end
legend(legStr);
xlabel("X [mm");
ylabel("\phi")

figure();
for i=1:10
    plot(xVec, db(phi(:,i)));
    hold on
    legStr(i) = sprintf("T =%d", dv.time.timeVec(i));
end
legend(legStr);
xlabel("X [mm");
ylabel("\phi [db]")

stdTail = std(phi(14:end,:), 0, 1);
[peakPower, maxIdx] = max(phi, [], 1);

SNR = peakPower./stdTail;

figure();
plot(timeVec, stdTail);
xlabel("AO Duration")
ylabel("Background STD")

figure();
plot(timeVec, SNR);
xlabel("AO Duration")
ylabel("SNR")

figure();
for i=1:5
    plot(dv.s2DVars.grid.firstZero, muEff(:,i));
    hold on
    legStr(i) = sprintf("T =%d", dv.time.timeVec(i));
end
legend(legStr)
xlabel("X [mm");
ylabel("\phi [db]")

%% Load saved Data
resAO  = cell(1,10);
varsAO = cell(1,10);
xVec = dv.s2DVars.grid.firstVec;
for i=1:9:10
    rawDataPath = sprintf("D:/Results/13-May-2020 16-17-01-1Rep DV/TimeResults/T=%d/AOResults/AO-R-1-X-%.2f-Results.mat", timeVec(i), xVec(maxIdx(i)));
    varsPath    = sprintf("D:/Results/13-May-2020 16-17-01-1Rep DV/TimeResults/T=%d/AO-Vars.mat", timeVec(i));
    resAO{i}  = load(rawDataPath);
    tmpVars   = load(varsPath);
    
    tmpVars.uVars.fileSystem.saveRawData        = false;
    tmpVars.uVars.fileSystem.saveNetSignal      = false;
    tmpVars.uVars.fileSystem.saveDemultiplexed  = false;
    tmpVars.uVars.fileSystem.saveReshapedSignal = false;
    tmpVars.uVars.fileSystem.saveFFT            = false;
    tmpVars.uVars.fileSystem.savePhiChCmplx     = false;

    tmpVars.uVars.fileSystem.saveResults = false;
    tmpVars.uVars.fileSystem.saveFigs    = false;

    varsAO{i} = tmpVars;
end
%% Analyze saved Data
ao = acoustoOptics();

i = 10; %select time frame
rawData = ao.loadRawDataToAO(resAO{i});
tmpVars.uVars.fileSystem.extProject = false;
varsAO{i}.uVars.figs.zIdx  = 29;
varsAO{i}.uVars.figs.quant = 500;
varsAO{i}.uVars.figs.displayFullFFT      = true;
varsAO{i}.uVars.figs.validStruct.FFT     = true;
varsAO{i}.uVars.figs.validStruct.qAvgFFT = true;
ao.setMeasVars(varsAO{i}.uVars);
loadedDataRes  = ao.analyseLoadedData();

zIdx  = 29;
quant = 50;
figure()
for j=1:4
    plot(ao.graphics.fBarShift, abs(squeeze(loadedDataRes.fftRes(quant, j, :, zIdx)))); hold on;
end
legend("ch1", "ch2", "ch3", "ch4")
title(sprintf("T = %d, Quant = %d, Z = %d", timeVec(i), quant, zIdx))



figure()
for j=1:4
    plot(ao.graphics.fBarShift, fftshift(fftRes,2)); hold on;
end
legend("ch1", "ch2", "ch3", "ch4")
title(sprintf("Quant averaged: T = %d, Z = %d", timeVec(i), zIdx))

sig_temp=0;
for j=1:4
   sig_temp=sig_temp+(squeeze(loadedDataRes.qAvgFFT(1, j, :, zIdx))).^2;
end
%%
fftRes = abs(squeeze(loadedDataRes.qAvgFFT(1, :, :, zIdx)));
fftResAvg = mean(fftRes,1);
fftResAvgShift = fftshift(fftResAvg);

envDC = 100e3;
envUS = 30e3;
 
samplesInFreq = length(fBar);
fBar = ao.measVars.algo.freq.frequencyBar;
fBarShift = ao.measVars.algo.freq.frequencyBarShifted;
fUS  = ao.measVars.algo.usSignal.fSin;
fs   = ao.measVars.algo.digitizer.fs;

fIdxPos = fUS*samplesInFreq/fs+1;
fIdxNeg = (fUS+(fs/2))*samplesInFreq/fs+1;
DcIdx   = 1;

figure();
plot(fBar, fftResAvg); hold on;
plot(fBar(fIdxPos), fftResAvg(fIdxPos), 'g+')
plot(fBar(fIdxNeg), fftResAvg(fIdxNeg), 'g+')

% if mod(samplesInFreq,2)
%     midComp = floor(samplesInFreqm/2)+1;
% else
%     midComp = floor(samplesInFreqm/2);
% end

shiftFactor = floor(samplesInFreq/2);

fIdxPosShift =  fIdxPos + shiftFactor;
fIdxNegShift =  fIdxNeg - shiftFactor;
DcIdxShift   =  DcIdx + shiftFactor;

envUsIdx = sum(abs(fBar-fUS) < envUS);
envDcIdx = sum(abs(fBar) < envDC);

fitIdxShift = [ 1:(fIdxNegShift - envUsIdx),...
              (fIdxNegShift+envUsIdx):(DcIdxShift-envDcIdx),...
              (DcIdxShift+envDcIdx):(fIdxPosShift-envUsIdx), ...
              (fIdxPosShift+envUsIdx) : samplesInFreq];

fitIdx     = [ envDcIdx:(fIdxPos-envUsIdx),...
              (fIdxPos+envUsIdx):(fIdxNeg-envUsIdx),...
              (fIdxNeg+envUsIdx):(samplesInFreq-envDcIdx)];
          
linearft = fittype({'1',...
                    'x',...
                    '(1/2)*(3*x^2 - 1)',...
                    '(1/2)*(5*x^3 - 3*x)',...
                    '(1/8)*(35*x^4 - 30*x^2 + 3)',...
                    '(1/8)*(63*x^5 - 70*x^3 + 15*x)',...
                    '(1/16)*(231*x^6 - 315*x^4 + 105*x^2 - 5)',...
                    '(1/16)*(429*x^7 - 693*x^5 + 315*x^3 - 35*x)',...
                    '(1/128)*(6435*x^8 - 12012*x^6 + 6930*x^4 - 1260*x^2 + 35)',...
                    '(1/128)*(12155*x^9 - 25740*x^7 + 18018*x^5 - 4620*x^3 + 315*x)',...
                    '(1/256)*(46189*x^10 - 109395*x^8 + 90090*x^6 - 30030*x^4 + 3465*x^2 - 63)'});
f = fit(fBar(fitIdx)', fftResAvg(fitIdx)', linearft);
fshift = fit(fBarShift(fitIdxShift)', fftResAvgShift(fitIdxShift)', linearft);

fmul = fit(fBar(fitIdx)', fftRes(:,fitIdx)', linearft);
fittedData = f.a + ...
             f.b*fBar + ...
             f.c*(1/2)*(3*fBar.^2-1) + ...
             f.d*(1/2)*(5*fBar.^3-3*fBar) + ...
             f.e*(1/8)*(35*fBar.^4-30*fBar.^2+3) + ...
             f.f*(1/8)*(63*fBar.^5-70*fBar.^3+15*fBar)+ ...
             f.g*(1/16)*(231*fBar.^6-315*fBar.^4+105*fBar.^2-5) + ...
             f.h*(1/16)*(429*fBar.^7-693*fBar.^5+315*fBar.^3-35*fBar) + ...
             f.k*(1/128)*(6435*fBar.^8 - 12012*fBar.^6 + 6930*fBar.^4 - 1260*fBar.^2 + 35) + ...
             f.l*(1/128)*(12155*fBar.^9 - 25740*fBar.^7 + 18018*fBar.^5 - 4620*fBar.^3 + 315*fBar) + ...
             f.m*(1/256)*(46189*fBar.^10 - 109395*fBar.^8 + 90090*fBar.^6 - 30030*fBar.^4 + 3465*fBar.^2 - 63);
         
fittedDataShift = fshift.a + ...
                  fshift.b*fBarShift + ...
                  fshift.c*(1/2)*(3*fBarShift.^2-1) + ...
                  fshift.d*(1/2)*(5*fBarShift.^3-3*fBarShift) + ...
                  fshift.e*(1/8)*(35*fBarShift.^4-30*fBarShift.^2+3) + ...
                  fshift.f*(1/8)*(63*fBarShift.^5-70*fBarShift.^3+15*fBarShift)+ ...
                  fshift.g*(1/16)*(231*fBarShift.^6-315*fBarShift.^4+105*fBarShift.^2-5) + ...
                  fshift.h*(1/16)*(429*fBarShift.^7-693*fBarShift.^5+315*fBarShift.^3-35*fBarShift) + ...
                  fshift.k*(1/128)*(6435*fBarShift.^8 - 12012*fBarShift.^6 + 6930*fBarShift.^4 - 1260*fBarShift.^2 + 35) + ...
                  fshift.l*(1/128)*(12155*fBarShift.^9 - 25740*fBarShift.^7 + 18018*fBarShift.^5 - 4620*fBarShift.^3 + 315*fBarShift) + ...
                  fshift.m*(1/256)*(46189*fBarShift.^10 - 109395*fBarShift.^8 + 90090*fBarShift.^6 - 30030*fBarShift.^4 + 3465*fBarShift.^2 - 63);         

figure()
subplot(1,2,1)
plot(fBar, fftResAvg); hold on;
plot(fBar, fittedData);
subplot(1,2,2)
plot(fBar, fftResAvg - fittedData)

figure()
subplot(1,2,1)
plot(fBarShift, fftResAvgShift); hold on;
plot(fBarShift, fittedDataShift);
subplot(1,2,2)
plot(fBarShift, fftResAvgShift - fittedDataShift)
%%
sig = sqrt(sig_temp);

Sig = fft(sig);

figure; 
plot(abs(Sig))

Sig(1:30)=0;
Sig(end-30:end) = 0;
sig2 = ifft(Sig);

figure; 
plot(ao.graphics.fBarShift, abs(sig2))
