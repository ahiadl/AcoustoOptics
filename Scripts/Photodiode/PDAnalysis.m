%% -------------------
% 19-08-21
% --------------------
close all
analysis = Analysis();

usIdxLow    = 24;
usIdxHigh   = 56;
scanIdxLow  = 1;
scanIdxHigh = 60;

numElemUS   = usIdxHigh   - usIdxLow +1;
numElemScan = scanIdxHigh - scanIdxLow +1;

sigScanIdx = 18;
sigUsIdx   = 39;

avgNum = 5;

%% PMT
resScanPMT = load('D:\Photodiode\EqualPhotonos-PDProject\PD-Scan\19-Aug-2021 17-19-02-PMT-scan2D-600V\2DScan-Results.mat');
varsScanPMT = load('D:\Photodiode\EqualPhotonos-PDProject\PD-Scan\19-Aug-2021 17-19-02-PMT-scan2D-600V\2DScan-Vars.mat');

pmtMeasIdx = 1;

depthAxCntr = varsScanPMT.grid.depthCntr*1e3;
depthAxZero = varsScanPMT.grid.depthZero*1e3;
scanAx      = varsScanPMT.grid.scanZero;

depthAxCntrCrop = depthAxCntr(1:numElemScan);
depthAxZeroCrop = depthAxCntrCrop - mean(depthAxCntrCrop);
scanAxCrop  = scanAx(1:numElemScan);

PMTMeas       = squeeze(resScanPMT.phiNorm);
PMTMeasLog    = squeeze(resScanPMT.phiLog);
PMTMeasAvg    = mean(PMTMeas(:,:, 1:avgNum), 3);
PMTMeasAVGLog = db(PMTMeasAvg);

pmtSigScan       = analysis.norm(PMTMeas(:, sigUsIdx, pmtMeasIdx));
pmtSigScanLog    = db(pmtSigScan);
pmtSigScanAvg    = analysis.norm(PMTMeasAvg(:, sigUsIdx));
pmtSigScanLogAvg = db(pmtSigScanAvg);

pmtSigUs       = analysis.norm(PMTMeas(sigScanIdx, :,pmtMeasIdx));
pmtSigUsLog    = db(pmtSigUs);
pmtSigUsAvg    = analysis.norm(PMTMeasAvg(sigScanIdx, :));
pmtSigUsLogAvg = db(pmtSigUsAvg);

PMTim    = analysis.norm2D(PMTMeas(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh, pmtMeasIdx));
PMTimLog = PMTMeasLog(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh, pmtMeasIdx);
PMTimAvg = analysis.norm2D(PMTMeasAvg(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh));
PMTimAvgLog = db(PMTimAvg);

figure(1);
subplot(2,2,1)
imagesc(depthAxZeroCrop, scanAxCrop, PMTim);
axis equal tight
title("PMT Measurement")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar
subplot(2,2,2)
imagesc(depthAxZeroCrop, scanAxCrop,PMTimLog);
axis equal tight
title("PMT Measurement (db)")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar
subplot(2,2,3)
imagesc(depthAxZeroCrop, scanAxCrop,PMTimAvg);
axis equal tight
title("Average PMT Measurement ")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar
subplot(2,2,4)
imagesc(depthAxZeroCrop, scanAxCrop,PMTimAvgLog);
axis equal tight
title("Average PMT Measurement (db)")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar

figure(2);
subplot(2,2,1)
plot(scanAx, pmtSigScan);
title("PMT Penetration Scan")
xlabel("Scan [mm]")
ylabel("Fluence [AU]")
subplot(2,2,2)
plot(scanAx,pmtSigScanLog);
title("PMT Penetration Scan (db)")
xlabel("Scan [mm]")
ylabel("Fluence [dB]")
subplot(2,2,3)
plot(scanAx,pmtSigScanAvg);
title("Average PMT Penetration Scan")
xlabel("Scan [mm]")
ylabel("Fluence [AU]")
subplot(2,2,4)
plot(scanAx,pmtSigScanLogAvg);
title("Average PMT Penetration Scan (db)")
xlabel("Scan [mm]")
ylabel("Fluence [dB]")

figure(3);
subplot(2,2,1)
plot(depthAxZero,pmtSigUs);
title("PMT Ultrasound Axis")
xlabel("US [mm]")
ylabel("Fluence [AU]")
subplot(2,2,2)
plot(depthAxZero, pmtSigUsLog);
title("PMT Ultrasound Axis (db)")
xlabel("US [mm]")
ylabel("Fluence [dB]")
subplot(2,2,3)
plot(depthAxZero, pmtSigUsAvg);
title("Averaged PMT Ultrasound Axis")
xlabel("US [mm]")
ylabel("Fluence [AU]")
subplot(2,2,4)
plot(depthAxZero, pmtSigUsLogAvg);
title("Averaged PMT Ultrasound Axis (db)")
xlabel("US [mm]")
ylabel("Fluence [dB]")

%% PD
resScanPD = load('D:\Photodiode\EqualPhotonos-PDProject\PD-Scan\19-Aug-2021 18-17-18-PD-scan2D-600V\2DScan-Results.mat');
varsScanPD = load('D:\Photodiode\EqualPhotonos-PDProject\PD-Scan\19-Aug-2021 18-17-18-PD-scan2D-600V\2DScan-Vars.mat');

pdMeasIdx = 1;

depthAxCntr = varsScanPD.grid.depthCntr*1e3;
scanAx      = varsScanPD.grid.scanZero;

depthAxCntrCrop = depthAxCntr(1:numElemScan);
depthAxZeroCrop = depthAxCntrCrop - mean(depthAxCntrCrop);
scanAxCrop  = scanAx(1:numElemScan);

PDMeas       = squeeze(resScanPD.phiNorm);
PDMeasLog    = squeeze(resScanPD.phiLog);
PDMeasAvg    = mean(PDMeas(:,:, 1:avgNum), 3);
PDMeasAVGLog = db(PDMeasAvg);

pdSigScan       = analysis.norm(PDMeas(:, sigUsIdx, pdMeasIdx));
pdSigScanLog    = db(pdSigScan);
pdSigScanAvg    = analysis.norm(PDMeasAvg(:, sigUsIdx));
pdSigScanLogAvg = db(pdSigScanAvg);

pdSigUs       = analysis.norm(PDMeas(sigScanIdx, :,pdMeasIdx));
pdSigUsLog    = db(pdSigUs);
pdSigUsAvg    = analysis.norm(PDMeasAvg(sigScanIdx, :));
pdSigUsLogAvg = db(pdSigUsAvg);

PDim    = analysis.norm2D(PDMeas(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh, pdMeasIdx));
PDimLog = PDMeasLog(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh, pdMeasIdx);
PDimAvg = analysis.norm2D(PDMeasAvg(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh));
PDimAvgLog = db(PDimAvg);

figure(4);
subplot(2,2,1)
imagesc(depthAxZeroCrop, scanAxCrop, PDim);
axis equal tight
title("PD Measurement")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar
subplot(2,2,2)
imagesc(depthAxZeroCrop, scanAxCrop,PDimLog);
axis equal tight
title("PD Measurement (db)")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar
subplot(2,2,3)
imagesc(depthAxZeroCrop, scanAxCrop,PDimAvg);
axis equal tight
title("Average PD Measurement ")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar
subplot(2,2,4)
imagesc(depthAxZeroCrop, scanAxCrop,PDimAvgLog);
axis equal tight
title("Average PD Measurement (db)")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar

figure(5);
subplot(2,2,1)
plot(scanAx, pdSigScan);
title("PD Penetration Scan")
xlabel("Scan [mm]")
ylabel("Fluence [AU]")
subplot(2,2,2)
plot(scanAx,pdSigScanLog);
title("PD Penetration Scan (db)")
xlabel("Scan [mm]")
ylabel("Fluence [dB]")
subplot(2,2,3)
plot(scanAx,pdSigScanAvg);
title("Average PD Penetration Scan")
xlabel("Scan [mm]")
ylabel("Fluence [AU]")
subplot(2,2,4)
plot(scanAx,pdSigScanLogAvg);
title("Average PD Penetration Scan (db)")
xlabel("Scan [mm]")
ylabel("Fluence [dB]")

figure(6);
subplot(2,2,1)
plot(depthAxZero,pdSigUs);
title("PD Ultrasound Axis")
xlabel("US [mm]")
ylabel("Fluence [AU]")
subplot(2,2,2)
plot(depthAxZero, pdSigUsLog);
title("PD Ultrasound Axis (db)")
xlabel("US [mm]")
ylabel("Fluence [dB]")
subplot(2,2,3)
plot(depthAxZero, pdSigUsAvg);
title("Averaged PD Ultrasound Axis")
xlabel("US [mm]")
ylabel("Fluence [AU]")
subplot(2,2,4)
plot(depthAxZero, pdSigUsLogAvg);
title("Averaged PD Ultrasound Axis (db)")
xlabel("US [mm]")
ylabel("Fluence [dB]")

%% Combined

minLin = min([min(PMTim(:)), min(PDim(:))]);
maxLin = max([max(PMTim(:)), max(PDim(:))]);

% minVec = mink(PMTimLog(:), 10);
% minVec = mink(PDimLog(:), 10);

minLog = min([min(PMTimLog(:)), min(PDimLog(:))]);
maxLog = max([max(PMTimLog(:)), max(PDimLog(:))]);

climLin = [minLin, maxLin];
climLog = [minLog, maxLog];

figure(7)
subplot(2,2,1)
imagesc(depthAxZeroCrop, scanAxCrop, PMTim, climLin);
axis equal tight
title("PMT Measurement")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar
subplot(2,2,2)
imagesc(depthAxZeroCrop, scanAxCrop, PDim, climLin);
axis equal tight
title("PD Measurement")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar
subplot(2,2,3)
imagesc(depthAxZeroCrop, scanAxCrop, PMTimLog, climLog);
axis equal tight
title("PMT Measurement (db)")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar
subplot(2,2,4)
imagesc(depthAxZeroCrop, scanAxCrop, PDimLog, climLog);
axis equal tight
title("PD Measurement (db)")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar

minLin = min([min(PMTimAvg(:)), min(PDimAvg(:))]);
maxLin = max([max(PMTimAvg(:)), max(PDimAvg(:))]);

minVecPMT = unique(mink(PMTimAvgLog(:), 10));
minVecPD = unique(mink(PDimAvgLog(:), 10));

minLog = min([minVecPMT(2), minVecPD(2)]);
maxLog = max([max(PMTimAvgLog(:)), max(PDimAvgLog(:))]);

climLin = [minLin, maxLin];
climLog = [minLog, maxLog];

figure(8)
subplot(2,2,1)
imagesc(depthAxZeroCrop, scanAxCrop,PMTimAvg, climLin);
axis equal tight
title("Average PMT Measurement ")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar
subplot(2,2,2)
imagesc(depthAxZeroCrop, scanAxCrop,PDimAvg, climLin);
axis equal tight
title("Average PD Measurement ")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar
subplot(2,2,3)
imagesc(depthAxZeroCrop, scanAxCrop, PMTimAvgLog, climLog);
axis equal tight
title("Average PMT Measurement (db)")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar
subplot(2,2,4)
imagesc(depthAxZeroCrop, scanAxCrop, PDimAvgLog, climLog);
axis equal tight
title("Average PD Measurement (db)")
xlabel("US [mm]")
ylabel("Scan [mm]")
colorbar

figure(9)
subplot(2,2,1)
plot(scanAx, pmtSigScan);hold on
plot(scanAx, pdSigScan);
title("Penetration Scan")
xlabel("Scan [mm]")
ylabel("Fluence [AU]")
legend('PMT', 'PD')
subplot(2,2,2)
plot(scanAx,pmtSigScanLog); hold on
plot(scanAx,pdSigScanLog);
title("Penetration Scan (db)")
xlabel("Scan [mm]")
ylabel("Fluence [dB]")
legend('PMT', 'PD')
subplot(2,2,3)
plot(scanAx,pmtSigScanAvg); hold on
plot(scanAx,pdSigScanAvg);
title("Average Penetration Scan")
xlabel("Scan [mm]")
ylabel("Fluence [AU]")
legend('PMT', 'PD')
subplot(2,2,4)
plot(scanAx,pmtSigScanLogAvg); hold on
plot(scanAx,pdSigScanLogAvg);
title("Average Penetration Scan (db)")
xlabel("Scan [mm]")
ylabel("Fluence [dB]")
legend('PMT', 'PD')

figure(10);
subplot(2,2,1)
plot(depthAxZero,pmtSigUs); hold on
plot(depthAxZero,pdSigUs); 
title("PD Ultrasound Axis")
xlabel("US [mm]")
ylabel("Fluence [AU]")
legend('PMT', 'PD')
subplot(2,2,2)
plot(depthAxZero, pmtSigUsLog); hold on
plot(depthAxZero, pdSigUsLog);
title("PD Ultrasound Axis (db)")
xlabel("US [mm]")
ylabel("Fluence [dB]")
legend('PMT', 'PD')
subplot(2,2,3)
plot(depthAxZero, pmtSigUsAvg); hold on
plot(depthAxZero, pdSigUsAvg);
title("Averaged PD Ultrasound Axis")
xlabel("US [mm]")
ylabel("Fluence [AU]")
legend('PMT', 'PD')
subplot(2,2,4)
plot(depthAxZero, pmtSigUsLogAvg); hold on
plot(depthAxZero, pdSigUsLogAvg);
title("Averaged PD Ultrasound Axis (db)")
xlabel("US [mm]")
ylabel("Fluence [dB]")
legend('PMT', 'PD')

%% PMT vs. Gain
% close all
% clear all
clc

analysis = Analysis();

load('D:\Photodiode\EqualPhotonos-PDProject\SNR-Vs-Gain\expVars-noRawData.mat');
idx = 8;
% % resTotal = pmtRes.resTotal;
% expVars = pmtRes.expVars;

aoLen = size(resTotal(1).rawPhi, 2);
reps = size(resTotal(1).rawPhi, 1);
numOfGain = length(expVars.gain);
phiPMT = zeros(numOfGain, reps, aoLen);

for i=1:numOfGain
   phiPMT(i,:,:) = resTotal(i).rawPhi; 
   pmtSNRMat(i,:) = resTotal(i).SNR;
   pmtBkdAvgMat(i,:) = resTotal(i).avgBkg;
   pmtBkgStdMat(i,:) = resTotal(i).stdBkg;
   pmtMaxMat(i,:) = resTotal(i).maxVal;
end

avgMax = mean(pmtMaxMat,2);
avgMaxNorm = avgMax / mean(avgMax);
avgBkgStd = mean(pmtBkgStdMat,2);
avgBkgStdNorm = avgBkgStd./ mean(avgBkgStd);
pmtAvgSNR = avgMaxNorm./avgBkgStdNorm;
pmtAvgSNRNorm = pmtAvgSNR./mean(pmtAvgSNR);

maxMatAC = pmtMaxMat - pmtBkdAvgMat;
voltage = expVars.gain;

m         = 0.007;
n         = 0.4;
log10Gain = voltage*m+n;
gain      = 10.^log10Gain;

pmtSNR        = expVars.SNRAvg;
pmtSNRstd     = expVars.SNRstd;
pmtNormSNR    = expVars.NormSNRAvg;
pmtNormSNRstd = expVars.NormSNRStd;


pmtStdMax     = std(maxMatAC,0,2);
pmtAvgMax     = mean(maxMatAC, 2);
pmtNormAvgMax = maxMatAC./repmat(pmtAvgMax, 1, reps);

pmtStdStdBkg     = std(pmtBkgStdMat,0, 2);
pmtAvgStdBkg     = mean(pmtBkgStdMat,2);
pmtNormAvgStdBkg = pmtBkgStdMat./repmat(pmtAvgStdBkg, 1, reps);

pmtSigStability = pmtAvgMax ./ pmtStdMax;

pmtNormMidSNR = pmtSNRMat./repmat(pmtSNR', 1, reps);


figure(11);
ax = axes();
plot(ax, voltage, avgMaxNorm);hold on
plot(ax, voltage, avgBkgStdNorm)
plot(ax, voltage, pmtAvgSNRNorm)
xlabel('Voltage')
ylabel('Normalized Value')
legend('Signal', 'Noise', 'SNR')
title("PMT Normalized Statistics - Averaged")
% set(ax, 'XScale', 'log')

figure(12);
ax = axes();
errorbar(ax, gain, pmtSNR, pmtSNRstd); hold on
errorbar(ax, gain, pmtNormSNR, pmtNormSNRstd);
set(ax, 'XScale', 'log')
% set(ax, 'YScale', 'log')
title("PMT SNR vs. Gain")

figure(13)
plot(pmtNormAvgMax(idx,:)); hold on
plot(pmtNormAvgStdBkg(idx,:))
plot(pmtNormMidSNR(idx,:))
title("SNR Components - 100 Repetitive Measurements");
legend("Sig Max", "Bkg std", "SNR")

% figure();
% ax = axes();
% plot(gain, pmtSigStability)
% set(ax, 'XScale', 'log')

%% PD: SNR vs. Ref
% clear all
% close all
clc;

load('D:\Photodiode\EqualPhotonos-PDProject\SNR-Vs-Ref\19-Aug-2021 14-23-39-SNR-Vs-Ref-ExpVars.mat')

pdSNR        = expVars.SNRAvg;
pdSNRstd     = expVars.SNRstd;
pdNormSNR    = expVars.NormSNRAvg;
pdNormSNRstd = expVars.NormSNRStd;

aoLen = size(resTotal(1).rawPhi, 2);
reps = size(resTotal(1).rawPhi, 1);
theta = expVars.thetaVec;
numOfTheta = length(theta);

phiPMT = zeros(numOfTheta, reps, aoLen);

for i=1:numOfTheta
   pdPhi(i,:,:)     = resTotal(i).rawPhi; 
   pdSNRMat(i,:)    = resTotal(i).SNR;
   pdBkdAvgMat(i,:) = resTotal(i).avgBkg;
   pdBkgStdMat(i,:) = resTotal(i).stdBkg;
   pdMaxMat(i,:)    = resTotal(i).maxVal;
end

sigIntensity = expVars.netSignal;
minIntensity = expVars.minIntensity;
maxIntensity = expVars.maxIntensity;
Ispan = maxIntensity - minIntensity;

thetaVec = flip(theta)-theta(end);
intensityVec = ((sin(deg2rad(thetaVec))).^2)*Ispan;
refFactor = intensityVec/sigIntensity;

pdAvgMax        = mean(pdMaxMat,2);
pdAvgMaxNorm    = pdAvgMax / mean(pdAvgMax);
pdAvgBkgStd     = mean(pdBkgStdMat,2);
pdAvgBkgStdNorm = pdAvgBkgStd./ mean(pdAvgBkgStd);
pdAvgSNR        = pdAvgMaxNorm./pdAvgBkgStdNorm;
pdAvgSNRNorm    = pdAvgSNR./mean(pdAvgSNR);

h = figure(15); clf(h)
ax = axes(h);
plot(ax, refFactor, pdAvgMaxNorm);hold on
plot(ax, refFactor, pdAvgBkgStdNorm)
plot(ax, refFactor, pdAvgSNRNorm)
set(ax, 'XScale', 'log')
xlabel('Voltage')
ylabel('Normalized Value')
legend('Signal', 'Noise', 'SNR')
title("PD Normalized Statistics - Averaged")


h = figure(16); clf(h)
ax = axes(h);
errorbar(ax, refFactor, pdSNR, pdSNRstd); hold on
errorbar(ax, refFactor, pdNormSNR, pdNormSNRstd);
set(ax, 'XScale', 'log')

%%
% close all
load('D:\Photodiode\EqualPhotonos-PDProject\PD-peakSNR - 5-times\23-Aug-2021 15-17-31-PD-peakSNR - 5-times-ExpVars.mat');

data = expData.data;

SNRmat = data.SNR;
for i = 1:5
    SNRVec(100*(i-1)+1 : i*100) = SNRmat(i,:);
end

avgCumSNR = cumsum(SNRVec)./(1:length(SNRVec));
for i = 1:400
    SNRSubSet(i) =  mean(SNRVec(i:i+100));
end

figure();
plot(SNRVec); hold on
plot(avgCumSNR); 
plot(101:500, SNRSubSet);
legend("SNR", "CumAvgSNR", "Sub Set SNR")


% figure()
% plot(

pd2 = load("D:\Photodiode\EqualPhotonos-PDProject\PD-peakSNR - 5-times\23-Aug-2021 15-05-31-PD-peakSNR - 5-times-22.50.mat");
pd1 = load("D:\Photodiode\EqualPhotonos-PDProject\PD-peakSNR - 5-times\23-Aug-2021 14-59-32-PD-peakSNR - 5-times-22.50.mat");

aoVars = expData.vars.ao;
fBar = aoVars.measVars.algo.freq.frequencyBarShifted;

fft1Sig = pd1.data.res(24).unFittedFFTShift(:,40);
fft1Bkg = pd1.data.res(24).unFittedFFTShift(:,97);
phi1 = pd1.data.res(24).phi;
SNR1 = SNRmat(2,24);
BKG1 = data.bkg(2,24,:);

fft2Sig = pd2.data.res(11).unFittedFFTShift(:,41);
fft2Bkg = pd2.data.res(11).unFittedFFTShift(:,97);
phi2    = pd2.data.res(11).phi;
SNR2= SNRmat(3,11);

phi1 = pd1.data.res(24).phi;
phi2    = pd2.data.res(11).phi;
phi3    = pd1.data.res(57).phi;
phi4    = pd1.data.res(88).phi;

figure();
plot(phi1); hold on
plot(phi2);
plot(phi3);
plot(phi4);

legend("124", "157", "188", "211")

stdVec = data.normStd';
maxVec = data.normMax';
SNRVec = data.normSNR';

stdVec = stdVec(:);
maxVes = maxVec(:);
SNRVec = SNRVec(:);

figure()
plot(SNRVec); hold on
plot(maxVes); 
plot(stdVec);
legend("SNR", "Sig", "Std Bkg")


figure();
subplot(2,2,1)
stem(phi1)
title("SNR = 514")
ylim([0, 1.6e-4])
subplot(2,2,2)
title("SNR = 85.43")
stem(phi2)
ylim([0, 1.6e-4])
subplot(2,2,3)
plot(fBar,fft1Sig); hold on;
plot(fBar, fft1Bkg);
legend("signal", "bkg")
subplot(2,2,4)
plot(fBar,fft2Sig); hold on;
plot(fBar, fft2Bkg);
legend("signal", "bkg")

%% PMT dev vs PD dev

pmtRes = load('D:\Photodiode\EqualPhotonos-PDProject\PMT-20MhzSR-5HzWL-rawPhi-575\24-Aug-2021 19-47-02-PMT-20MhzSR-5HzWL-rawPhi-575-ExpVars.mat');
pdRes = load('D:\Photodiode\EqualPhotonos-PDProject\PD-20MhzSR-5HzWL-rawPhi-coarseScan\24-Aug-2021 19-19-05-PD-20MhzSR-5HzWL-rawPhi-coarseScan-ExpVars.mat');

pmtVars = pmtRes.expData.vars;
pmtData = pmtRes.expData.data;

pdVars = pdRes.expData.vars;
pdData = pdRes.expData.data;


SNRPMT = pmtData.SNRAvgMeas;

SNRFactor = pdData.SNRAvgMeas / SNRPMT;
SNRstdFactor = pdData.SNRstdMeas/SNRPMT;


stdVec = pmtData.normStd';
maxVec = pmtData.normMax';
SNRVec = pmtData.normSNR';

stdVec = stdVec(:);
maxVes = maxVec(:);
SNRVec = SNRVec(:);

figure()
plot(SNRVec); hold on
plot(maxVes); 
plot(stdVec);
legend("SNR", "Sig", "Std Bkg")


figure()
ax = axes();
errorbar(pdVars.pwr.refFactor, SNRFactor, SNRstdFactor)
set(ax, 'XScale', 'log')
xlabel("Reference Factor")
ylabel("SNR Gain")


%%

pdRes = load('D:\Photodiode\EqualPhotonos-PDProject\PD-20MhzSR-5HzWL-rawPhi-NoDonut\25-Aug-2021 14-46-35-PD-20MhzSR-5HzWL-rawPhi-NoDonut-ExpVars.mat');

% pmtVars = pmtRes.expData.vars;
% pmtData = pmtRes.expData.data;

pdVars = pdRes.expData.vars;
pdData = pdRes.expData.data;


SNRPMT = 144;

SNRFactor = pdData.SNRAvgMeas / SNRPMT;
SNRstdFactor = pdData.SNRstdMeas/SNRPMT;

figure()
ax = axes();
errorbar(pdVars.pwr.refFactor, SNRFactor, SNRstdFactor)
set(ax, 'XScale', 'log')
xlabel("Reference Factor")
ylabel("SNR Gain")


%%
resPMT = load("D:\Photodiode\EqualPhotonos-PDProject\PMT-20MhzSR-5HzWL-Polarized-SM\26-Aug-2021 10-05-41-PD-20MhzSR-5HzWL-Polarized-SM-ExpVars");
resPD  = load("D:\Photodiode\EqualPhotonos-PDProject\PD-20MhzSR-5HzWL-Polarized-SM-thetaScan\26-Aug-2021 10-59-27-PD-20MhzSR-5HzWL-Polarized-SM-thetaScan-ExpVars");

pmtSNR = resPMT.expData.data.SNRAvgMeas;
pdSNR  = resPD.expData.data.SNRAvgMeas;
pdSNRstd = resPD.expData.data.SNRstdMeas;

snrGain = pdSNR./pmtSNR;
snrGainStd = pdSNRstd/pmtSNR;

refAx = resPD.expData.vars.pwr.refFactor;
refAx = 1.0e+04 * [0.0020    0.0457    0.2304    0.5382    0.9319    1.3642    1.7828    2.1373];

figure();
ax = axes();
errorbar(ax, refAx, snrGain, snrGainStd)
set(ax, 'XScale', 'log')

% ref = [50, 100, 250, 500, 1000, 2500, 5000, 7500, 9000, 10000, 11000, 12500, 15000, 17500]


%%

cwPMT = load("D:\Photodiode\InStabilityCheck\PMT-2s-CW\01-Sep-2021 09-02-15-PMT-2s-CW-ExpVars.mat");

cwPMTData = cwPMT.expData.data;

stdVec = cwPMTData.normStd';
maxVec = cwPMTData.normMax';
SNRVec = cwPMTData.normSNR';

stdVec = stdVec(:);
maxVes = maxVec(:);
SNRVec = SNRVec(:);

figure()
% plot(SNRVec); hold on
plot(maxVes); 
% plot(stdVec);
% legend("SNR", "Sig", "Std Bkg")


%%
pdRes = load("D:\Photodiode\InStabilityCheck\PD-20MhzSR-0.5HzWL-Polarized-SM-withUltrasound-OnSpot-0.002Q-8s\26-Aug-2021 19-19-41-PD-20MhzSR-0.5HzWL-Polarized-SM-withUltrasound-OnSpot-0.002Q-8s-ExpVars.mat");

pdResData = pdRes.expData.data;

stdVec = pdResData.normStd';
maxVec = pdResData.normMax';
SNRVec = pdResData.normSNR';

stdVec = stdVec(:);
maxVes = maxVec(:);
SNRVec = SNRVec(:);

figure()
plot(SNRVec); hold on
plot(maxVes); 
plot(stdVec);
legend("SNR", "Sig", "Std Bkg")


%%

pdRes = load("D:\Photodiode\EqualPhotonos-PDProject\PD-20MhzSR-5HzWL-Polarized-SM-thetaScan-Fine\26-Aug-2021 17-13-22-PD-20MhzSR-5HzWL-Polarized-SM-thetaScan-Fine-ExpVars.mat");

pdData = pdRes.expData.data;

pdRes.expData.vars.uVars.refFactor

%%

pdRes = load("D:\Photodiode\SmallReferenceSpot\PD-2sT-5KhzRR-0.5HzWL-20MHzFs-1msQ-50reps-SNR-vs-Ref-Theta-4mW\02-Sep-2021 16-26-09-PD-2sT-5KhzRR-0.5HzWL-20MHzFs-1msQ-50reps-SNR-vs-Ref-Theta-4mW-ExpVars.mat");


data = pdRes.expData.data;
vars = pdRes.expData.vars;
j=1;
i=1;

numOfRef  = vars.measLen;
numOfReps = vars.reps;

for j = 1:numOfRef
    for i=1:numOfReps
        rawStd(i,j) = pdRes.loopsData(j).res(i).std;
    end
end



fftBeforeSig = pdRes.expData.data.fftResSig(10)
% fft

figure();
plot(rawStd(:))
legend

idx = 19;
stdVec = data.normStd(idx, :);
maxVec = data.normMax(idx, :);
SNRVec = data.normSNR(idx, :);

figure()
plot(SNRVec); hold on
plot(maxVec); 
plot(stdVec);

%%
close all
analysis = Analysis();

usIdxLow    = 24;
usIdxHigh   = 56;

scanIdxLow  = 1;
scanIdxHigh = 73;

numElemUS   = usIdxHigh   - usIdxLow +1;
numElemScan = scanIdxHigh - scanIdxLow +1;

sigScanIdx = 16;
sigUsIdx   = 36;

avgNum = 5;

%% PMT
% resScanPMT = load('D:\Results\05-Sep-2021 12-54-44-Scan2D-PMT\2DScan-Results.mat');
% varsScanPMT = load('D:\Results\05-Sep-2021 12-54-44-Scan2D-PMT\2DScan-Vars.mat');

pmtMeasIdx = 9;

depthAxCntr = varsScanPMT.grid.depthCntr*1e3;
depthAxZero = varsScanPMT.grid.depthZero*1e3;
scanAx      = varsScanPMT.grid.scanZero;

depthAxCntrCrop = depthAxCntr(1:numElemScan);
depthAxZeroCrop = depthAxCntrCrop - mean(depthAxCntrCrop);
scanAxCrop  = scanAx(1:numElemScan);

PMTMeas       = squeeze(resScanPMT.phiNorm);
PMTMeasLog    = squeeze(resScanPMT.phiLog);
PMTMeasAvg    = mean(PMTMeas(:,:, 1:avgNum), 3);
PMTMeasAVGLog = db(PMTMeasAvg);

pmtSigScan       = analysis.norm(PMTMeas(:, sigUsIdx, pmtMeasIdx));
pmtSigScanLog    = db(pmtSigScan);
pmtSigScanAvg    = analysis.norm(PMTMeasAvg(:, sigUsIdx));
pmtSigScanLogAvg = db(pmtSigScanAvg);

pmtSigUs       = analysis.norm(PMTMeas(sigScanIdx, :,pmtMeasIdx));
pmtSigUsLog    = db(pmtSigUs);
pmtSigUsAvg    = analysis.norm(PMTMeasAvg(sigScanIdx, :));
pmtSigUsLogAvg = db(pmtSigUsAvg);

PMTim    = analysis.norm2D(PMTMeas(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh, pmtMeasIdx));
PMTimLog = PMTMeasLog(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh, pmtMeasIdx);
PMTimAvg = analysis.norm2D(PMTMeasAvg(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh));
PMTimAvgLog = db(PMTimAvg);
% PD
% resScanPD = load('D:\Results\05-Sep-2021 10-48-37-Scan2D-PD\2DScan-Results.mat');
% varsScanPD = load('D:\Results\05-Sep-2021 10-48-37-Scan2D-PD\2DScan-Vars.mat');

pdMeasIdx = 3;

depthAxCntr = varsScanPD.grid.depthCntr*1e3;
scanAx      = varsScanPD.grid.scanZero;

depthAxCntrCrop = depthAxCntr(1:numElemScan);
depthAxZeroCrop = depthAxCntrCrop - mean(depthAxCntrCrop);
scanAxCrop  = scanAx(1:numElemScan);

phi = squeeze(resScanPD.rawPhi);
sizeVec = size(phi);
tmp = reshape(phi, sizeVec(1)*sizeVec(2), []);
maxVec  = max(tmp, [], 1);
minVec  = min(tmp, [], 1);
spanVec = maxVec-minVec;

maxMat  = repmat(permute(maxVec, [1,3,2]), sizeVec(1),sizeVec(2), 1);
minMat  = repmat(permute(minVec, [1,3,2]), sizeVec(1),sizeVec(2), 1);
spanMat = repmat(permute(spanVec, [1,3,2]), sizeVec(1),sizeVec(2), 1);

phiNorm = (phi - minMat)./spanMat;
phiLog  = db(phiNorm);

PDMeas       = squeeze(phiNorm);
PDMeasLog    = squeeze(phiLog);
PDMeasAvg    = mean(PDMeas(:,:, 1:avgNum), 3);
PDMeasAVGLog = db(PDMeasAvg);

pdSigScan       = analysis.norm(PDMeas(:, sigUsIdx, pdMeasIdx));
pdSigScanLog    = db(pdSigScan);
pdSigScanAvg    = analysis.norm(PDMeasAvg(:, sigUsIdx));
pdSigScanLogAvg = db(pdSigScanAvg);

pdSigUs       = analysis.norm(PDMeas(sigScanIdx, :,pdMeasIdx));
pdSigUsLog    = db(pdSigUs);
pdSigUsAvg    = analysis.norm(PDMeasAvg(sigScanIdx, :));
pdSigUsLogAvg = db(pdSigUsAvg);

PDim    = analysis.norm2D(PDMeas(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh, pdMeasIdx));
PDimLog = PDMeasLog(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh, pdMeasIdx);
PDimAvg = analysis.norm2D(PDMeasAvg(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh));
PDimAvgLog = db(PDimAvg);

% Combined

minLin = min([min(PMTim(:)), min(PDim(:))]);
maxLin = max([max(PMTim(:)), max(PDim(:))]);

% minVec = mink(PMTimLog(:), 10);
% minVec = mink(PDimLog(:), 10);

minLog = min([min(PMTimLog(:)), min(PDimLog(:))]);
maxLog = max([max(PMTimLog(:)), max(PDimLog(:))]);

climLin = [minLin, maxLin];
climLog = [minLog, maxLog];

% figure(7)
% subplot(2,2,1)
% imagesc(depthAxZeroCrop, scanAxCrop, PMTim, climLin);
% axis equal tight
% title("PMT Measurement")
% xlabel("US [mm]")
% ylabel("Scan [mm]")
% colorbar
% subplot(2,2,2)
% imagesc(depthAxZeroCrop, scanAxCrop, PDim, climLin);
% axis equal tight
% title("PD Measurement")
% xlabel("US [mm]")
% ylabel("Scan [mm]")
% colorbar
% subplot(2,2,3)
% imagesc(depthAxZeroCrop, scanAxCrop, PMTimLog, climLog);
% axis equal tight
% title("PMT Measurement (db)")
% xlabel("US [mm]")
% ylabel("Scan [mm]")
% colorbar
% subplot(2,2,4)
% imagesc(depthAxZeroCrop, scanAxCrop, PDimLog, climLog);
% axis equal tight
% title("PD Measurement (db)")
% xlabel("US [mm]")
% ylabel("Scan [mm]")
% colorbar

minLin = min([min(PMTimAvg(:)), min(PDimAvg(:))]);
maxLin = max([max(PMTimAvg(:)), max(PDimAvg(:))]);

minVecPMT = unique(mink(PMTimAvgLog(:), 10));
minVecPD = unique(mink(PDimAvgLog(:), 10));

minLog = min([minVecPMT(2), minVecPD(2)]);
maxLog = max([max(PMTimAvgLog(:)), max(PDimAvgLog(:))]);

climLin = [minLin, maxLin];
climLog = [minLog, maxLog];

% figure(8)
% subplot(2,2,1)
% imagesc(depthAxZeroCrop, scanAxCrop,PMTimAvg, climLin);
% axis equal tight
% title("Average PMT Measurement ")
% xlabel("US [mm]")
% ylabel("Scan [mm]")
% colorbar
% subplot(2,2,2)
% imagesc(depthAxZeroCrop, scanAxCrop,PDimAvg, climLin);
% axis equal tight
% title("Average PD Measurement ")
% xlabel("US [mm]")
% ylabel("Scan [mm]")
% colorbar
% subplot(2,2,3)
% imagesc(depthAxZeroCrop, scanAxCrop, PMTimAvgLog, climLog);
% axis equal tight
% title("Average PMT Measurement (db)")
% xlabel("US [mm]")
% ylabel("Scan [mm]")
% colorbar
% subplot(2,2,4)
% imagesc(depthAxZeroCrop, scanAxCrop, PDimAvgLog, climLog);
% axis equal tight
% title("Average PD Measurement (db)")
% xlabel("US [mm]")
% ylabel("Scan [mm]")
% colorbar

figure(pmtMeasIdx)
subplot(2,2,1)
plot(scanAx, pmtSigScan);hold on
plot(scanAx, pdSigScan);
title("Penetration Scan")
xlabel("Scan [mm]")
ylabel("Fluence [AU]")
legend('PMT', 'PD')
subplot(2,2,2)
plot(scanAx,pmtSigScanLog); hold on
plot(scanAx,pdSigScanLog);
title("Penetration Scan (db)")
xlabel("Scan [mm]")
ylabel("Fluence [dB]")
legend('PMT', 'PD')
subplot(2,2,3)
plot(scanAx,pmtSigScanAvg); hold on
plot(scanAx,pdSigScanAvg);
title("Average Penetration Scan")
xlabel("Scan [mm]")
ylabel("Fluence [AU]")
legend('PMT', 'PD')
subplot(2,2,4)
plot(scanAx,pmtSigScanLogAvg); hold on
plot(scanAx,pdSigScanLogAvg);
title("Average Penetration Scan (db)")
xlabel("Scan [mm]")
ylabel("Fluence [dB]")
legend('PMT', 'PD')

figure(10);
subplot(2,2,1)
plot(depthAxZero,pmtSigUs); hold on
plot(depthAxZero,pdSigUs); 
title("PD Ultrasound Axis")
xlabel("US [mm]")
ylabel("Fluence [AU]")
legend('PMT', 'PD')
subplot(2,2,2)
plot(depthAxZero, pmtSigUsLog); hold on
plot(depthAxZero, pdSigUsLog);
title("PD Ultrasound Axis (db)")
xlabel("US [mm]")
ylabel("Fluence [dB]")
legend('PMT', 'PD')
subplot(2,2,3)
plot(depthAxZero, pmtSigUsAvg); hold on
plot(depthAxZero, pdSigUsAvg);
title("Averaged PD Ultrasound Axis")
xlabel("US [mm]")
ylabel("Fluence [AU]")
legend('PMT', 'PD')
subplot(2,2,4)
plot(depthAxZero, pmtSigUsLogAvg); hold on
plot(depthAxZero, pdSigUsLogAvg);
title("Averaged PD Ultrasound Axis (db)")
xlabel("US [mm]")
ylabel("Fluence [dB]")
legend('PMT', 'PD')

