close all
% clear all
clc
addpath(genpath(pwd))

%%
loadData = true;
saveFlag = false;

analysis = Analysis();

usIdxLow    = 22;
usIdxHigh   = 52;
scanIdxLow  = 1;
scanIdxHigh = 90;

numElemUS   = usIdxHigh   - usIdxLow +1;
numElemScan = scanIdxHigh - scanIdxLow +1;

sigScanIdx = 16;
sigUsIdx   = 37;

echoIdxs = [102:116, 117:187];
emIdxs = [188:251];

cutIdxs = [emIdxs, echoIdxs];
cutIdxLen = length(emIdxs)+length(echoIdxs);

avgNum = 5;


%% PMT - load and process
if loadData
    resScanPMT = load('..\Measurements\PDFinal-Exp\05-Sep-2021 12-54-44-Scan2D-PMT\2DScan-Results.mat');
    varsScanPMT = load('..\Measurements\PDFinal-Exp\05-Sep-2021 12-54-44-Scan2D-PMT\2DScan-Vars.mat');
end

pmtMeasIdx = 9;

depthAxCntr = varsScanPMT.grid.depthCntr*1e3;
depthAxZero = varsScanPMT.grid.depthZero*1e3;
scanAx      = varsScanPMT.grid.scanZero;
scanAx2     =  varsScanPMT.grid.scanVec;
depthAxCntrCrop = depthAxCntr(1:numElemScan);
depthAxZeroCrop = depthAxCntrCrop - mean(depthAxCntrCrop)*1500/1430;
scanAxCrop  = scanAx(1:numElemScan);

phi = squeeze(resScanPMT.rawPhi);
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

PMTMeas       = squeeze(phiNorm);
PMTMeasLog    = squeeze(phiLog);
PMTMeasAvg    = mean(PMTMeas(:,:, 1:avgNum), 3);
PMTMeasAVGLog = db(PMTMeasAvg);

pmtSigScan       = analysis.norm(PMTMeas(:, sigUsIdx, pmtMeasIdx));
pmtSigScanLog    = db(pmtSigScan);
pmtSigScanAvg    = analysis.norm(PMTMeasAvg(:, sigUsIdx));
pmtSigScanLogAvg = db(pmtSigScanAvg);

% figure();
% imagesc(PDMeas(:,:,1))
% figure();
% plot(PDMeas(20,:,1))

idxChopVec = [103:113, 141:155, 192:202, ];
PMTMeasClean = PMTMeas;
PMTMeasClean(:,idxChopVec,:) = [];

figure();
subplot(1,2,1)
plot(PMTMeasClean(:,:,1)')
subplot(1,2,2)
imagesc(PMTMeasClean(:,:,1))

tailIdx = 50;
pmtScanSNRMat = squeeze(max(PMTMeasClean(:,:,5),[], 2) ./ std(PMTMeasClean(:, tailIdx:end, 5), 0, 2));
pmtScanSNRMatAvg = mean(pmtScanSNRMat, 2);

figure();
plot(pmtScanSNRMatAvg)

pmtSigUs       = analysis.norm(PMTMeas(sigScanIdx, :,pmtMeasIdx));
pmtSigUsLog    = db(pmtSigUs);
pmtSigUsAvg    = analysis.norm(PMTMeasAvg(sigScanIdx, :));
pmtSigUsLogAvg = db(pmtSigUsAvg);

pmtSigUsFin                = pmtSigUs;
pmtSigUsFin(cutIdxs)       = [];
% pmtSigUsLogFin             = pmtSigUsLog;
% pmtSigUsLogFin(cutIdxs)    = [];
pmtSigUsLogFin             = circshift(pmtSigUsLog,8);
pmtSigUsLogFin(91:end)     = [];
pmtSigUsAvgFin             = pmtSigUsAvg;
pmtSigUsAvgFin(cutIdxs)    = [];
% pmtSigUsLogAvgFin          = pmtSigUsLogAvg;
% pmtSigUsLogAvgFin(cutIdxs) = [];
pmtSigUsLogAvgFin          = circshift(pmtSigUsLogAvg,8);
pmtSigUsLogAvgFin(91:end)  = [];

keppIdxsLen = length(pmtSigUs) - cutIdxLen;
depthAxZeroFin = depthAxZero (1:keppIdxsLen);
depthAxZeroFin = depthAxZero (1:90) - 50;

PMTim    = analysis.norm2D(PMTMeas(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh, pmtMeasIdx));
PMTimLog = PMTMeasLog(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh, pmtMeasIdx);
PMTimAvg = analysis.norm2D(PMTMeasAvg(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh));
PMTimAvgLog = db(PMTimAvg);

%% PD - load and process
if loadData
    resScanPD = load('..\Measurements\PDFinal-Exp\05-Sep-2021 10-48-37-Scan2D-PD\2DScan-Results.mat');
    varsScanPD = load('..\Measurements\PDFinal-Exp\05-Sep-2021 10-48-37-Scan2D-PD\2DScan-Vars.mat');
end
pdMeasIdx = 3;

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

% figure();
% imagesc(PDMeas(:,:,1))
% figure();
% plot(PDMeas(20,:,1))

idxChopVec = [103:113, 141:155, 192:202, ];
PDMeasClean = PDMeas;
PDMeasClean(:,idxChopVec,:) = [];

figure();
subplot(1,2,1)
plot(PDMeasClean(:,:,1)')
subplot(1,2,2)
imagesc(PDMeasClean(:,:,1))

tailIdx = 50;
pdScanSNRMat = squeeze(max(PDMeasClean(:,:,2),[], 2) ./ std(PDMeasClean(:, tailIdx:end, 2), 0, 2));
pdScanSNRMatAvg = mean(pdScanSNRMat, 2);

figure();
plot(pdScanSNRMatAvg)

pdSigUs       = analysis.norm(PDMeas(sigScanIdx, :,pdMeasIdx));
pdSigUsLog    = db(pdSigUs);
pdSigUsAvg    = analysis.norm(PDMeasAvg(sigScanIdx, :));
pdSigUsLogAvg = db(pdSigUsAvg);

pdSigUsFin                = pdSigUs;
pdSigUsFin(cutIdxs)          = [];
% pdSigUsLogFin             = pdSigUsLog;
% pdSigUsLogFin(cutIdxs)    = [];
pdSigUsLogFin             = circshift(pdSigUsLog, 8);
pdSigUsLogFin(91:end)     = [];
pdSigUsAvgFin             = pdSigUsAvg;
pdSigUsAvgFin(cutIdxs)    = [];
pdSigUsLogAvgFin          = pdSigUsLogAvg;
pdSigUsLogAvgFin(cutIdxs) = [];
pdSigUsLogAvgFin          = circshift(pdSigUsLogAvg, 8);
pdSigUsLogAvgFin(91:end)  = [];

PDim    = analysis.norm2D(PDMeas(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh, pdMeasIdx));
PDimLog = PDMeasLog(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh, pdMeasIdx);
PDimAvg = analysis.norm2D(PDMeasAvg(scanIdxLow:scanIdxHigh,usIdxLow:usIdxHigh));
PDimAvgLog = db(PDimAvg);

%
minLin = min([min(PMTim(:)), min(PDim(:))]);
maxLin = max([max(PMTim(:)), max(PDim(:))]);

minVecPMT = unique(mink(PMTimLog(:), 10));
minVecPD = unique(mink(PDimLog(:), 10));

minLog = min([minVecPMT(2), minVecPD(2)]);
maxLog = max([max(PMTimLog(:)), max(PDimLog(:))]);

climLin = [minLin, maxLin];
climLog = [minLog, maxLog];

%
minLinAvg = min([min(PMTimAvg(:)), min(PDimAvg(:))]);
maxLinAvg = max([max(PMTimAvg(:)), max(PDimAvg(:))]);

minVecPMTAvg = unique(mink(PMTimAvgLog(:), 10));
minVecPDAvg = unique(mink(PDimAvgLog(:), 10));

minLogAvg = min([minVecPMTAvg(2), minVecPDAvg(2)]);
maxLogAvg = max([max(PMTimAvgLog(:)), max(PDimAvgLog(:))]);

climLinAvg = [minLinAvg, maxLinAvg];
climLogAvg = [minLogAvg, maxLogAvg];

figure();
subplot(2,2,1)
imagesc(PDMeas(:,:,pdMeasIdx))
subplot(2,2,2)
imagesc(PDMeasLog(:,:,pdMeasIdx))
subplot(2,2,3)
imagesc(PMTMeas(:,:,pmtMeasIdx))
subplot(2,2,4)
imagesc(PMTMeasLog(:,:,pmtMeasIdx))

%% SNR gain vs. depth comparison
snrGainScan = pdScanSNRMatAvg./pmtScanSNRMatAvg;

figure(); 
plot(snrGainScan);

avgSNRGain = mean(snrGainScan(1:42));
snrGainScanNorm = snrGainScan/avgSNRGain;

figure();
ax1 = subplot(1,2,1);
plot(scanAx, pdScanSNRMatAvg)
hYL1 = ylabel("SNR PD");
hXL1 = xlabel("Scan[mm]");
yyaxis right
plot(scanAx, pmtScanSNRMatAvg)
hYL2 = ylabel("SNR PMT");
ax2 = subplot(1,2,2);
plot(scanAx, snrGainScanNorm)
hXL2 = xlabel("Scan[mm]");
hYL3 = ylabel("Normalized SNR Gain");

set(ax1, 'FontSize', 18, 'TickLabelInterpreter', 'latex')
set(ax2, 'FontSize', 18, 'TickLabelInterpreter', 'latex')

set(hXL1, 'FontSize', 18, 'Interpreter', 'latex')
set(hXL2, 'FontSize', 18, 'Interpreter', 'latex')

set(hYL1, 'FontSize', 18, 'Interpreter', 'latex')
set(hYL2, 'FontSize', 18, 'Interpreter', 'latex')
set(hYL3, 'FontSize', 18, 'Interpreter', 'latex')


%% fit Scan
% figure(); 
% plot(pdSigScanAvg)

axPhi = scanAx(17:70) - scanAx(17)+6;
phiPD = pdSigScanAvg(17:70);

musp = 1.5; % mm^-1;
b = 0.1806;
mua = b^2/(3*musp) * 10;


%% Figures Parameters
textFontSize  = 18;
labelFontSize = 14;
axesFontSize  = 14;
textFactor = 1.7;

%% PD Single
hPreSpace  = 5;
hPostSpace = 15;
vPreSpace  = 5;
vPostSpace = 10;

if saveFlag
    close all;
end

hFig = figure();
ax = axes(hFig);
him = imagesc(ax, 'XData', depthAxZeroCrop,...
                  'YData', scanAxCrop,...
                  'CData', PDimLog);
axis equal 
axis tight
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hXL = xlabel(ax, "US (Z) [mm]");
set(hXL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hYL = ylabel(ax, "Scan (Y) [mm]");
set(hYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hTxt = text(depthAxZeroCrop(1), scanAxCrop(end-8)  ,'(b)');
set(hTxt, 'FontSize', textFontSize, 'Interpreter', 'latex', 'color', 'w')

hCB = colorbar();
caxis(ax, climLog)
set(hCB, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hCYL = ylabel(hCB, "Fluence(dB)");
set(hCYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
set(hCB, 'Units', 'pixels');
cbPos = get(hCB, 'Position');
cbWidth = cbPos(3) + (axesFontSize + labelFontSize)*1.7; 
cbHeight = cbPos(4);

set(ax, 'Units', 'pixels');
axPos = get(ax, 'Position');
axWidth = axPos(3); 
axHeight = axPos(4);
axLeft   = hPreSpace + (labelFontSize + axesFontSize)*1.7;
axBottom = -axHeight/4 + vPreSpace + (axesFontSize + labelFontSize)*1.7;
set(ax, 'Position', [axLeft, axBottom, axWidth, axHeight])

axBottom = 0;
figWidth  = hPreSpace + (labelFontSize + axesFontSize)*1.7 + axWidth  + cbWidth  + hPostSpace;
figHeight = vPreSpace + (labelFontSize + axesFontSize)*1.7  + axHeight/2     + vPostSpace;
figPos = get(hFig, 'Position');
set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figHeight])

%%
if saveFlag
    filename = '..\..\Publications\Homodyne AO\Figures\FluenceMap-PD-Single.pdf';
    set(hFig, 'Units', 'centimeters');
    figPosCM = get(hFig, 'Position');
    hFig.PaperUnits = 'centimeters';
    hFig.PaperType = '<custom>';
    hFig.PaperSize = [figPosCM(3), figPosCM(4)];
    hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
    print(hFig, filename, '-dpdf', '-r0') 
end
%% PMT Single
hPreSpace  = 5;
hPostSpace = 15;
vPreSpace  = 5;
vPostSpace = 10;

if saveFlag
    close all;
end

hFig = figure();
ax = axes(hFig);
him = imagesc(ax, 'XData', depthAxZeroCrop,...
                  'YData', scanAxCrop,...
                  'CData', PMTimLog);
axis equal 
axis tight
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hYL = ylabel("Scan (Y) [mm]");
set(hYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hTxt = text(depthAxZeroCrop(1), scanAxCrop(end-8)  ,'(a)');
set(hTxt, 'FontSize', textFontSize, 'Interpreter', 'latex', 'color', 'w')
set(ax,'xtick',[])

hCB = colorbar();
caxis(ax, climLog)
set(hCB, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hCYL = ylabel(hCB, "Fluence(dB)");
set(hCYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')

set(ax, 'Units', 'pixels');
axPos = get(ax, 'Position');
axWidth = axPos(3); 
axHeight = axPos(4);
axLeft   = hPreSpace + (labelFontSize + axesFontSize)*textFactor;
axBottom = -axHeight/4+vPreSpace;
set(ax, 'Position', [axLeft, axBottom, axWidth, axHeight])

figHeight = figHeight - (labelFontSize + axesFontSize)*textFactor;
set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figHeight])
%%
if saveFlag
    filename = '..\..\Publications\Homodyne AO\Figures\FluenceMap-PMT-Single.pdf';
    set(hFig, 'Units', 'centimeters');
    figPosCM = get(hFig, 'Position');
    hFig.PaperUnits = 'centimeters';
    hFig.PaperType = '<custom>';
    hFig.PaperSize = [figPosCM(3), figPosCM(4)];
    hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
    print(hFig, filename, '-dpdf', '-r0') 
end
%% PD Avg
hPreSpace  = 5;
hPostSpace = 15;
vPreSpace  = 5;
vPostSpace = 5;

if saveFlag
    close all;
end

hFig = figure();
ax = axes(hFig);
him = imagesc(ax, 'XData', depthAxZeroCrop,...
                  'YData', scanAxCrop,...
                  'CData', PDimAvgLog);
axis equal 
axis tight
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hXL = xlabel(ax, "US (Z) [mm]");
set(hXL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hYL = ylabel(ax, "Scan (Y) [mm]");
set(hYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hTxt = text(depthAxZeroCrop(1), scanAxCrop(end-8)  ,'(b)');
set(hTxt, 'FontSize', textFontSize, 'Interpreter', 'latex', 'color', 'w')

hCB = colorbar();
caxis(ax, climLogAvg)
set(hCB, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hCYL = ylabel(hCB, "Fluence(dB)");
set(hCYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
set(hCB, 'Units', 'pixels');
cbPos = get(hCB, 'Position');
cbWidth = cbPos(3) + (axesFontSize + labelFontSize)*textFactor; 
cbHeight = cbPos(4);

set(ax, 'Units', 'pixels');
axPos = get(ax, 'Position');
axWidth = axPos(3); 
axHeight = axPos(4);
axLeft   = hPreSpace + (labelFontSize + axesFontSize)*textFactor;
axBottom = -axHeight/4 + vPreSpace + (axesFontSize + labelFontSize)*textFactor;
set(ax, 'Position', [axLeft, axBottom, axWidth, axHeight])

axBottom = 0;
figWidth  = hPreSpace + (labelFontSize + axesFontSize)*textFactor + axWidth  + cbWidth  + hPostSpace;
figHeight = vPreSpace + (labelFontSize + axesFontSize)*textFactor  + axHeight/2     + vPostSpace;
figPos = get(hFig, 'Position');
set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figHeight])
%%
if saveFlag
    filename = '..\..\Publications\Homodyne AO\Figures\FluenceMap-PD-Avg.pdf';
    set(hFig, 'Units', 'centimeters');
    figPosCM = get(hFig, 'Position');
    hFig.PaperUnits = 'centimeters';
    hFig.PaperType = '<custom>';
    hFig.PaperSize = [figPosCM(3), figPosCM(4)];
    hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
    print(hFig, filename, '-dpdf', '-r0') 
end
%% PMT Avg
hPreSpace  = 5;
hPostSpace = 15;
vPreSpace  = 5;
vPostSpace = 10;

textFactor = 1.7;

if saveFlag
    close all;
end


hFig = figure();
ax = axes(hFig);
him = imagesc(ax, 'XData', depthAxZeroCrop,...
                  'YData', scanAxCrop,...
                  'CData', PMTimAvgLog);
axis equal 
axis tight
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hYL = ylabel("Scan (Y) [mm]");
set(hYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hTxt = text(depthAxZeroCrop(1), scanAxCrop(end-8)  ,'(a)');
set(hTxt, 'FontSize', textFontSize, 'Interpreter', 'latex', 'color', 'w')
set(ax,'xtick',[])

hCB = colorbar();
caxis(ax, climLogAvg)
set(hCB, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hCYL = ylabel(hCB, "Fluence(dB)");
set(hCYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
set(hCB, 'Units', 'pixels');
cbPos = get(hCB, 'Position');
cbWidth = cbPos(3) + (axesFontSize + labelFontSize)*textFactor; 
cbHeight = cbPos(4);

set(ax, 'Units', 'pixels');
axPos = get(ax, 'Position');
axWidth = axPos(3); 
axHeight = axPos(4);
axLeft   = hPreSpace + (labelFontSize + axesFontSize)*textFactor;
axBottom = -axHeight/4+vPreSpace;
set(ax, 'Position', [axLeft, axBottom, axWidth, axHeight])

axBottom = 0;

figHeight = figHeight - (labelFontSize + axesFontSize)*textFactor;
set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figHeight])
%%
if saveFlag
    filename = '..\..\Publications\Homodyne AO\Figures\FluenceMap-PMT-Avg.pdf';
    set(hFig, 'Units', 'centimeters');
    figPosCM = get(hFig, 'Position');
    hFig.PaperUnits = 'centimeters';
    hFig.PaperType = '<custom>';
    hFig.PaperSize = [figPosCM(3), figPosCM(4)];
    hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
    print(hFig, filename, '-dpdf', '-r0') 
end
%% FFT Comparison
textFontSize  = 28;
labelFontSize = 24;
axesFontSize  = 18;
textFactor    = 1.7;


% hFig = figure();
% ax = axes(hFig);
% plot(ax, scanAx2, pmtSigScanLog); hold on
% plot(ax, scanAx2, pdSigScanLog);

sigUsIdx = 37;
bkgIdx = 27;
pmtMeasIdx = 9;
pdMeasIdx  = 3;
if loadData
    aoResPMTPeak = load('..\Measurements\PDFinal-Exp\AO-R-9-Y-7.50-Results-PMT.mat');
    str = sprintf("../Measurements/PDFinal-Exp/05-Sep-2021 12-54-44-Scan2D-PMT/AOResults/AO-R-%d-Y-%.2f-Results.mat", pmtMeasIdx, bkgIdx);
    aoResPMTBkg = load(str);
    aoResPDPeak = load('..\Measurements\PDFinal-Exp\AO-R-3-Y-7.50-Results-PD.mat');
    str = sprintf('../Measurements/PDFinal-Exp/05-Sep-2021 10-48-37-Scan2D-PD/AOResults/AO-R-%d-Y-%.2f-Results.mat', pdMeasIdx, bkgIdx);
    aoResPDBkg = load(str);
end
algo = Algo();
% varsScanPD.ao.extVars.algo.useGPU = false;
% varsScanPMT.ao.extVars.algo.useGPU = false;
% uVarsAlgo = algo.createUserVars();
% uVarsAlgo.export

uVarsAlgoPD                     = varsScanPD.ao.extVars.algo;
uVarsAlgoPD.fSqnc               = uVarsAlgoPD.fTrain;
uVarsAlgoPD.cycPerPulse         = uVarsAlgoPD.cycInPulse;
uVarsAlgoPD.fs                  = uVarsAlgoPD.fExtClk;
uVarsAlgoPD.fgClk               = uVarsAlgoPD.fSclk;
uVarsAlgoPD.sClkDcyc            = uVarsAlgoPD.extClkDcyc;
uVarsAlgoPD.frameTime           = uVarsAlgoPD.quantTime;
uVarsAlgoPD.useFrame            = uVarsAlgoPD.useQuant;
uVarsAlgoPD.analyzeSingleCh     = false;
uVarsAlgoPD.contSpeckleAnalysis = false;
uVarsAlgoPD.cutArtfct           = false;
uVarsAlgoPD.contHadamard        = false;
uVarsAlgoPD.artfctIdxVec        = [];
uVarsAlgoPD.export.meas = false;
uVarsAlgoPD.export.signal = false;
uVarsAlgoPD.export.deMul = false;
uVarsAlgoPD.export.reshaped = false;
uVarsAlgoPD.export.fft = false;
uVarsAlgoPD.export.usCompCmplx = false;

uVarsAlgoPMT                     = varsScanPMT.ao.extVars.algo;
uVarsAlgoPMT.fSqnc               = uVarsAlgoPMT.fTrain;
uVarsAlgoPMT.cycPerPulse         = uVarsAlgoPMT.cycInPulse;
uVarsAlgoPMT.fs                  = uVarsAlgoPMT.fExtClk;
uVarsAlgoPMT.fgClk               = uVarsAlgoPMT.fSclk;
uVarsAlgoPMT.sClkDcyc            = uVarsAlgoPMT.extClkDcyc;
uVarsAlgoPMT.frameTime           = uVarsAlgoPMT.quantTime;
uVarsAlgoPMT.useFrame            = uVarsAlgoPMT.useQuant;
uVarsAlgoPMT.analyzeSingleCh     = false;
uVarsAlgoPMT.contSpeckleAnalysis = false;
uVarsAlgoPMT.cutArtfct           = false;
uVarsAlgoPMT.contHadamard        = false;
uVarsAlgoPMT.artfctIdxVec        = [];
uVarsAlgoPMT.export.meas = false;
uVarsAlgoPMT.export.signal = false;
uVarsAlgoPMT.export.deMul = false;
uVarsAlgoPMT.export.reshaped = false;
uVarsAlgoPMT.export.fft = false;
uVarsAlgoPMT.export.usCompCmplx = false;
% figure();

%---------------------
%PD Reprocessing
%---------------------
algo.setVars(uVarsAlgoPD);
varsAlgoPD = algo.getVars();

% Peak
algo.curRes.qAvgChFFT = fftshift(aoResPDPeak.qAvgChFFT,2);
algo.extractPhiFromFFT();
resPDPeak     = algo.curRes;
% subplot(2,2,3); stem(resPDPeak.phi); title("PD Peak");
resPDPeakNorm = resPDPeak.fittedChAvgFFT(:, sigUsIdx)-1;
maxPDPeak = max(resPDPeakNorm);
minPDPeak = min(resPDPeakNorm);
spanPDPeak = maxPDPeak-minPDPeak;
resPDPeakNorm = analysis.norm(resPDPeakNorm);


%Bag
algo.curRes.qAvgChFFT = fftshift(aoResPDBkg.qAvgChFFT,2);
algo.extractPhiFromFFT();
resPDBkg     = algo.curRes;
% subplot(2,2,4); stem(resPDBkg.phi); title("PD Bkg");
resPDBkgNorm = (resPDBkg.fittedChAvgFFT(:, sigUsIdx))-1;
minPDBkg     = min(resPDBkgNorm);
spanPDBkg    = maxPDPeak - minPDPeak;
resPDBkgNorm = (resPDBkgNorm - minPDBkg) ./ spanPDPeak;


%---------------------
%PMT Reprocessing
%---------------------
algo.setVars(uVarsAlgoPMT);
varsAlgoPMT = algo.getVars();

%Peak
algo.curRes.qAvgChFFT = fftshift(aoResPMTPeak.qAvgChFFT,2);
algo.extractPhiFromFFT();
resPMTPeak     = algo.curRes;
% subplot(2,2,1); stem(resPMTPeak.phi); title("PMT Peak");
resPMTPeakNorm = resPMTPeak.fittedChAvgFFT(:, sigUsIdx)-1;
maxPMTPeak   = max(resPMTPeakNorm);
minPMTPeak   = min(resPMTPeakNorm);
spanPMTPeak  = maxPMTPeak-minPMTPeak;
resPMTPeakNorm = (resPMTPeakNorm - minPMTPeak) ./ spanPMTPeak;

% resPMTFFTPeak  = resPMTPeak.fittedFFT(:, sigUsIdx);
% resPMTPeakNorm = analysis.norm(resPMTFFTPeak);

%Bkg
algo.curRes.qAvgChFFT = fftshift(aoResPMTBkg.qAvgChFFT,2);
algo.extractPhiFromFFT();
resPMTBkg = algo.curRes;
% subplot(2,2,2); stem(resPMTBkg.phi); title("PMT Bkg");
resPMTBkgNorm = (resPMTBkg.fittedChAvgFFT(:, sigUsIdx))-1;
minPMTBkg     = min(resPMTBkgNorm);
spanPMTBkg    = maxPMTPeak - minPMTPeak;
resPMTBkgNorm = (resPMTBkgNorm - minPMTPeak) ./ spanPMTBkg;

% figure();
% plot(fBarPMT*1e-6, resPMTBkgNorm); hold on
% plot(fBarPD*1e-6,  resPDBkgNorm);

% resPMTFFTBkg = resPMTBkg.fittedFFT(:, sigUsIdx);
% resPMTBkgNorm = analysis.norm(resPMTFFTBkg);

%---------------------------------
%Frquency Domain Post Processing
%---------------------------------
fBarPMT = varsAlgoPMT.freq.fBar;
fBarPD  = varsAlgoPD.freq.fBar;

resPMTPeakLog = db(resPMTPeakNorm, 'power');
resPMTBkgLog  = db(resPMTBkgNorm,  'power');
resPDPeakLog  = db(resPDPeakNorm,  'power');
resPDBkgLog   = db(resPDBkgNorm,   'power');

threshold = -40;

resPMTBkgLog (resPMTBkgLog<threshold) = -inf;
idxsInf = find(resPMTBkgLog==-inf);
for i=1:length(idxsInf)
    resPMTBkgLog(idxsInf(i)) = 0.5*(resPMTBkgLog(idxsInf(i)-1)+resPMTBkgLog(idxsInf(i)+1));
end

resPDBkgLog (resPDBkgLog<threshold) = -inf;
idxsInf = find(resPDBkgLog==-inf);
for i=1:length(idxsInf)
    resPDBkgLog(idxsInf(i)) = 0.5*(resPDBkgLog(idxsInf(i)-1)+resPDBkgLog(idxsInf(i)+1));
end

resPMTPeakLog (resPMTPeakLog<threshold) = -inf;
idxsInf = find(resPMTPeakLog==-inf);
for i=1:length(idxsInf)
    resPMTPeakLog(idxsInf(i)) = 0.5*(resPMTPeakLog(idxsInf(i)-1)+resPMTPeakLog(idxsInf(i)+1));
end

resPDPeakLog (resPDPeakLog<threshold) = -inf;
idxsInf = find(resPDPeakLog==-inf);
for i=1:length(idxsInf)
    resPDPeakLog(idxsInf(i)) = 0.5*(resPDPeakLog(idxsInf(i)-1)+resPDPeakLog(idxsInf(i)+1));
end

figure();
subplot(2,2,1)
plot(fBarPMT/1e6,resPMTPeakLog); hold on
plot(fBarPD/1e6,resPDPeakLog);
title("Peak db")
ylim([threshold, 0])
subplot(2,2,2)
plot(fBarPMT/1e6,resPMTBkgLog); hold on
plot(fBarPD/1e6,resPDBkgLog);
title("Bkg db")
ylim([threshold, 0])
subplot(2,2,3)
plot(fBarPMT/1e6,resPMTPeakNorm); hold on
plot(fBarPD/1e6,resPDPeakNorm);
title("Peak Lin")
subplot(2,2,4)
plot(fBarPMT/1e6,resPMTBkgNorm); hold on
plot(fBarPD/1e6,resPDBkgNorm);
title("Bkg Lin")

minPeak = mink([resPMTPeakLog', resPDPeakLog'],1);
minPeak(minPeak == -inf) = [];
maxPeak = 0;

minBkg = mink([resPMTBkgLog', resPDBkgLog'],1);
minBkg(minBkg <-100) = [];
maxBkg = 0;

minY = min ([minPeak, minBkg]);
maxY = max([maxPeak, maxBkg]);

%% FFT Comparison Peak
hPreSpace  = 5;
hPostSpace = 40;
vPreSpace  = 5;
vPostSpace = 5;

hFig = figure();
ax = axes(hFig);
plot(fBarPMT/1e6, resPMTPeakNorm); hold on
plot(fBarPD/1e6, resPDPeakNorm);
xlim([-5, 5])
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hXL = xlabel(ax, "f [MHz]");
set(ax, 'XTick', [-5, -3.75, -2.5, -1.25, 0, 1.25, 2.5 3.75, 5])
set(hXL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hYL = ylabel(ax, "Spectrum [AU]");
set(hYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hTxt = text(-5, 0.9  ,'(a)');
set(hTxt, 'FontSize', textFontSize, 'Interpreter', 'latex', 'color', 'k')

legPos = [0.809878259495187,0.803571428571429,0.152884613751219,0.106746031746032];

hLeg = legend(ax, "PMT", "PD");
set(hLeg, 'Position', legPos, 'FontSize', 10)

set(ax, 'Units', 'pixels');
axPos = get(ax, 'Position');

figPos = get(hFig, 'Position');
figWidth = figPos(3) - hPostSpace;
set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figPos(4)])
figPos = get(hFig, 'Position');

%%
if saveFlag
    filename = '..\..\Publications\Homodyne AO\Figures\FFT-Comparison-Peak.pdf';
    set(hFig, 'Units', 'centimeters');
    figPosCM = get(hFig, 'Position');
    hFig.PaperUnits = 'centimeters';
    hFig.PaperType = '<custom>';
    hFig.PaperSize = [figPosCM(3), figPosCM(4)];
    hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
    print(hFig, filename, '-dpdf', '-r0')
end

%% FFT Comparison Bkg
hPreSpace  = 10;
hPostSpace = 0;    
vPreSpace  = 5;
vPostSpace = 5;

if saveFlag
    close all;
end

hFig = figure();
ax = axes(hFig);
plot(fBarPMT/1e6, resPMTBkgNorm); hold on
plot(fBarPD/1e6, resPDBkgNorm);

xlim([-5, 5])
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hXL = xlabel(ax, "f [MHz]");
set(hXL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hTxt = text(-5, 0.9  ,'(b)');
set(hTxt, 'FontSize', textFontSize, 'Interpreter', 'latex', 'color', 'k')
set(ax, 'ytick', []);
set(ax, 'XTick', [-5, -3.75, -2.5, -1.25, 0, 1.25, 2.5 3.75, 5])

legPos = [0.778794770375394,0.813095238095238,0.177257523189819,0.100793650793651];

hLeg = legend("PMT", "PD");
set(hLeg, 'Position', legPos, 'FontSize', 10)

set(ax, 'Units', 'pixels');
set(ax, 'Position', [hPreSpace, axPos(2), axPos(3), axPos(4)])

figWidth  = figPos(3) - (labelFontSize + axesFontSize)*textFactor-hPostSpace;
set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figPos(4)])
%%
if saveFlag
    filename = '..\..\Publications\Homodyne AO\Figures\FFT-Comparison-bkg.pdf';
    set(hFig, 'Units', 'centimeters');
    figPosCM = get(hFig, 'Position');
    hFig.PaperUnits = 'centimeters';
    hFig.PaperType = '<custom>';
    hFig.PaperSize = [figPosCM(3), figPosCM(4)];
    hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
    print(hFig, filename, '-dpdf', '-r0')
end


%% FFT Comparison Peak - Log
hPreSpace  = 5;
hPostSpace = 40;
vPreSpace  = 5;
vPostSpace = 5;

hFig = figure();
ax = axes(hFig);
plot(fBarPMT/1e6,resPMTPeakLog); hold on
plot(fBarPD/1e6,resPDPeakLog);
xlim([-5, 5])
% ylim(ax, [minY, maxY])
set(ax, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
hXL = xlabel(ax, "f [MHz]");
set(ax, 'XTick', [-5, -3.75, -2.5, -1.25, 0, 1.25, 2.5 3.75, 5])
set(hXL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hYL = ylabel(ax, "Power Spectrum [dB]");
set(hYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
% hTxt = text(-5, -45.094  ,'(a)');
% set(hTxt, 'FontSize', textFontSize, 'Interpreter', 'latex', 'color', 'k')

hLeg = legend("PMT", "PD");
set(hLeg, 'Location', 'northeast', 'FontSize', 12)

set(ax, 'Units', 'pixels');
axPos = get(ax, 'Position');

figPos = get(hFig, 'Position');
figWidth = figPos(3) - hPostSpace;
set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figPos(4)])
figPos = get(hFig, 'Position');

%%
if saveFlag
    filename = '..\..\Publications\Homodyne AO\Figures\FFT-Comparison-Peak-Log.pdf';
    set(hFig, 'Units', 'centimeters');
    figPosCM = get(hFig, 'Position');
    hFig.PaperUnits = 'centimeters';
    hFig.PaperType = '<custom>';
    hFig.PaperSize = [figPosCM(3), figPosCM(4)];
    hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
    print(hFig, filename, '-dpdf', '-r0')
end

%% FFT Comparison Bkg - Log
hPreSpace  = 10;
hPostSpace = 40;    
vPreSpace  = 5;
vPostSpace = 5;

% if saveFlag
%     close all;
% end

hFig = figure();
ax = axes(hFig);
plot(fBarPMT/1e6,resPMTBkgLog); hold on
plot(fBarPD/1e6,resPDBkgLog);

xlim([-5, 5])
% ylim(ax, [minY, maxY])
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hXL = xlabel(ax, "f [MHz]");
set(hXL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hYL = ylabel(ax, "Power Spectrum [dB]");
set(hYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hTxt = text(-5, -27.1326  ,'(b)');
set(hTxt, 'FontSize', textFontSize, 'Interpreter', 'latex', 'color', 'k')
% set(ax, 'ytick', []);
set(ax, 'XTick', [-5, -3.75, -2.5, -1.25, 0, 1.25, 2.5 3.75, 5])


hLeg = legend("PMT", "PD");
set(hLeg, 'Location', 'southeast', 'FontSize', 12)

set(ax, 'Units', 'pixels');
% set(ax, 'Position', [hPreSpace, axPos(2), axPos(3), axPos(4)])

% figWidth  = figPos(3) - (labelFontSize + axesFontSize)*textFactor-hPostSpace;
% set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figPos(4)])

figPos = get(hFig, 'Position');
figWidth = figPos(3) - hPostSpace;
set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figPos(4)])
figPos = get(hFig, 'Position');

%%
if saveFlag
    filename = '..\..\Publications\Homodyne AO\Figures\FFT-Comparison-bkg-Log.pdf';
    set(hFig, 'Units', 'centimeters');
    figPosCM = get(hFig, 'Position');
    hFig.PaperUnits = 'centimeters';
    hFig.PaperType = '<custom>';
    hFig.PaperSize = [figPosCM(3), figPosCM(4)];
    hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
    print(hFig, filename, '-dpdf', '-r0')
end


%%
textFontSize  = 28;
labelFontSize = 24;
axesFontSize  = 24;
textFactor = 1.7;

tmp = mink([pmtSigScanLog', pdSigScanLog'],3);
tmp(tmp == -inf) = [];

minY = min ([min([pmtSigUsLogFin, pdSigUsLogFin]), tmp]);
maxY = max ([max([pmtSigUsLogFin, pdSigUsLogFin]), max([pmtSigScanLog,pdSigScanLog])]);

%% US Comparison
hPreSpace  = 5;
hPostSpace = 40;
vPreSpace  = 5;
vPostSpace = 5;

if saveFlag
    close all;
end

hFig = figure();
ax = axes(hFig);
plot(ax, depthAxZeroFin, pmtSigUsLogFin); hold on
plot(ax, depthAxZeroFin, pdSigUsLogFin);
ylim(ax, [minY, maxY])
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hXL = xlabel(ax, "US (Z) [mm]");
set(hXL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hYL = ylabel(ax, "Fluence Rate [dB]");
set(hYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hTxt = text(depthAxZeroFin(1), -60  ,'(a)');
set(hTxt, 'FontSize', textFontSize, 'Interpreter', 'latex', 'color', 'k')

legend("PMT", "PD")

set(ax, 'Units', 'pixels');
axPos = get(ax, 'Position');

figPos = get(hFig, 'Position');
figWidth = figPos(3) - hPostSpace;
set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figPos(4)])
figPos = get(hFig, 'Position');
%%
if saveFlag
    filename = '..\..\Publications\Homodyne AO\Figures\US-Comparison.pdf';
    set(hFig, 'Units', 'centimeters');
    figPosCM = get(hFig, 'Position');
    hFig.PaperUnits = 'centimeters';
    hFig.PaperType = '<custom>';
    hFig.PaperSize = [figPosCM(3), figPosCM(4)];
    hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
    print(hFig, filename, '-dpdf', '-r0')
end
%% Scan Comparison
hPreSpace  = 10;
hPostSpace = 0;    
vPreSpace  = 5;
vPostSpace = 5;

if saveFlag
    close all;
end

% figure()
% plot(pmtSigScanLog); hold on
% plot(pdSigScanLog)

hFig = figure();
ax = axes(hFig);
plot(ax, scanAx(15:end)-scanAx(15), pmtSigScanLog(15:end)); hold on
plot(ax, scanAx(15:end)-scanAx(15), pdSigScanLog(15:end));
ylim(ax, [minY, maxY])
xlim(ax, [0, 53])
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hXL = xlabel(ax, "Scan (Y) [mm]");
set(hXL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hTxt = text(scanAx(1), -60  ,'(b)');
set(hTxt, 'FontSize', textFontSize, 'Interpreter', 'latex', 'color', 'k')
set(ax, 'ytick', []);

legend("PMT", "PD")

set(ax, 'Units', 'pixels');
set(ax, 'Position', [hPreSpace, axPos(2), axPos(3), axPos(4)])

figWidth  = figPos(3) - (labelFontSize + axesFontSize)*textFactor-hPostSpace;
set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figPos(4)])
%%
if saveFlag
    filename = '..\..\Publications\Homodyne AO\Figures\Scan-Comparison.pdf';
    set(hFig, 'Units', 'centimeters');
    figPosCM = get(hFig, 'Position');
    hFig.PaperUnits = 'centimeters';
    hFig.PaperType = '<custom>';
    hFig.PaperSize = [figPosCM(3), figPosCM(4)];
    hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
    print(hFig, filename, '-dpdf', '-r0')
end

%%
tmp = mink([pmtSigScanLogAvg', pdSigScanLogAvg'],3);
tmp(tmp == -inf) = [];

minY = min ([min([pmtSigUsLogAvgFin, pdSigUsLogAvgFin]), tmp]);
maxY = max ([max([pmtSigUsLogAvgFin, pdSigUsLogAvgFin]), max([pmtSigScanLogAvg', pdSigScanLogAvg'])]);
%% US Comparison - Avg
hPreSpace  = 5;
hPostSpace = 40;
vPreSpace  = 5;
vPostSpace = 5;

if saveFlag
    close all;
end

hFig = figure();
ax = axes(hFig);
plot(ax, depthAxZeroFin, pmtSigUsLogAvgFin); hold on
plot(ax, depthAxZeroFin, pdSigUsLogAvgFin);

set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hXL = xlabel(ax, "US (Z) [mm]");
set(hXL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hYL = ylabel(ax, "Fluence [AU]");
set(hYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hTxt = text(depthAxZeroFin(1), -10  ,'(a)');
set(hTxt, 'FontSize', textFontSize, 'Interpreter', 'latex', 'color', 'k')
ylim(ax, [minY, maxY])
legend("PMT", "PD")

set(ax, 'Units', 'pixels');
axPos = get(ax, 'Position');

figPos = get(hFig, 'Position');
figWidth = figPos(3) - hPostSpace;
set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figPos(4)])
figPos = get(hFig, 'Position');
%%
if saveFlag
    filename = '..\..\Publications\Homodyne AO\Figures\US-Comparison-Avg.pdf';
    set(hFig, 'Units', 'centimeters');
    figPosCM = get(hFig, 'Position');
    hFig.PaperUnits = 'centimeters';
    hFig.PaperType = '<custom>';
    hFig.PaperSize = [figPosCM(3), figPosCM(4)];
    hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
    print(hFig, filename, '-dpdf', '-r0')
end

%% Scan Comparison - AVG
hPreSpace  = 10;
hPostSpace = 0;    
vPreSpace  = 5;
vPostSpace = 5;

if saveFlag
    close all;
end

% figure();
% plot(pmtSigScanLogAvg); hold on
% plot(pdSigScanLogAvg);

hFig = figure();
ax = axes(hFig);
plot(ax, scanAx(18:end) - scanAx(18), pmtSigScanLogAvg(18:end)); hold on
plot(ax, scanAx(18:end) - scanAx(18), pdSigScanLogAvg(18:end));

set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hXL = xlabel(ax, "Scan (Y) [mm]");
set(hXL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hTxt = text(0, -55  ,'(b)');
set(hTxt, 'FontSize', textFontSize, 'Interpreter', 'latex', 'color', 'k')
set(ax, 'ytick', []);
ylim(ax, [minY, maxY])
legend("PMT", "PD")

set(ax, 'Units', 'pixels');
set(ax, 'Position', [hPreSpace, axPos(2), axPos(3), axPos(4)])

figWidth  = figPos(3) - (labelFontSize + axesFontSize)*textFactor-hPostSpace;
set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figPos(4)])
%%
if saveFlag
    filename = '..\..\Publications\Homodyne AO\Figures\Scan-Comparison-Avg.pdf';
    set(hFig, 'Units', 'centimeters');
    figPosCM = get(hFig, 'Position');
    hFig.PaperUnits = 'centimeters';
    hFig.PaperType = '<custom>';
    hFig.PaperSize = [figPosCM(3), figPosCM(4)];
    hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
    print(hFig, filename, '-dpdf', '-r0')
end
%% SNR Vs. Ref
if loadData
    snrVsRefRes =load("../Measurements/PDFinal-Exp/PD-2sT-5KhzRR-5HzWL-20MHzFs-1msQ-20reps-SNR-vs-Ref-Fine/05-Sep-2021 17-07-49-PD-2sT-5KhzRR-5HzWL-20MHzFs-1msQ-20reps-SNR-vs-Ref-Fine-ExpVars.mat");
end
data = snrVsRefRes.expData.data;
vars = snrVsRefRes.expData.vars;

SNRPMT   = 150;
avgSNRPD = data.SNRAvgMeas;
stdSNRPD = data.SNRstdMeas;

avgSNRGain = avgSNRPD./SNRPMT;
stdSNRGain = stdSNRPD./SNRPMT;

avgSNRGain(end) = [];
stdSNRGain(end) = [];

refAx = vars.pwr.refFactor;
intAx = vars.pwr.Ivec;
refAx(end) = [];

avgSNRGainLog = log(avgSNRGain);
refAxLog      = log( refAx);

fitModel = fit([refAxLog(1:7)]', [avgSNRGainLog(1:7)]', 'poly1');
fitCurve = fitModel.p1*refAxLog + fitModel.p2;

figure();
plot(refAxLog, avgSNRGainLog); hold on;
plot(refAxLog, fitCurve)


%% SNR Vs. Ref - linear
hPreSpace  = 5;
hPostSpace = 40;
vPreSpace  = 5;
vPostSpace = 5;

if saveFlag
    close all;
end

hFig = figure();
ax = axes(hFig);
errorbar(ax, refAx, avgSNRGain, stdSNRGain)

set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hXL = xlabel(ax, "$\gamma$");
set(hXL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hYL = ylabel(ax, "SNR Gain [AU]");
set(hYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
% hTxt = text(depthAxZeroFin(1), -10  ,'(a)');
% set(hTxt, 'FontSize', textFontSize, 'Interpreter', 'latex', 'color', 'k')
xlim(ax, [20,32000])
set(ax, 'xScale', 'log')
% set(ax, 'yScale', 'log')
ylim(ax, [-0.1, 4.9])
xlim (ax, [3e1, 3.5e4])
% set(ax, 'Units', 'pixels');
% axPos = get(ax, 'Position');

% figPos = get(hFig, 'Position');
% figWidth = figPos(3) - hPostSpace;
% set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figPos(4)])
% figPos = get(hFig, 'Position');
%%
if saveFlag
    filename = '..\..\Publications\Homodyne AO\Figures\SNR-vs.Ref.pdf';
    set(hFig, 'Units', 'centimeters');
    figPosCM = get(hFig, 'Position');
    hFig.PaperUnits = 'centimeters';
    hFig.PaperType = '<custom>';
    hFig.PaperSize = [figPosCM(3), figPosCM(4)];
    hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
    print(hFig, filename, '-dpdf', '-r0')
end

%% Figures Parameters
textFontSize  = 18;
labelFontSize = 20;
axesFontSize  = 20;
textFactor = 1.7;

%% SNR Vs. Ref - ylog
hFig = figure();
ax = axes(hFig);
errorbar(ax, refAx, avgSNRGain, stdSNRGain); hold on
plot(exp(refAxLog), exp(fitCurve))

set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
hXL = xlabel(ax, "$\gamma$");
set(hXL, 'FontSize', labelFontSize, 'Interpreter', 'latex')
hYL = ylabel(ax, "Optical Gain");
set(hYL, 'FontSize', labelFontSize, 'Interpreter', 'latex')

xlim(ax, [20,32000])
set(ax, 'xScale', 'log')
set(ax, 'yScale', 'log')
ylim(ax, [-0.1, 4.9])
xlim (ax, [3e1, 3.5e4])
set(ax, 'YTick', [0.03,0.1,0.3,1,3])
str = sprintf("y=%.2fx%.2f", fitModel.p1, fitModel.p2);
hTxt = text(10^3, 0.2  ,str);
set(hTxt, 'FontSize', 18, 'Interpreter', 'latex', 'color', 'k')

hLeg = legend("Measured", "Fit");
set(hLeg, 'Location', 'northwest', 'FontSize', 14)

if saveFlag
    filename = '..\..\Publications\Homodyne AO\Figures\SNR-vs.Ref-ylog.pdf';
    set(hFig, 'Units', 'centimeters');
    figPosCM = get(hFig, 'Position');
    hFig.PaperUnits = 'centimeters';
    hFig.PaperType = '<custom>';
    hFig.PaperSize = [figPosCM(3), figPosCM(4)];
    hFig.PaperPosition = [0, 0, figPosCM(3)+1.5, figPosCM(4)+0.75];
    print(hFig, filename, '-dpdf', '-r0')
end

