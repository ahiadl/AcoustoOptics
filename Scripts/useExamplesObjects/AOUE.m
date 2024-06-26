% close all;
% clear all;
% clc;
% 
% ao = acoustoOptics();
% ao.init();
% % stages = stages('COM3');
% stages.connect();

% stages.assignStagesAxes(['X','Y'], [1,2]);
% stages.moveStageAxisAbs('X', 101);
% stages.moveStageAxisAbs('Y', 88);
%%
uVars = ao.createUserVars();

%-------------------------
% Algo and Peripheral Vars
%-------------------------
% US Signal
uVars.ao.cycPerPulse = 1;      %[#] 
uVars.ao.fSin        = 1.25e6; %[Hz]              
uVars.ao.fSqnc       = 5e3;    %[Hz]
uVars.ao.frameTime   = 0.002;  %[s]

% Sampling Clk
uVars.ao.fs       = 20e6;   %[Hz]
uVars.ao.sClkDcyc = 50;     %[%]
uVars.ao.fgClk    = 100e6;  %[S/s]

% Digitizer
uVars.ao.timeToSample   = 1; %[s]
uVars.ao.channels       = 1;     %[#]
uVars.ao.extCropDef     = false;
uVars.ao.extCropSamples = 0;
uVars.ao.measTimeLimit   = 1;

% Frequency
uVars.ao.envDC = 100e3;
uVars.ao.envUS = 50e3;

% Geometry & Length
uVars.ao.c               = 1500; %[m/s]
uVars.ao.distFromPhantom = 0;    %[m]

% General Operation
uVars.ao.useFrame        = false;
uVars.ao.useGPU          = true;
uVars.ao.useHadamard     = false;
uVars.ao.analyzeSingleCh = true;
uVars.ao.acCoupling      = true;
uVars.ao.keepSplits      = false;

uVars.ao.displayBuildUp  = false;
uVars.ao.displayBUEvery  = 16;

uVars.ao.useCalibration  = false;
uVars.ao.timeToSampleCal = 2; %[s]

uVars.ao.cutArtfct    = false;
uVars.ao.artfctIdxVec = 200:220;

uVars.ao.calcMuEff  = false;
uVars.ao.muEffModel = 'Uniform'; % Uniform, Point, TwoFibers
uVars.ao.muEffIdxs  = 200:220;

uVars.ao.noiseStd   = [];

uVars.ao.useVirtualData      = false;
uVars.ao.virtualDataNoiseSTD = 0;
uVars.ao.limitByN            = true;
uVars.ao.N                   = 10;
uVars.ao.dispTimeTable       = true; % relevant only with GUI

uVars.ao.skipParamsCheck = false;

uVars.ao.uploadToTelegram    = false;
uVars.ao.telegramChatID      = '-512325870';

% fGen
uVars.ao.usPower = 100; % [%]

% IO
uVars.ao.IOPort = 1;
uVars.ao.IOLine = 4;

% Export
uVars.ao.exportData.meas        = false;
uVars.ao.exportData.signal      = false;
uVars.ao.exportData.deMul       = false;
uVars.ao.exportData.reshaped    = false;
uVars.ao.exportData.fft         = false;
uVars.ao.exportData.usCompCmplx = false;

%----------------------------
% Figures Variables
%---------------------------
uVars.figs.ch        = 1;
uVars.figs.depthIdx  = 146;
uVars.figs.frame     = 1;

uVars.figs.reopenFigures  = false;
uVars.figs.displayFullFFT = true;
uVars.figs.FFTenv         = 250e3; %Hz

uVars.figs.singleChIdx = 1;
uVars.figs.depthAxType = 'Normal';
uVars.figs.intExt      = 'int';

uVars.figs.validStruct.extClk   = false;
uVars.figs.validStruct.usSignal = false;

uVars.figs.validStruct.cropped  = false;

uVars.figs.validStruct.meas     = false;
uVars.figs.validStruct.signal   = false;
uVars.figs.validStruct.deMul    = false;
uVars.figs.validStruct.reshaped = false;

uVars.figs.validStruct.rawFFT         = false;
uVars.figs.validStruct.calibration    = false;
uVars.figs.validStruct.fittingModel   = false;
uVars.figs.validStruct.fittedPowerFFT = false;
uVars.figs.validStruct.finalFFT       = true;

uVars.figs.validStruct.phi    = true;
uVars.figs.validStruct.rawPhi = false;
uVars.figs.validStruct.phiLog = false;
%----------------------%
% File System Variables
%----------------------%
uVars.fileSystem.saveMeas       = true;
uVars.fileSystem.saveSignal     = false;
uVars.fileSystem.saveDeMul      = false;
uVars.fileSystem.saveReshaped   = false;
uVars.fileSystem.saveFFT        = false;
uVars.fileSystem.savePhiChCmplx = false;

uVars.fileSystem.saveResults = true;
uVars.fileSystem.saveFigs    = false; %not implemented yet

uVars.fileSystem.dirPath     = "D:\Basics";
uVars.fileSystem.projName    = "PMT-NoFrame-NoWLScan-1Ch-Tmp";
uVars.fileSystem.resDirName  = "Results";

uVars.fileSystem.extProject     = false;
uVars.fileSystem.useExtVarsPath = false;
uVars.fileSystem.extVarsPath    = [];

ao.setVars(uVars);
ao.configPeripherals();
aoVars = ao.getVars();
res    = ao.runAcoustoOptics();
%%
algoVars = aoVars.measVars.algo;
signal   = res.export.signal;
reShaped = squeeze(res.export.reshaped);
figure(); 
subplot(2,2,1)
plot(signal)
subplot(2,2,2)
imagesc(reShaped)
subplot(2,2,3)
plot(reShaped(:,33)); hold on
plot(reShaped(:, 34))
subplot(2,2,4)
plot(reShaped(:,35)); hold on
plot(reShaped(:, 36))
plot(reShaped(:, 37))

%% Single Point AO
% %% Live Acousto Optics
% uVars.ao.limitByN            = true;
% uVars.ao.N                   = 3;
% uVars.fileSystem.saveRawData = false;
% uVars.fileSystem.projName = "Example Live AO";
% 
% ao.setMeasVars(uVars);
% ao.configPeripherals();
% aoVars = ao.getVars();
% res    = ao.liveAcoustoOptics();

%% Loaded Data
clc
path = "D:\MuEff\Uniform\GradedScatteringSet\Focused\Phantom-1";

data = ao.loadData(path);


%% Recalc Loaded Data
% close all;
uVarsReCalc = data.vars.uVarsNew;
uVarsReCalc.ao.acCoupling = true;
uVarsReCalc.ao.noiseStd = 6e-3;
uVarsReCalc.figs.depthAxType = 'Index';

uVarsReCalc.figs.ch        = 3;

uVarsReCalc.figs.validStruct.rawFFT         = true;
uVarsReCalc.figs.validStruct.calibration    = false;
uVarsReCalc.figs.validStruct.fittingModel   = true;
uVarsReCalc.figs.validStruct.fittedPowerFFT = true;
uVarsReCalc.figs.validStruct.finalFFT       = true;

uVarsReCalc.figs.validStruct.phi    = true;
uVarsReCalc.figs.validStruct.rawPhi = true;
uVarsReCalc.figs.validStruct.phiLog = true;

resNew = ao.reCalcData(uVarsReCalc);


% figure();
% plot(res.analysis.phi.phiLaplace)

% figure(); 
% subplot(2,2,1); 
% plot(phi)
% title("phi")
% subplot(2,2,2); 
% plot(normTypes.phiNorm3)
% title("phi Norm")
% subplot(2,2,3)
% plot(normTypes.phiLog3)


%%
close all;

algo = Algo();
algoVars = vars.extVars.algo;

algoVars.channels = 1;
algoVars.noiseStd = 6e-3;

algo.setVars(algoVars);

depthAxis = -1*(vars.measVars.algo.len.depthVec*1e3 - 104);
% depthAxis = (vars.measVars.algo.len.depthVec*1e3);
numOfPos  = vars.measVars.algo.samples.numOfPos;

hFig1 = figure();
set(hFig1, 'color', 'white');   
ax1 = subplot(1,2,1);
set(ax1, 'color', 'white');
hP1 = plot(depthAxis, zeros(1, numOfPos));
% hP1 = plot( zeros(1, numOfPos));
hTit1 = title("AO Reconstruction After 0s");
xlabel("Depth[mm]");
ylabel("Signal [Linear - AU]")
xlim([0, 100]);
ylim([0,1.01])
hText = text(55, 0.9, "SNR:");
ax2 = subplot(1,2,2);
set(ax2, 'color', 'white');
hP2 = plot(depthAxis, zeros(1, numOfPos));
hTit2 = title("AO Reconstruction After 0s");
xlabel("Depth[mm]");
ylabel("Signal [dB]")
xlim([0, 100]);
ylim([-70,1])

set(hFig1, 'Position', [197,428,1206,365]);

for i=1:length(res.splitRes)
    curRes = res.splitRes(i);
    curRes.frameAvgPowerFFT = curRes.frameAvgPowerFFT(1,:,:);
    if  i==1
        accRes = curRes;
        accRes.frameAvgPowerFFT = accRes.frameAvgPowerFFT(1,:,:);
    else
        accRes.frameAvgPowerFFT = (accRes.frameAvgPowerFFT*(i-1) +  curRes.frameAvgPowerFFT) /i;
    end

    resSingle(i) = algo.reconExtFFTData(curRes);
    phiSingle = normMatf(resSingle(i).phiNorm(1:104));
    SNRSingle(i) = 1/std(phiSingle(1:40));  

    resTotal(i)  = algo.reconExtFFTData(accRes);
    phiTotal  = normMatf(resTotal(i).phiNorm(1:104));
    SNRTotal  = 1/std(phiTotal(1:40)); 

    set(hP1, 'YData',  phiTotal);
    set(hP2, 'YData',  db(resTotal(i).phiNorm));
    set(hTit1, 'String', sprintf( "AO Reconstruction After %d[s]",i));
    set(hTit2, 'String', sprintf( "AO Reconstruction After %d[s]",i));
    set(hText, 'String', sprintf("SNR: %.2f", SNRTotal));
    drawnow();
    pause(0.01);
    frame(i) = getframe(hFig1);
end


skip = 1;
for i = 1:skip:128
    im = frame2im(frame(i)); 
    [imind,cm] = rgb2ind(im,256);
    if i == 1
      imwrite(imind, cm, "transPhantom.gif", 'gif', 'DelayTime', 0.5, 'Loopcount', inf); 
    else 
      imwrite(imind, cm, "transPhantom.gif", 'gif', 'DelayTime', 0.5, 'WriteMode', 'append'); 
    end 
end

mean(SNRSingle)





