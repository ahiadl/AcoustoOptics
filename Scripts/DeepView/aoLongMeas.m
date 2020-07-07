close all;
clear all;
clc;

ao = acoustoOptics();
ao.init();
stages = stages('COM3');
stages.connect();
stages.assignStagesAxes(['X', 'Y'], [2,1]);
stages.moveStageAxisAbs('X', 100)
stages.moveStageAxisAbs('Y', 85)
%%
uVars = ao.uVarsCreate();

%-------------------------
% Algo and Peripheral Vars
%-------------------------
uVars.ao.c                 = 1500;   %[m/s]
uVars.ao.fSin              = 1.25e6; %[Hz]              
uVars.ao.fTrain            = 15.9e3; %[Hz]
uVars.ao.cycInPulse        = 1;      %[#] 
uVars.ao.channels          = 4; %[#]

uVars.ao.fExtClk           = 5e6;    %[Hz]
uVars.ao.extClkDcyc        = 50;     %[%]
uVars.ao.fSclk             = 100e6;  %[S/s]

uVars.ao.timeToSample      = 2; %[s]
uVars.ao.useQuant          = true;
uVars.ao.quantTime         = 0.002; %[s]

uVars.ao.phantomDepth      = 4.6e-2; %[m]
uVars.ao.distFromPhantom   = 7.2e-2; %[m]
uVars.ao.envDC             = 100e3;
uVars.ao.envUS             = 78e3;
uVars.ao.envHar            = 200e3;

uVars.ao.fastAnalysis      = false;
uVars.ao.useHadamard       = false;
uVars.ao.useGPU            = false;

uVars.ao.useVirtualData   = false;
uVars.ao.limitByN         = false;
uVars.ao.N                = 10;
uVars.ao.dispTimeTable    = true; % relevant only with GUI

uVars.ao.IOPort = 1;
uVars.ao.IOLine = 4;

%----------------------------
% Figures Variables
%---------------------------
uVars.figs.ch    = 1;
uVars.figs.zIdx  = 20;
uVars.figs.quant = 1;

uVars.figs.reopenFigures = true;
uVars.figs.plotPhiInd = true; 
uVars.figs.displayFullFFT = true;
uVars.figs.FFTenv = 250e3; %Hz

uVars.figs.intExt = 'int';

uVars.figs.validStruct.extClk         = false;
uVars.figs.validStruct.usSignal       = false;
uVars.figs.validStruct.fullSignal     = false;
uVars.figs.validStruct.measSamples    = false;
uVars.figs.validStruct.netSignal      = false;
uVars.figs.validStruct.deMul          = false;
uVars.figs.validStruct.reshapedSignal = false;
uVars.figs.validStruct.qAvgChFFT      = false;
uVars.figs.validStruct.unFittedFFT    = true;
uVars.figs.validStruct.fittedFFT      = false;
uVars.figs.validStruct.phi            = true;

%----------------------%
% File System Variables
%----------------------%
uVars.fileSystem.saveRawData        = false;
uVars.fileSystem.saveNetSignal      = false;
uVars.fileSystem.saveDemultiplexed  = false;
uVars.fileSystem.saveReshapedSignal = false;
uVars.fileSystem.saveFFT            = false;
uVars.fileSystem.savePhiChCmplx     = false;

uVars.fileSystem.saveResults = false;
uVars.fileSystem.saveFigs    = false; %not implemented yet

uVars.fileSystem.dirPath     = ".\Scripts\objectUseExamples\examplesData";
uVars.fileSystem.projName    = "Example AO";
uVars.fileSystem.resDirName  = "Results";

uVars.fileSystem.extProject     = false;
uVars.fileSystem.useExtVarsPath = false;
uVars.fileSystem.extVarsPath    = [];

%% Single Point AO
uVars.ao.timeToSample      = 2;

ao.setMeasVars(uVars);
ao.configPeripherals();
aoVars = ao.getVars();
res2    = ao.runAcoustoOptics();

uVars.ao.timeToSample      = 10;

ao.setMeasVars(uVars);
ao.configPeripherals();
aoVars = ao.getVars();
res10    = ao.runAcoustoOptics();

idx = 20;
tailStart = 37;

SNR2 = res2.phi(idx)/std(res2.phi(tailStart:end));
SNR10 = res10.phi(idx)/std(res10.phi(tailStart:end));

factorSNR = SNR10/SNR2;


%% Live Acousto Optics
uVars.ao.limitByN                   = true;
uVars.ao.N                          = 10;
uVars.fileSystem.saveRawData        = false;
uVars.fileSystem.projName = "Example Live AO";

uVars.ao.timeToSample      = 2;

ao.setMeasVars(uVars);
ao.configPeripherals();
aoVars = ao.getVars();
res2x10    = ao.liveAcoustoOptics();

uVars.ao.timeToSample      = 10;

ao.setMeasVars(uVars);
ao.configPeripherals();
aoVars = ao.getVars();
res10x10    = ao.liveAcoustoOptics();
%%
idx = 20;
tailStart = 37;

for i=1:10
    STD2x10(i)  = mean(res2x10(i).phi(tailStart:end).^2);
    SNR2x10(i)  = res2x10(i).phi(idx)^2/ STD2x10(i);
    STD10x10(i) = mean(res10x10(i).phi(tailStart:end).^2);
    SNR10x10(i) = res10x10(i).phi(idx)^2/STD10x10(i);
    stdFactor(i) = STD2x10(i)/STD10x10(i);
end

SNR2x10tot = mean(SNR2x10);
SNR10x10tot = mean(SNR10x10);

STD2x10tot = mean(STD2x10);
STD10x10tot = mean(STD10x10);

factorSNRx10 = SNR10x10tot/SNR2x10tot;
factorSTDx10 = mean(stdFactor);
%%
uVars.ao.limitByN                   = true;
uVars.ao.N                          = 10;
uVars.fileSystem.saveRawData        = false;
uVars.fileSystem.projName = "Example Live AO";

uVars.ao.timeToSample      = 20;

ao.setMeasVars(uVars);
ao.configPeripherals();
aoVars = ao.getVars();
res20x10    = ao.liveAcoustoOptics();
%%
idx = 20;
tailStart = 37;

for i=1:10
    STD2x10(i)  = std(res2x10(i).phi(tailStart:end).^2);
    SNR2x10(i)  = res2x10(i).phi(idx)^2/ STD2x10(i);
    STD20x10(i) = std(res20x10(i).phi(tailStart:end).^2);
    SNR20x10(i) = res20x10(i).phi(idx)^2/STD20x10(i);
    stdFactor(i) = STD2x10(i)/STD20x10(i);
end

SNR2x10tot = mean(SNR2x10);
SNR20x10tot = mean(SNR20x10);

STD2x10tot = mean(STD2x10);
STD20x10tot = mean(STD20x10);

factorSNRx10 = SNR20x10tot/SNR2x10tot;
factorSTDx10 = mean(stdFactor);

%%
uVars.ao.limitByN                   = true;
uVars.ao.N                          = 5;
uVars.fileSystem.saveRawData        = false;
uVars.fileSystem.projName = "Example Live AO";

uVars.ao.timeToSample      = 200;

ao.setMeasVars(uVars);
ao.configPeripherals();
aoVars = ao.getVars();
res200x10 = ao.liveAcoustoOptics();
res200x5  = res200x10;
%%
idx = 20;
tailStart = 37;

for i=1:5
    STD2x5(i)    = std(res2x10(i).phi(tailStart:end).^2);
    SNR2x5(i)    = res2x10(i).phi(idx)^2/ STD2x5(i);
    STD200x5(i)  = std(res200x5(i).phi(tailStart:end).^2);
    SNR200x5(i)  = res200x5(i).phi(idx)^2/STD200x5(i);
    stdFactor(i) = STD2x5(i)/STD200x5(i);
end

SNR2x5tot   = mean(SNR2x5);
SNR200x5tot = mean(SNR200x5);

STD2x5tot   = mean(STD2x5);
STD200x5tot = mean(STD200x5);

factorSNRx5 = SNR200x5tot/SNR2x5tot;
factorSTDx5 = mean(stdFactor);

%% Loaded Data
uVars.fileSystem.saveResults = false;
uVars.fileSystem.saveFigs    = false;

rawDataPath = ".\Scripts\objectUseExamples\examplesData\07-May-2020 19-47-42-Example AO\AO-Results.mat";
varsPath    = ".\Scripts\objectUseExamples\examplesData\07-May-2020 19-47-42-Example AO\AO-Vars.mat";

rawData     = ao.loadRawDataToAO(rawDataPath);
vars        = load(varsPath);

%Important so no additional saving will occur
vars.uVars.fileSystem.saveRawData        = false;
vars.uVars.fileSystem.saveNetSignal      = false;
vars.uVars.fileSystem.saveDemultiplexed  = false;
vars.uVars.fileSystem.saveReshapedSignal = false;
vars.uVars.fileSystem.saveFFT            = false;
vars.uVars.fileSystem.savePhiChCmplx     = false;

vars.uVars.fileSystem.saveResults = false;
vars.uVars.fileSystem.saveFigs    = false; 

ao.setMeasVars(vars.uVars);
loadedDataRes  = ao.analyseLoadedData();



