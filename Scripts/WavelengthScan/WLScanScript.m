close all;
% clear all;
clc;

pdexp = PDExp(ao);
stagesT = pdexp.stagesT;
%% Config
uVars = pdexp.createUserVars();

%AO
uVars.quant           = 0.002;
uVars.timeToSample    = 32;
uVars.distFromPhantom = 6.5e-2;
uVars.fs              = 20e6;
uVars.useGPU          = true;
uVars.exportRawData   = true;
uVars.idxToPlot       = 40;
uVars.useHadamard     = false;

uVars.fwl = 0.5;

uVars.usePhiPP = false;
uVars.refSNR   = 144.69;
% Extract Net Signal
uVars.peakEnv = 22;
uVars.peakIdx = 38;
uVars.cutIdxs = [100:116, 136:157, 188:204 ];
uVars.bkgIdx = 78;
% Power
uVars.netSignal = 92e-9;
uVars.minPwr    = 1.66e-6 ;
uVars.maxPwr    = 2.8e-3;

% Measurement
uVars.reps       = 250;
uVars.stabTime   = 0;
uVars.subSetSize = 50;
uVars.measType   = 'Point';

uVars.minTheta   = 58.8;
uVars.thetaTrend = 'Negative';

uVars.thetaVec  = [58.8, 55, 50, 45, 40, 35, 30, 25];
uVars.refFactor = [];
uVars.refFactor = [];
uVars.voltage   = [575,600,625]; 
uVars.rotate    = true;

uVars.logScale  = true;

% File System
uVars.measName = 'PD-2s-CW';
uVars.dir      = 'D:\Photodiode\InStabilityCheck';
uVars.saveFlag = true;
uVars.saveN    = 0;

uVars.tgUpdate = true;

pdexp.setVars(uVars);
pdExpVars = pdexp.getVars();

fprintf("Done Setting Vars\n")
%%
% fsin = pdExpVars.ao.measVars.algo.freq.fUS;
% samplesPerTrain = pdExpVars.ao.measVars.algo.usSignal.SclkSamplesInTrain;
% dt  = 1/pdExpVars.ao.measVars.algo.extClk.fSclk;
% 
% tVec = dt*(0:1:(samplesPerTrain-1));
% sig = sin(2*pi*fsin.*tVec);
% 
% % figure();
% % plot(tVec, sig)
% 
% clkData = pdExpVars.ao.measVars.algo.extClk.data;
% 
% dataCh1 = ao.fGen.getChDataStruct();
% dataCh2 = ao.fGen.getChDataStruct();
% 
% dataCh1.data = sig;
% dataCh2.data = clkData;
% dataCh1.dataLen = length(dataCh1.data);
% dataCh2.dataLen = length(dataCh2.data);
% 
% ao.measVars.fGen = ao.extVars.fGen; 
% 
% ao.fGen.reset();
% ao.fGen.setProperties(pdExpVars.ao.extVars.fGen.ch{1}, pdExpVars.ao.extVars.fGen.ch{2});
% ao.fGen.setData(dataCh1, dataCh2);
% ao.fGen.configChannel(1);
% ao.fGen.configChannel(2);



%% Run
res = pdexp.measure();