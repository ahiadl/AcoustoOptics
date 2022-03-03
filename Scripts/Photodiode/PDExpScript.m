close all;
clear all;
clc;

pdexp = PDExp();
stagesT = pdexp.stagesT;


%% soft Init
% pdexp =  PDExp(ao, stagesT);
%% Config
uVars = pdexp.createUserVars();

%AO
uVars.frame           = 0.002;
uVars.timeToSample    = 2;
uVars.fs              = 20e6;
uVars.fSqnc           = 5e3;
uVars.distFromPhantom = 10.5e-2;
uVars.useGPU          = true;
uVars.idxToPlot       = 35;
uVars.useHadamard     = true;
uVars.exportMeas      = false;
% (Other AO variables are configurable via GUI)

uVars.fwl = 0.5;

uVars.usePhiPP    = true;
uVars.refSNR      = 80;
uVars.updateAOGui = true;

% Extract Net Signal
uVars.peakEnv = 29;
uVars.peakIdx = 34;
uVars.cutIdxs = [120:140, 165:185, 210:230];
uVars.bkgIdx  = 150;

% Power
uVars.netSignal = 660e-9;
uVars.minPwr    = 4.8e-6;
uVars.maxPwr    = 2.9e-3;

% Measurement
uVars.reps       = 30;
uVars.stabTime   = 0;
uVars.subSetSize = 10; %calculate dynamic average SNR on subset
uVars.measType   = 'Dual';

uVars.minTheta   = 56.25;
uVars.thetaTrend = 'Negative';

uVars.thetaVec  = [56.25, 50:-5:15, 11.25];
% uVars.thetaVec  = [];
% uVars.refFactor = [50, 100, 500];
uVars.refFactor = [];
uVars.voltage   = [575,600,625]; 
uVars.rotate    = true;

uVars.logScale  = true;

% File System
uVars.measName = 'PD-2sT-5KhzRR-5HzWL-20MHzFs-1msQ-20reps-SNR-vs-Ref-Fine';
uVars.dir      = 'D:\TMP';
uVars.saveFlag = false;
uVars.saveN    = 0;

uVars.tgUpdate = true;

pdexp.setVars(uVars);
pdExpVars = pdexp.getVars();

fprintf("Done Setting Vars\n")


% Run
res = pdexp.measure();
