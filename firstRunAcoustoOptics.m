close all;
clear all;
clc;

%should affect fGen algo & digitizer

uVars = acoustoOptics.uVarsCreate();

uVars.c                 = 1500;
uVars.fSin              = 1.25e6;              
uVars.fTrain            = 19e3;
uVars.cycInPulse        = 4; 
uVars.channels          = 4; %update in digitizer
uVars.phantomDepth      = 4.2e-2;
uVars.distFromPhantom   = 0;
uVars.fExtClk           = 5e6; %fs
uVars.fSclk             = 100e6;     %update in fGen
uVars.timeToSample      = 0.25;
uVars.extClkDcyc        = 50; % [%]
uVars.IOPort            = 1;
uVars.IOLine            = 4;
uVars.useGPU            = true; %algo, digitizer;
uVars.exportRawData     = true;

uVars.gReq.ch  = 4;
uVars.gReq.pos = 1;
uVars.gReq.intExt = 'int';
names = fieldnames(uVars.gReq.validStruct);
for i=1:length(names)
    uVars.gReq.validStruct.(names{i}) = true;
end

acoustoOptics = acoustoOptics();

acoustoOptics.init();

acoustoOptics.setMeasVars(uVars);

% acoustoOptics.prepareGraphics();

acoustoOptics.configPeripherals();

tic
acoustoOptics.measureAndAnlayse();
toc

acoustoOptics.fullMeasureAndAnalyse(uVars);

timeTable = acoustoOptics.getTimeTable();

% acoustoOptics.fGen.disableAllOutputs()
% clear acoustoOptics



