close all;
clear all;

uVars = Consistency.uVarsCreate();

uVarsAcoustoOptics = acoustoOptics.uVarsCreate();

uVarsAcoustoOptics.c                 = 1500;
uVarsAcoustoOptics.fSin              = 1.25e6;              
uVarsAcoustoOptics.fTrain            = 19e3;
uVarsAcoustoOptics.cycInPulse        = 4; 
uVarsAcoustoOptics.channels          = 4; %update in digitizer
uVarsAcoustoOptics.phantomDepth      = 4.2e-2;
uVarsAcoustoOptics.distFromPhantom   = 0;
uVarsAcoustoOptics.fExtClk           = 5e6; %fs
uVarsAcoustoOptics.fSclk             = 100e6;     %update in fGen
uVarsAcoustoOptics.timeToSample      = 6;
uVarsAcoustoOptics.extClkDcyc        = 50; % [%]
uVarsAcoustoOptics.IOPort            = 1;
uVarsAcoustoOptics.IOLine            = 1;
uVarsAcoustoOptics.useHadamard       = false;
uVarsAcoustoOptics.fastAnalysis      = false;
uVarsAcoustoOptics.useGPU            = false; %algo, digitizer;
uVarsAcoustoOptics.exportRawData     = true;

uVarsAcoustoOptics.gReq.ch  = 4;
uVarsAcoustoOptics.gReq.pos = 1;
uVarsAcoustoOptics.gReq.intExt = 'int';
names = fieldnames(uVarsAcoustoOptics.gReq.validStruct);
for i=1:length(names)
    uVarsAcoustoOptics.gReq.validStruct.(names{i}) = true;
end

uVars.acoustoOptics = uVarsAcoustoOptics;

uVars.stages.xPos = 0;
uVars.stages.yPos = 0;

uVars.fileSystem.scanName = 'Test';
uVars.fileSystem.resDirPath = 'C:\Users\sahiadl.EED\OneDrive - Technion\Masters\AcoustoOpticSystem\Code\Tests';
uVars.fileSystem.saveFullData = false;
uVars.fileSystem.saveReducedData = false;
uVars.fileSystem.saveFigs = false;
uVars.fileSystem.saveResults  = false;

uVars.scan.startTime = 0.5;
uVars.scan.stride = 0.5;
uVars.scan.endTime = 2;
uVars.scan.useQuant = false;
uVars.scan.quantTime = 0.25;
uVars.scan.numOfSets = 4;

uVars.gReq.ch  = 1;
uVars.gReq.pos = 1;
uVars.gReq.intExt = 'int';
names = fieldnames(uVars.gReq.validStruct);
for i=1:length(names)
    uVars.gReq.validStruct.(names{i}) = true;
end

acoustoOptics = acoustoOptics();
acoustoOptics.init();
stages = stages('COM3');

acoustoOptics.setMeasVars(uVarsAcoustoOptics);
acoustoOptics.configPeripherals();
resAO = acoustoOptics.measureAndAnlayse();

consistency = Consistency(acoustoOptics, stages, []);
consistency.init();
consistency.startScan(uVars);