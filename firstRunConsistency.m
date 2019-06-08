close all;
clear all;
clc;

%should affect fGen algo & digitizer
uVars.acoustoOptics.fGen.ch{1}.amp = 1;
uVars.acoustoOptics.fGen.ch{1}.bias = 0;
uVars.acoustoOptics.fGen.ch{1}.Sclk = 100e6;
uVars.acoustoOptics.fGen.ch{1}.daqFilter = 50;

uVars.acoustoOptics.fGen.ch{2}.amp = 1;
uVars.acoustoOptics.fGen.ch{2}.bias = 0;
uVars.acoustoOptics.fGen.ch{2}.Sclk = 100e6;
uVars.acoustoOptics.fGen.ch{2}.daqFilter = 50;

uVars.acoustoOptics.algo.fSin = 1.25e6;              
uVars.acoustoOptics.algo.fTrain = 19e3;
uVars.acoustoOptics.algo.cycInPulse = 4;
uVars.acoustoOptics.algo.channels = 4;

%should affect digitizer & algo
uVars.acoustoOptics.algo.c = 1550;  
uVars.acoustoOptics.algo.phantomDepth = 4.2e-2;
uVars.acoustoOptics.algo.distFromPhantom = 0;
uVars.acoustoOptics.algo.bufferSizeBytes = 8*2^20;
uVars.acoustoOptics.algo.bytesPerSample = 2;
uVars.acoustoOptics.algo.preTriggerSamples = 14;
uVars.acoustoOptics.algo.fExtClk = 5e6; %fs
uVars.acoustoOptics.algo.fSclk   = 100e6;
uVars.acoustoOptics.algo.timeToSample = 20;
uVars.acoustoOptics.algo.extClkDcyc = 50; % [%]

uVars.acoustoOptics.digitizer.mode = 'TS';
uVars.acoustoOptics.digitizer.channels = 4;
uVars.acoustoOptics.digitizer.preTriggerSamples = 14;
uVars.acoustoOptics.digitizer.voltsRange = 2;

%should affect IO
uVars.acoustoOptics.IO.port = 1;
uVars.acoustoOptics.IO.line = 4;
uVars.acoustoOptics.IO.mode = 0; %outputOnly

uVars.acoustoOptics.graphics.internal = true;
uVars.acoustoOptics.graphics.ch  = 1;
uVars.acoustoOptics.graphics.pos = 1;
uVars.acoustoOptics.graphics.mask = ones(1,9);

uVars.acoustoOptics.general.useGPU = true;
uVars.acoustoOptics.general.exportRawData = false;

uVars.stages.xPos = 126;
uVars.stages.yPos = 75;

uVars.fileSystem.scanName = 'Test';
uVars.fileSystem.resDirPath = 'C:\Users\sahiadl.EED\OneDrive - Technion\Masters\AcoustoOpticSystem\Code\Tests';
uVars.fileSystem.saveFullData = true;
uVars.fileSystem.saveReducedData = false;
uVars.fileSystem.saveFigs = true;

uVars.scan.startTime = 2;
uVars.scan.stride = 1;
uVars.scan.endTime = 4;
uVars.scan.useQuant = true;
uVars.scan.speckleTime = 0.25;
uVars.scan.numOfSets = 4;

consistency = consistency();
consistency.init();
consistency.startScan(uVars);



% acoustoOptics = acoustoOptics();
% 
% acoustoOptics.init();
% 
% acoustoOptics.setMeasVars(uVars);
% 
% acoustoOptics.prepareGraphics();
% 
% acoustoOptics.configPeripherals();
% 
% tic
% acoustoOptics.measureAndAnlayse();
% toc
% 
% timeTable = acoustoOptics.getTimeTable();
% clear acoustoOptics



