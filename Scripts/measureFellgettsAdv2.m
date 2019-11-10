close all;
clear all;
clc;

%----------------------------------------------------------%
% In this file: Trying to understand why results from measureFellgettsAdv.m
% didn't correspond to theory.
% 1. init Acousto Optics element
% 2. Comparing with ultrasound naive and multiplexed analysis
% 3. Comparing without ultrasound naive and multiplexed analysis
%----------------------------------------------------------%

ao = acoustoOptics();
ao.init();

uVars = acoustoOptics.uVarsCreate();

uVars.c                 = 1445;
uVars.fSin              = 1.25e6;              
uVars.fTrain            = 16.5e3;
uVars.cycInPulse        = 4; 
uVars.channels          = 4; %update in digitizer
uVars.phantomDepth      = 4.6e-2;
uVars.distFromPhantom   = 7.15e-2;
uVars.fExtClk           = 5e6; %fs
uVars.fSclk             = 100e6;     %update in fGen
uVars.timeToSample      = 2;
uVars.extClkDcyc        = 50; % [%]
uVars.IOPort            = 1;
uVars.IOLine            = 4;
uVars.useGPU            = true; %algo, digitizer;
uVars.fastAnalysis      = false;
uVars.useQuant          = true;
uVars.quantTime         = 0.002;
uVars.useHadamard       = true;

uVars.gReq.ch  = 4;
uVars.gReq.pos = 1;
uVars.gReq.intExt = 'int';
names = fieldnames(uVars.gReq.validStruct);
for i=1:length(names)
    uVars.gReq.validStruct.(names{i}) = false;
end

uVars.gReq.validStruct.usSignal = true;
uVars.gReq.validStruct.phi = true;


uVars.exportRawData.rawData       = false;
uVars.exportRawData.netSignal     = false;
uVars.exportRawData.deMultiplexed = false;
uVars.exportRawData.reshape       = false;
uVars.exportRawData.FFT           = false;

uVars.fs.scanName        = [];
uVars.fs.resDirPath      = 'D:\Results';
uVars.fs.saveFullData    = false;
uVars.fs.saveFigs        = false;
uVars.fs.saveResults     = false;

repPerConf = 30;
%%

%no Hadamard:
uVars.useHadamard = false;
uVars.useGPU = true;
name = sprintf("Comparison-woHad-withUltrasound-4Cyc-17n");
uVars.fs.scanName = name;

ao.setMeasVars(uVars);
numOfPosNaive = ao.algo.samples.numOfPos;
ao.configPeripherals();

fprintf(name);
for k = 1:repPerConf
    resNaiveUSon{k} = ao.measureAndAnlayse();
    disp(k);
    pause(0.01);
end

%with Hadamard:
uVars.useHadamard = true;
uVars.useGPU = false;
name = sprintf("Comparison-wHad-withUltrasound-4Cyc-17n");
uVars.fs.scanName = name;
fprintf(name);

ao.setMeasVars(uVars);
numOfPosHad = ao.algo.samples.numOfPos;
ao.configPeripherals();

fprintf(name);
for k = 1:repPerConf
    resHadUSon{k} = ao.measureAndAnlayse();
    disp(k);
    pause(0.01);
end

%%
f = msgbox('Disconnect the transducer','Done');
waitfor(f)
fprintf("Runnins no ultrasound measurement");

%%
%no Hadamard:
uVars.useHadamard = false;
uVars.useGPU = true;
name = sprintf("Comparison-woHad-withUltrasound-4Cyc-17n");
uVars.fs.scanName = name;

ao.setMeasVars(uVars);
numOfPosNaive = ao.algo.samples.numOfPos;
ao.configPeripherals();

fprintf(name);
for k = 1:repPerConf
    resNaiveUSoff{k} = ao.measureAndAnlayse();
    disp(k);
    pause(0.01);
end

%with Hadamard:
uVars.useHadamard = true;
uVars.useGPU = false;
name = sprintf("Comparison-wHad-withUltrasound-4Cyc-17n");
uVars.fs.scanName = name;
fprintf(name);

ao.setMeasVars(uVars);
numOfPosHad = ao.algo.samples.numOfPos;
ao.configPeripherals();

fprintf(name);
for k = 1:repPerConf
    resHadUSoff{k} = ao.measureAndAnlayse();
    disp(k);
    pause(0.01);
end

%% Analysis

resNaiveUSon
resHadUSon

resNaiveUSoff
resHadUSoff

figure()
plot(resNaiveUSon{1}.phi); hold on;
plot(resHadUSon{1}.phi); hold on;

figure()
plot(resNaiveUSoff{1}.phi); hold on;
plot(resHadUSoff{1}.phi); hold on;

for i=1:30
    phiNaiveUSoff(i,:) = resNaiveUSoff{i}.phi;
    phiHadUSoff(i,:)   = resHadUSoff{i}.phi;
    phiNaiveUSon(i,:)  = resNaiveUSon{i}.phi;
    phiHadUSon(i,:)    = resHadUSon{i}.phi;
end


phiNaiveUSoffSTD = std(phiNaiveUSoff,0,1);
phiHadUSoffSTD   = std(phiHadUSoff,0,1);
phiNaiveUSonSTD  = std(phiNaiveUSon,0,1);
phiHadUSonSTD    = std(phiHadUSon,0,1);


figure();
plot(phiNaiveUSonSTD./phiHadUSonSTD); hold on
plot(phiNaiveUSoffSTD./phiHadUSoffSTD); hold on 
plot(resHadUSon{1}.phi*100);

save("./matFiles/FellgettsAdv2.mat", '-v7.3')


