close all;
clear all;
clc;

%%
instrreset
ao  = acoustoOptics();
ao.init();
stages = stages('PI');
stages.connect();
%%
aoS2D = scan2DAO(ao, stages);

%%
oa = optoAcoustic();

%%
Scope_visa_address = 'USB0::0x1AB1::0x04CE::DS1ZA172115085::0::INSTR';
scope = Rigol_DS1074z('ni',Scope_visa_address);
peUS =  pulseEchoUS(scope, stages);

%%
handles.aoS2D = aoS2D;
handles.stages = stages;
handles.oa = oa;
handles.peUS = peUS;
intSys = IntSys(handles);

%%
uVars = intSys.uVarsCreate();

uVars.grid.aoScanStart   = 70; 
uVars.grid.aoScanStride  = 1;
uVars.grid.aoScanEnd     = 80;

uVars.grid.oaScanStart  = 145;
uVars.grid.oaScanStride = 1;
uVars.grid.oaScanEnd    = 145;

uVars.grid.aoScanAxis = 'X'; % should be constantlu defined
uVars.grid.oaScanAxis = 'Z'; % should be constantlu defined
uVars.grid.depthAxis  = 'Y';

% Align the stages to requested directions
uVars.stages.order = 'XYZ';
uVars.stages.stagesOwnedByParent = false;

%-------------------------
% General Vars
%-------------------------
uVars.c              = 1500;   %[m/s]
uVars.imageWidth     = 30e-3;
uVars.distFromCenter = 9.2e-2;

%-------------------------
% Acousto Optics Vars
%-------------------------
uVars.aoS2D.ao.ao.fSin             = 1.25e6; %[Hz]              
uVars.aoS2D.ao.ao.fTrain           = 15.9e3; %[Hz]
uVars.aoS2D.ao.ao.cycInPulse       = 1;      %[#] 
uVars.aoS2D.ao.ao.channels         = 4; %[#]

uVars.aoS2D.ao.ao.fExtClk          = 5e6;    %[Hz]
uVars.aoS2D.ao.ao.extClkDcyc       = 50;     %[%]
uVars.aoS2D.ao.ao.fSclk            = 100e6;  %[S/s]

uVars.aoS2D.ao.ao.timeToSample     = 0.5; %[s]
uVars.aoS2D.ao.ao.useQuant         = true;
uVars.aoS2D.ao.ao.quantTime        = 0.002; %[s]

uVars.aoS2D.ao.ao.phantomDepth     = 4.6e-2; %[m]
uVars.aoS2D.ao.ao.distFromPhantom  = 7.3e-2; %[m]

uVars.aoS2D.ao.ao.envDC            = 100e3;
uVars.aoS2D.ao.ao.envUS            = 78e3;
uVars.aoS2D.ao.ao.envHar           = 200e3;

uVars.aoS2D.ao.ao.fastAnalysis     = false;
uVars.aoS2D.ao.ao.useHadamard      = true;
uVars.aoS2D.ao.ao.useGPU           = false;

uVars.aoS2D.ao.ao.useVirtualData   = false;
uVars.aoS2D.ao.ao.limitByN         = false;
uVars.aoS2D.ao.ao.N                = 10;
uVars.aoS2D.ao.ao.dispTimeTable    = true; % relevant only with GUI

uVars.aoS2D.ao.ao.IOPort = 1;
uVars.aoS2D.ao.ao.IOLine = 4;

%-------------------------
%Opto Acoustics Vars
%-------------------------

uVars.oa.activeChannels  = 1:256; %Mask
uVars.oa.numOfSamples    = 2030;
uVars.oa.numOfAvg        = 66;
uVars.oa.numOfFrames     = 1;

uVars.oa.fs = 40e6;
uVars.oa.geometry = 'circular';
uVars.oa.imageWidth = 30e-3; %[m]
uVars.oa.BPmode = 3;

%--------------------------
%Pulse Echo Ultrasound Vars
%--------------------------
uVars.peus.usFreq        = 1.25e6;
uVars.peus.signalChannel = 1;
uVars.peus.trigCh        = 2;
uVars.peus.avg           = 128;
uVars.peus.usRepRate     = 100;
uVars.peus.pulseIntVal   = 12;

%---------------------%
%Figures Vars
%---------------------%
uVars.aoS2D.figs.intExt        = 'int';
uVars.aoS2D.figs.reopenFigures = false;

uVars.figs.aoScanAxType   = 'Normal'; %Default: Normal [start : end]. Options: 'Center' [-span/2, span/2], 'Zero' [0, len], 'Index' [1:length(axis)]
uVars.figs.aoDepthAxType  = 'Normal'; 
uVars.figs.oaScanAxType   = 'Normal';

%Scan 2D Graphics
uVars.aoS2D.figs.depthIdx      = 20;

uVars.aoS2D.figs.validStruct.curMainAxis     = false;
uVars.aoS2D.figs.validStruct.curMainAxisAvg  = false;
uVars.aoS2D.figs.validStruct.curMainAxisLog  = false;
uVars.aoS2D.figs.validStruct.curMainPlane    = true;
uVars.aoS2D.figs.validStruct.curMainPlaneAvg = false;
uVars.aoS2D.figs.validStruct.curMainPlaneLog = false;
d
uVars.aoS2D.figs.normColorsToPlane = true;

% AO Graphics
uVars.aoS2D.ao.figs.depthIdx = 1;
uVars.aoS2D.ao.figs.ch       = 1;
uVars.aoS2D.ao.figs.quant    = 1;

uVars.aoS2D.ao.plotPhiInd          = false; 
uVars.aoS2D.ao.figs.displayFullFFT = true;
uVars.aoS2D.ao.figs.FFTenv         = 250e3;

%AO int/ext settings are dictated by the scan2D int/ext flag.
uVars.aoS2D.ao.figs.validStruct.extClk         = false;
uVars.aoS2D.ao.figs.validStruct.usSignal       = false;
uVars.aoS2D.ao.figs.validStruct.fullSignal     = false;
uVars.aoS2D.ao.figs.validStruct.measSamples    = false;
uVars.aoS2D.ao.figs.validStruct.netSignal      = false;
uVars.aoS2D.ao.figs.validStruct.deMul          = false;
uVars.aoS2D.ao.figs.validStruct.reshapedSignal = false;
uVars.aoS2D.ao.figs.validStruct.qAvgChFFT      = false;
uVars.aoS2D.ao.figs.validStruct.unFittedFFT    = false;
uVars.aoS2D.ao.figs.validStruct.fittedFFT      = false;
uVars.aoS2D.ao.figs.validStruct.phi            = true;
uVars.aoS2D.ao.figs.validStruct.rawPhi         = false;

uVars.oa.figs.validStruct.sinogram = true;
uVars.oa.figs.validStruct.recon    = true;

uVars.peus.figs.validStruct.rawData = true;
uVars.peus.figs.validStruct.recon   = true;
%----------------------%
% File System Variables
%----------------------%
uVars.fileSystem.saveResults = false;

uVars.aoS2D.ao.fileSystem.saveRawData        = false;
uVars.aoS2D.ao.fileSystem.saveNetSignal      = false;
uVars.aoS2D.ao.fileSystem.saveDemultiplexed  = false;
uVars.aoS2D.ao.fileSystem.saveReshapedSignal = false;
uVars.aoS2D.ao.fileSystem.saveFFT            = false;

%%
intSys.setUserVars(uVars);
intSys.configureScan();
%%
res = intSys.scan();