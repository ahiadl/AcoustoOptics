close all;
clear all;
clc;

instrreset
ao  = acoustoOptics();
ao.init();
stages = stages('COM3');
stages.connect();
s2D = scan2DAO(ao, stages);
handles.s2D = s2D;
handles.stages = stages;
s3D = scan3DAO(handles);

%%
%-------------------------
% Scan 3D Variables
%-------------------------
uVars = s3D.uVarsCreate();

uVars.grid.repeats   = 1;

uVars.grid.firstStart   = 80;
uVars.grid.firstStride  = 1;
uVars.grid.firstEnd     = 90;

uVars.grid.secondStart  = 100;
uVars.grid.secondStride = -1;
uVars.grid.secondEnd    = 99;

% This will inform the stages in what order to scan
% and the grid will be calculated accordingly
uVars.grid.firstAxis  = 'Y';
uVars.grid.secondAxis = 'X';

% Align the stages to requested directions
uVars.stages.firstAxStageId      = 1;
uVars.stages.secondAxStageId     = 2;
uVars.stages.stagesOwnedByParent = false;
%-------------------------
% Algo and Peripheral Vars
%-------------------------
uVars.ao.ao.c                = 1500;   %[m/s]
uVars.ao.ao.fSin             = 1.25e6; %[Hz]              
uVars.ao.ao.fTrain           = 15.9e3; %[Hz]
uVars.ao.ao.cycInPulse       = 1;      %[#] 
uVars.ao.ao.channels         = 4; %[#]

uVars.ao.ao.fExtClk          = 5e6;    %[Hz]
uVars.ao.ao.extClkDcyc       = 50;     %[%]
uVars.ao.ao.fSclk            = 100e6;  %[S/s]

uVars.ao.ao.timeToSample     = 2; %[s]
uVars.ao.ao.useQuant         = true;
uVars.ao.ao.quantTime        = 0.002; %[s]

uVars.ao.ao.phantomDepth     = 4.6e-2; %[m]
uVars.ao.ao.distFromPhantom  = 7.3e-2; %[m]

uVars.ao.ao.envDC            = 100e3;
uVars.ao.ao.envUS            = 78e3;
uVars.ao.ao.envHar           = 200e3;

uVars.ao.ao.fastAnalysis     = false;
uVars.ao.ao.useHadamard      = false;
uVars.ao.ao.useGPU           = false;

uVars.ao.ao.useVirtualData   = false;
uVars.ao.ao.limitByN         = false;
uVars.ao.ao.N                = 10;
uVars.ao.ao.dispTimeTable    = true; % relevant only with GUI

uVars.ao.ao.IOPort = 1;
uVars.ao.ao.IOLine = 4;

%---------------------%
%Figures Vars
%---------------------%
% scan3D Graphics
uVars.figs.depthIdx = 19;
uVars.figs.intExt   = 'int';
uVars.figs.validStruct.nav    = true;
uVars.figs.validStruct.navAvg = true;
uVars.figs.normColorsToPlane  = false;
uVars.figs.reopenFigures      = false;

uVars.figs.firstAxType  = 'Normal'; %Default: normal [start : end]. Options: 'cntr' [-span/2, span/2], 'zero' [0, len]
uVars.figs.depthAxType  = 'Normal'; %Default: normal [start : end]. Options: 'cntr' [-span/2, span/2], 'zero' [0, len]
uVars.figs.secondAxType = 'Normal';

%Scan 2D Graphics
uVars.s2D.figs.validStruct.curMainAxis     = true;
uVars.s2D.figs.validStruct.curMainAxisAvg  = true;
uVars.s2D.figs.validStruct.curMainPlane    = true;
uVars.s2D.figs.validStruct.curMainPlaneAvg = true;

% AO Graphics
uVars.ao.figs.zIdx  = 1;
uVars.ao.figs.ch    = 1;
uVars.ao.figs.quant = 1;

uVars.ao.plotPhiInd          = false; 
uVars.ao.figs.displayFullFFT = true;
uVars.ao.figs.FFTenv         = 250e3;

%AO int/ext settings are dictated by the scan2D int/ext flag.
uVars.ao.figs.validStruct.extClk         = false;
uVars.ao.figs.validStruct.usSignal       = false;
uVars.ao.figs.validStruct.fullSignal     = false;
uVars.ao.figs.validStruct.measSamples    = false;
uVars.ao.figs.validStruct.netSignal      = false;
uVars.ao.figs.validStruct.deMul          = false;
uVars.ao.figs.validStruct.reshapedSignal = false;
uVars.ao.figs.validStruct.qAvgChFFT      = true;
uVars.ao.figs.validStruct.unFittedFFT    = false;
uVars.ao.figs.validStruct.fittedFFT      = true;
uVars.ao.figs.validStruct.phi            = true;

%----------------------%
% File System Variables
%----------------------%
% scan3D Filesystem variables
uVars.fileSystem.saveResults  = false;
uVars.fileSystem.saveFigs     = false; %not implemented yet
uVars.fileSystem.dontSaveVars = false; % default: false

uVars.fileSystem.dirPath  = ".\Scripts\objectUseExamples\examplesData";
uVars.fileSystem.projName = "Example S3D";

uVars.fileSystem.extProject = false;

% Acousto Optics
uVars.ao.fileSystem.saveRawData        = false;
uVars.ao.fileSystem.saveNetSignal      = false;
uVars.ao.fileSystem.saveDemultiplexed  = false;
uVars.ao.fileSystem.saveReshapedSignal = false;
uVars.ao.fileSystem.saveFFT            = false;

%%
navRep = 1;
navAx = 'Z';
navIdx = 29;

s3D.setUserVars(uVars);
s3D.configureScan();
s3D.setNavVars(navAx, navIdx, navRep);
res = s3D.scan3D();

%% navigate


%% load data


