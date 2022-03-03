close all;
clear all;
clc;

ao  = acoustoOptics();
ao.init();
stages = stages('COM3');
stages.connect();
s2D = scan2DAO(ao, stages);
dv  = deepView(s2D);

%%
%-------------------------
% Scan 3D Variables
%-------------------------
uVars = dv.uVarsCreate();

uVars.grid.repeats   = 1;

uVars.grid.firstStart   = 102;
uVars.grid.firstStride  = -1;
uVars.grid.firstEnd     = 32;
uVars.grid.secondPos    = 86.5;

uVars.time.start  = 2;
uVars.time.stride = 198;
uVars.time.end    = 200;

% This will inform the stages in what order to scan
% and the grid will be calculated accordingly
uVars.grid.firstAxis  = 'X';
uVars.grid.secondAxis = 'Y';

% Align the stages to requested directions
uVars.stages.firstAxStageId  = 2;
uVars.stages.secondAxStageId = 1;

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
uVars.ao.ao.useHadamard      = true;
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
uVars.figs.intExt = 'int'; %Default: 'int', Optins: 'ext' must come with argument to extH struct.
uVars.figs.validStruct.deepView = true;

uVars.figs.depthIdx = 20;
uVars.figs.firstAxType = 'normal'; %Default: normal [start : end]. Options: 'cntr' [-span/2, span/2], 'zero' [0, len], ind [1:length(axis)]
uVars.figs.depthAxType = 'normal'; %Default: normal [start : end]. Options: 'cntr' [-span/2, span/2], 'zero' [0, len], ind [1:length(axis)]

uVars.figs.reopenFigures = false;

%Scan 2D Graphics
uVars.s2D.figs.validStruct.curMainAxis     = false;
uVars.s2D.figs.validStruct.curMainAxisAvg  = false;
uVars.s2D.figs.validStruct.curMainPlane    = true;
uVars.s2D.figs.validStruct.curMainPlaneAvg = false;

% AO Graphics
uVars.ao.figs.zIdx  = 20;
uVars.ao.figs.ch    = 1;
uVars.ao.figs.quant = 1;

uVars.ao.figs.plotPhiInd = true; 
uVars.ao.figs.displayFullFFT = true;
uVars.ao.figs.FFTenv = 250e3;

%AO int/ext settings are dictated by the scan2D int/ext flag.
uVars.ao.figs.validStruct.extClk         = false;
uVars.ao.figs.validStruct.usSignal       = false;
uVars.ao.figs.validStruct.fullSignal     = false;
uVars.ao.figs.validStruct.measSamples    = false;
uVars.ao.figs.validStruct.netSignal      = false;
uVars.ao.figs.validStruct.deMul          = false;
uVars.ao.figs.validStruct.reshapedSignal = false;
uVars.ao.figs.validStruct.qAvgChFFT      = true;
uVars.ao.figs.validStruct.unFittedFFT    = true;
uVars.ao.figs.validStruct.fittedFFT      = true;
uVars.ao.figs.validStruct.phi            = true;

%----------------------%
% File System Variables
%----------------------%
% scan3D Filesystem variables
uVars.fileSystem.saveResults  = true;
uVars.fileSystem.saveFigs     = false; %not implemented yet
uVars.fileSystem.dontSaveVars = false; % default: false

uVars.fileSystem.dirPath  = "D:\Results";
uVars.fileSystem.projName = "DV 2-200";

uVars.fileSystem.extProject     = false;

% Acousto Optics
uVars.ao.fileSystem.saveRawData        = false;
uVars.ao.fileSystem.saveNetSignal      = false;
uVars.ao.fileSystem.saveDemultiplexed  = false;
uVars.ao.fileSystem.saveReshapedSignal = false;
uVars.ao.fileSystem.saveFFT            = false;

%%
dv.setUserVars(uVars); 
dv.configureScan();
res = dv.scanDeepView();

%% navigate


%% load data


