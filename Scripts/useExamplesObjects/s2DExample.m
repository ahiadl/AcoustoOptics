close all;
clear all;
clc;

instrreset
ao  = acoustoOpticsNew();
ao.init();
stages = stages('Zaber', 'COM3');
stages.connect();
s2D = scan2DAONew(ao, stages);

%%
ao.setMeasLimit(0.2)

%-------------------------
% Scan 2D Variables
%-------------------------
uVars = s2D.createUserVars();

uVars.general.externalScan = false; % this means that scan2D is independant and dont have any parent object

uVars.grid.repeats = 2;

uVars.grid.scanStart   = 75;
uVars.grid.scanStride  = 1;
uVars.grid.scanEnd     = 80;

% The position on the non scanning axis.
uVars.grid.secondPos = 75;
uVars.grid.thirdPos  = 75;

% This will inform the stages in what order to scan
% and the grid will be calculated accordingly
uVars.grid.scanAxis   = 'X';
uVars.grid.secondAxis = 'Y'; %TODO: support 3D stages
uVars.grid.thirdAxis  = 'Z';

uVars.grid.depthAxis = 'Z';

uVars.stages.moveNonScanStages     = true;
uVars.stages.stagesOwnedByParent   = false;

uVars.general.externalScan = false;

%-------------------------
% Algo and Peripheral Vars
%-------------------------
% US Signal
uVars.ao.ao.cycPerPulse       = 1;      %[#] 
uVars.ao.ao.fSin              = 1.25e6; %[Hz]              
uVars.ao.ao.fSqnc             = 10e3; %[Hz]
uVars.ao.ao.frameTime         = 0.002; %[s]

% Sampling Clk
uVars.ao.ao.fs                = 5e6;    %[Hz]
uVars.ao.ao.sClkDcyc          = 50;     %[%]
uVars.ao.ao.fgClk             = 100e6;  %[S/s]

% Digitizer
uVars.ao.ao.timeToSample      = 0.4; %[s]
uVars.ao.ao.channels          = 2; %[#]

% Frequency
uVars.ao.ao.envDC             = 100e3;
uVars.ao.ao.envUS             = 78e3;

% Geometry & Length
uVars.ao.ao.c                 = 1500;   %[m/s]
uVars.ao.ao.distFromPhantom   = 0; %[m]

% General Operation
uVars.ao.ao.useFrame            = true;
uVars.ao.ao.useGPU              = true;
uVars.ao.ao.useHadamard         = true;
uVars.ao.ao.contHadamard        = false;
uVars.ao.ao.highResAO           = false;
uVars.ao.ao.analyzeSingleCh     = true;
uVars.ao.ao.contSpeckleAnalysis = false;

uVars.ao.ao.useVirtualData   = false;
uVars.ao.ao.limitByN         = true;
uVars.ao.ao.N                = 10;
uVars.ao.ao.dispTimeTable    = true; % relevant only with GUI

% IO
uVars.ao.ao.IOPort = 1;
uVars.ao.ao.IOLine = 4;

% Export
uVars.ao.ao.exportData.meas        = true;
uVars.ao.ao.exportData.signal      = true;
uVars.ao.ao.exportData.deMul       = true;
uVars.ao.ao.exportData.reshaped    = true;
uVars.ao.ao.exportData.fft         = true;
uVars.ao.ao.exportData.usCompCmplx = true;

%---------------------%
%Figures Vars
%---------------------%
%Scan 2D Graphics
uVars.figs.depthIdx      = 30;
uVars.figs.intExt        = 'int';
uVars.figs.useExtClims   = false;
uVars.figs.reopenFigures = false;

uVars.figs.firstAxType = 'Normal'; %Default: Normal [start : end]. Options: 'Center' [-span/2, span/2], 'Zero' [0, len], 'Index' [1:length(axis)]
uVars.figs.depthAxType = 'Normal';

uVars.figs.validStruct.curMainAxis     = true;
uVars.figs.validStruct.curMainAxisAvg  = true;
uVars.figs.validStruct.curMainPlane    = true;
uVars.figs.validStruct.curMainPlaneAvg = true;

uVars.figs.sepChIdx = 1;

% AO Graphics
uVars.figs.ch        = 1;
uVars.figs.depthIdx  = 30;
uVars.figs.frame     = 1;

uVars.figs.reopenFigures  = false;
uVars.figs.plotPhiInd     = true; 
uVars.figs.displayFullFFT = true;
uVars.figs.FFTenv         = 250e3; %Hz

uVars.figs.validStruct.extClk   = false;
uVars.figs.validStruct.usSignal = false;

uVars.figs.validStruct.meas     = false;
uVars.figs.validStruct.signal   = false;
uVars.figs.validStruct.deMul    = false;
uVars.figs.validStruct.reshaped = false;

uVars.figs.validStruct.unFittedFFTAllCh    = false;
uVars.figs.validStruct.unFittedFFTSingleCh = false;
uVars.figs.validStruct.fittedFFTSingleCh   = false;
uVars.figs.validStruct.fittedFFTAllCh      = false;
uVars.figs.validStruct.avgFFT              = false;
uVars.figs.validStruct.normFFT             = false;

uVars.figs.validStruct.phi    = true;
uVars.figs.validStruct.rawPhi = false;

%----------------------%
% File System Variables
%----------------------%
% Scan 2D FileSystem Variables
uVars.fileSystem.saveResults  = false;
uVars.fileSystem.saveFigs     = false; %not implemented yet
uVars.fileSystem.dontSaveVars = false; % default: false

uVars.fileSystem.dirPath     = ".\Scripts\objectUseExamples\examplesData";
uVars.fileSystem.projName    = "Example S2D";

uVars.fileSystem.extProject  = false;

% Acousto Optics
uVars.ao.fileSystem.saveMeas       = false;
uVars.ao.fileSystem.saveSignal     = false;
uVars.ao.fileSystem.saveDeMul      = false;
uVars.ao.fileSystem.saveReshaped   = false;
uVars.ao.fileSystem.saveFFT        = false;
uVars.ao.fileSystem.savePhiChCmplx = false;

%%
% close all;

s2D.setUserVars(uVars);
s2D.configureScan();
res = s2D.scan2D();
