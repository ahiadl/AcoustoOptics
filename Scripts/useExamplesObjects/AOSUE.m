close all;
clear all;
clc;

instrreset;
ao  = acoustoOptics();
ao.init();
stages = stages('Zaber', 'COM3');
stages.connect();

%%
s2D = scanAO(ao, stages);

%%
%-------------------------
% Scan 2D Variables
%-------------------------
uVars = s2D.createUserVars();

% Grid scan variables
uVars.grid.scan1Start  = 59.3256;
uVars.grid.scan1Stride = 0.3372;
uVars.grid.scan1End    = 68.1456;

uVars.grid.scan2Start  = 62.38;
uVars.grid.scan2Stride = 0;
uVars.grid.scan2End    = 0;

uVars.grid.scan3Start  = 111.65;
uVars.grid.scan3Stride = 0;
uVars.grid.scan3End    = 0;

uVars.grid.scan1Label  = 'X';
uVars.grid.scan2Label  = 'Y';
uVars.grid.scan3Label  = 'Z';
uVars.grid.depthLabel  = 'D';

uVars.grid.scanType   = 'Bi';

% General scan variables
uVars.general.stagesOwnedByParent = true;
uVars.general.moveNonScanStages = false;
uVars.general.externalScan = false;
uVars.general.keepUSOn     = false;
uVars.general.tgSend       = false;
uVars.general.chatID       = '-512325870';

%-------------------------
% Algo and Peripheral Vars
%-------------------------
% US Signal
uVars.ao.ao.cycPerPulse = 1;
uVars.ao.ao.fSin        = 1e6; 
uVars.ao.ao.fSqnc       = 4.07e3;
uVars.ao.ao.frameTime   = 0.05; 
uVars.ao.ao.usPower     = 100;

% Sampling Clk
uVars.ao.ao.fs                = 20e6;  %[Hz]
uVars.ao.ao.sClkDcyc          = 50;    %[%]
uVars.ao.ao.fgClk             = 100e6; %[S/s]

% Digitizer
uVars.ao.ao.timeToSample      = 0.4; %[s]
uVars.ao.ao.channels          = 2;   %[#]
uVars.ao.ao.measTimeLimit     = 1;

% Frequency
uVars.ao.ao.envDC             = 100e3;
uVars.ao.ao.envUS             = 78e3;

% Geometry & Length
uVars.ao.ao.c                 = 1500; %[m/s]
uVars.ao.ao.distFromPhantom   = 0;    %[m]

% General Operation
uVars.ao.ao.useFrame        = true;
uVars.ao.ao.useGPU          = true;
uVars.ao.ao.useHadamard     = true;
uVars.ao.ao.analyzeSingleCh = false;
uVars.ao.ao.acCoupling      = true;
uVars.ao.ao.useCalibration  = false;
uVars.ao.ao.keepSplits      = false;
uVars.ao.ao.externalIOCtrl  = false;
uVars.ao.ao.displayBuildUp  = true;
uVars.ao.ao.displayBUEvery  = 4;
uVars.ao.ao.dispTimeTable   = true; % relevant only with GUI
uVars.ao.ao.timeToSampleCal = 1;

uVars.ao.ao.cutArtfct    = true;
uVars.ao.ao.artfctIdxVec = 1:15;

uVars.ao.ao.calcMuEff  = false;
uVars.ao.ao.muEffModel = 'Uniform';
uVars.ao.ao.muEffIdxs  = 90:110;

uVars.ao.ao.extSNRDef   = true;
uVars.ao.ao.extPeakIdx  = 126;
uVars.ao.ao.extNoiseIdx = 30:40;

uVars.ao.ao.uploadToTelegram = false;
uVars.ao.ao.telegramChatID   = '-512325870';
uVars.ao.ao.tgBUEvery        = 16;

% Virtual Data
uVars.ao.useVirtualData       = false;
uVars.ao.virtualDataNoiseSTD  = [];

%Live AOI
uVars.ao.ao.limitByN       = true;
uVars.ao.ao.N              = 10;

% IO
uVars.ao.ao.IOPort = 1;
uVars.ao.ao.IOLine = 4;

% Export
uVars.ao.ao.exportData.meas        = false;
uVars.ao.ao.exportData.signal      = false;
uVars.ao.ao.exportData.deMul       = false;
uVars.ao.ao.exportData.reshaped    = false;
uVars.ao.ao.exportData.fft         = false;
uVars.ao.ao.exportData.usCompCmplx = false;

%---------------------%
%Figures Vars
%---------------------%
%Scan 2D Graphics
uVars.figs.reopenFigures       = true;
uVars.figs.sepChIdx            = 1;
uVars.figs.normDispPlaneColors = false;

% Axes Types
uVars.figs.scan1AxType = 'Normal';
uVars.figs.scan2AxType = 'Normal';
uVars.figs.scan3AxType = 'Normal';
uVars.figs.depthAxType = 'Index';

% Navigators
uVars.figs.varsNavLine.labels = {'X', 'Y', 'Z', 'D'};
uVars.figs.varsNavLine.idx    = [1, 1, 126];

uVars.figs.varsNavPlane.labels = {'X', 'Y', 'Z', 'D'};
uVars.figs.varsNavPlane.idx     = [1, 126];

uVars.figs.intExt = 'int';

uVars.figs.validStruct.scanPlane    = true;
uVars.figs.validStruct.scanPlaneLog = true;
uVars.figs.validStruct.navLine      = true;
uVars.figs.validStruct.navLineLog   = true;
uVars.figs.validStruct.navPlane     = true;
uVars.figs.validStruct.navPlaneLog  = true;

% Acousto Optics
uVars.ao.figs.ch        = 1;
uVars.ao.figs.depthIdx  = 146;
uVars.ao.figs.frame     = 1;

uVars.ao.figs.reopenFigures  = false;
uVars.ao.figs.displayFullFFT = true;
uVars.ao.figs.FFTenv         = 250e3; %Hz

uVars.ao.figs.singleChIdx = 1;
uVars.ao.figs.depthAxType = 'Normal';
uVars.ao.figs.intExt      = 'int';

uVars.ao.figs.validStruct.extClk   = false;
uVars.ao.figs.validStruct.usSignal = false;

uVars.ao.figs.validStruct.cropped  = false;

uVars.ao.figs.validStruct.meas     = false;
uVars.ao.figs.validStruct.signal   = false;
uVars.ao.figs.validStruct.deMul    = false;
uVars.ao.figs.validStruct.reshaped = false;

uVars.ao.figs.validStruct.rawFFT         = false;
uVars.ao.figs.validStruct.calibration    = false;
uVars.ao.figs.validStruct.fittingModel   = false;
uVars.ao.figs.validStruct.fittedPowerFFT = false;
uVars.ao.figs.validStruct.finalFFT       = true;

uVars.ao.figs.validStruct.phi    = true;
uVars.ao.figs.validStruct.rawPhi = false;
uVars.ao.figs.validStruct.phiLog = false;

%----------------------%
% File System Variables
%----------------------%
% Scan 2D FileSystem Variables
uVars.fileSystem.saveResults  = false;
uVars.fileSystem.saveFigs     = false; %not implemented yet
uVars.fileSystem.saveAO       = false;

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
close all;

s2D.setUserVars(uVars);
s2D.configureScan();
res = s2D.scan();
