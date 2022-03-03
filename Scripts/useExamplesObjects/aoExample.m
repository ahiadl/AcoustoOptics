close all;
clear all;
clc;

ao = acoustoOptics();
ao.init();
% stages = stages('COM3');
% stages.connect();

% stages.assignStagesAxes(['X','Y'], [1,2]);
% stages.moveStageAxisAbs('X', 101);
% stages.moveStageAxisAbs('Y', 88);
%%

ao.setMeasLimit(0.2)

uVars = ao.createUserVars();

%-------------------------
% Algo and Peripheral Vars
%-------------------------
% US Signal
uVars.ao.cycPerPulse       = 1;      %[#] 
uVars.ao.fSin              = 1.25e6; %[Hz]              
uVars.ao.fSqnc             = 10e3; %[Hz]
uVars.ao.frameTime         = 0.002; %[s]

% Sampling Clk
uVars.ao.fs                = 5e6;    %[Hz]
uVars.ao.sClkDcyc          = 50;     %[%]
uVars.ao.fgClk             = 100e6;  %[S/s]

% Digitizer
uVars.ao.timeToSample      = 0.1; %[s]
uVars.ao.channels          = 1; %[#]

% Frequency
uVars.ao.envDC             = 100e3;
uVars.ao.envUS             = 78e3;

% Geometry & Length
uVars.ao.c                 = 1500;   %[m/s]
uVars.ao.distFromPhantom   = 0; %[m]

% General Operation
uVars.ao.useFrame            = true;
uVars.ao.useGPU              = true;
uVars.ao.useHadamard         = false;
uVars.ao.contHadamard        = false;
uVars.ao.highResAO           = false;
uVars.ao.analyzeSingleCh     = true;
uVars.ao.contSpeckleAnalysis = false;

uVars.ao.cutArtfct           = false;
uVars.ao.artfctIdxVec        = [];

uVars.ao.useVirtualData   = false;
uVars.ao.limitByN         = true;
uVars.ao.N                = 10;
uVars.ao.dispTimeTable    = true; % relevant only with GUI

% IO
uVars.ao.IOPort = 1;
uVars.ao.IOLine = 4;

% Export
uVars.ao.exportData.meas        = true;
uVars.ao.exportData.signal      = true;
uVars.ao.exportData.deMul       = true;
uVars.ao.exportData.reshaped    = true;
uVars.ao.exportData.fft         = true;
uVars.ao.exportData.usCompCmplx = true;

%----------------------------
% Figures Variables
%---------------------------
uVars.figs.ch    = 1;
uVars.figs.depthIdx  = 30;
uVars.figs.frame = 1;

uVars.figs.reopenFigures = false;
uVars.figs.displayFullFFT = true;
uVars.figs.FFTenv = 250e3; %Hz

uVars.figs.singleChIdx = 1;
uVars.figs.depthAxType = 'Index';
uVars.figs.intExt = 'int';

uVars.figs.validStruct.extClk   = false;
uVars.figs.validStruct.usSignal = false;

uVars.figs.validStruct.meas     = false;
uVars.figs.validStruct.signal   = false;
uVars.figs.validStruct.deMul    = false;
uVars.figs.validStruct.reshaped = false;

uVars.figs.validStruct.unFittedFFTAllCh    = true;
uVars.figs.validStruct.unFittedFFTSingleCh = true;
uVars.figs.validStruct.fittedFFTSingleCh   = true;
uVars.figs.validStruct.fittedFFTAllCh      = true;
uVars.figs.validStruct.avgFFT              = true;
uVars.figs.validStruct.normFFT             = true;

uVars.figs.validStruct.phi    = true;
uVars.figs.validStruct.rawPhi = true;
%----------------------%
% File System Variables
%----------------------%
uVars.fileSystem.saveMeas       = false;
uVars.fileSystem.saveSignal     = false;
uVars.fileSystem.saveDeMul      = false;
uVars.fileSystem.saveReshaped   = false;
uVars.fileSystem.saveFFT        = false;
uVars.fileSystem.savePhiChCmplx = false;

uVars.fileSystem.saveResults = false;
uVars.fileSystem.saveFigs    = false; %not implemented yet

uVars.fileSystem.dirPath     = ".\Scripts\objectUseExamples\examplesData";
uVars.fileSystem.projName    = "Example AO";
uVars.fileSystem.resDirName  = "Results";

uVars.fileSystem.extProject     = false;
uVars.fileSystem.useExtVarsPath = false;
uVars.fileSystem.extVarsPath    = [];

ao.setVars(uVars);
ao.configPeripherals();
aoVars = ao.getVars();
res    = ao.runAcoustoOptics();
%% Single Point AO
% %% Live Acousto Optics
% uVars.ao.limitByN            = true;
% uVars.ao.N                   = 3;
% uVars.fileSystem.saveRawData = false;
% uVars.fileSystem.projName = "Example Live AO";
% 
% ao.setMeasVars(uVars);
% ao.configPeripherals();
% aoVars = ao.getVars();
% res    = ao.liveAcoustoOptics();
% 
% %% Loaded Data
% uVars.fileSystem.saveResults = false;
% uVars.fileSystem.saveFigs    = false;
% 
% rawDataPath = ".\Scripts\objectUseExamples\examplesData\07-May-2020 19-47-42-Example AO\AO-Results.mat";
% varsPath    = ".\Scripts\objectUseExamples\examplesData\07-May-2020 19-47-42-Example AO\AO-Vars.mat";
% 
% rawData     = ao.loadRawDataToAO(rawDataPath);
% vars        = load(varsPath);
% 
% %Important so no additional saving will occur
% vars.uVars.fileSystem.saveRawData        = false;
% vars.uVars.fileSystem.saveNetSignal      = false;
% vars.uVars.fileSystem.saveDemultiplexed  = false;
% vars.uVars.fileSystem.saveReshapedSignal = false;
% vars.uVars.fileSystem.saveFFT            = false;
% vars.uVars.fileSystem.savePhiChCmplx     = false;
% 
% vars.uVars.fileSystem.saveResults = false;
% vars.uVars.fileSystem.saveFigs    = false; 
% 
% ao.setMeasVars(vars.uVars);
% loadedDataRes  = ao.analyseLoadedData();



