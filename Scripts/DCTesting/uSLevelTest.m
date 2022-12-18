close all
clear all
clc

ao = acoustoOptics();
ao.init();

%%
ao.setMeasLimit(1)

uVars = ao.createUserVars();

%-------------------------
% Algo and Peripheral Vars
%-------------------------
% US Signal
uVars.ao.cycPerPulse       = 1;      %[#] 
uVars.ao.fSin              = 1.25e6; %[Hz]              
uVars.ao.fSqnc             = 5e3; %[Hz]
uVars.ao.frameTime         = 0.002; %[s]

% Sampling Clk
uVars.ao.fs                = 20e6;    %[Hz]
uVars.ao.sClkDcyc          = 50;     %[%]
uVars.ao.fgClk             = 100e6;  %[S/s]

% Digitizer
uVars.ao.timeToSample   = 256; %[s]
uVars.ao.channels       = 16; %[#]
uVars.ao.extCropDef     = false;
uVars.ao.extCropSamples = 0;

% Frequency
uVars.ao.envDC             = 100e3;
uVars.ao.envUS             = 78e3;

% Geometry & Length
uVars.ao.c                 = 1500;   %[m/s]
uVars.ao.distFromPhantom   = 0; %[m]

% General Operation
uVars.ao.useFrame            = true;
uVars.ao.useGPU              = true;
uVars.ao.useHadamard         = true;
uVars.ao.contHadamard        = false;
uVars.ao.highResAO           = false;
uVars.ao.analyzeSingleCh     = false;
uVars.ao.contSpeckleAnalysis = false;
uVars.ao.useCalibration      = false;
uVars.ao.acCoupling          = true;
uVars.ao.autoCalibration     = false;

uVars.ao.timeToSampleCal = 0;

uVars.ao.cutArtfct           = true;
uVars.ao.artfctIdxVec        = 1:20;

uVars.ao.calcMuEff = false;
uVars.ao.muEffModel = 'Uniform';
uVars.ao.muEffIdxz = [];

uVars.ao.useVirtualData      = false;
uVars.ao.virtualDataNoiseSTD = 0;
uVars.ao.limitByN            = true;
uVars.ao.N                   = 10;
uVars.ao.dispTimeTable       = true; % relevant only with GUI

uVars.ao.skipParamsCheck = false;

% fGen
uVars.ao.usPower = 100; 

% IO
uVars.ao.IOPort = 1;
uVars.ao.IOLine = 4;

% Export
uVars.ao.exportData.meas        = false;
uVars.ao.exportData.signal      = false;
uVars.ao.exportData.deMul       = false;
uVars.ao.exportData.reshaped    = false;
uVars.ao.exportData.fft         = false;
uVars.ao.exportData.usCompCmplx = false;

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

uVars.figs.validStruct.cropped  = false;

uVars.figs.validStruct.meas     = false;
uVars.figs.validStruct.signal   = false;
uVars.figs.validStruct.deMul    = false;
uVars.figs.validStruct.reshaped = false;

uVars.figs.validStruct.rawFFT         = false;
uVars.figs.validStruct.calibration    = false;
uVars.figs.validStruct.fittingModel   = false;
uVars.figs.validStruct.fittedPowerFFT = false;
uVars.figs.validStruct.finalFFT       = true;

uVars.figs.validStruct.phi    = true;
uVars.figs.validStruct.rawPhi = true;
uVars.figs.validStruct.phiLog = true;

%----------------------%
% File System Variables
%----------------------%
uVars.fileSystem.saveMeas       = false;
uVars.fileSystem.saveSignal     = false;
uVars.fileSystem.saveDeMul      = false;
uVars.fileSystem.saveReshaped   = false;
uVars.fileSystem.saveFFT        = false;
uVars.fileSystem.savePhiChCmplx = false;

uVars.fileSystem.saveResults = true;
uVars.fileSystem.saveFigs    = false; %not implemented yet

uVars.fileSystem.dirPath     = "D:\USLevel";
uVars.fileSystem.projName    = "USLvl-0";
uVars.fileSystem.resDirName  = "Results";

uVars.fileSystem.extProject     = false;
uVars.fileSystem.useExtVarsPath = false;
uVars.fileSystem.extVarsPath    = [];


%%
idxVec = 0:1:10;
usIntVec = 100*(0.75.^idxVec);

for i=1:length(idxVec)
    uVars.ao.usPower = usIntVec(i);
    uVars.fileSystem.projName = sprintf("USLvl-%d", idxVec(i));
    ao.setVars(uVars);
    ao.configPeripherals();
    aoVars = ao.getVars();
    res    = ao.runAcoustoOptics();
end





