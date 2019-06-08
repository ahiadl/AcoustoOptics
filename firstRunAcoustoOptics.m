close all;
clear all;
clc;

%should affect fGen algo & digitizer

uVars = acoustoOptics.uVarsCreate();

uVars.fSin              = 1.25e6;              
uVars.fTrain            = 19e3;
uVars.cycInPulse        = 4; 
uVars.channels          = 4; %update in digitizer
uVars.phantomDepth      = 4.2e-2;
uVars.distFromPhantom   = 0;
uVars.fExtClk           = 5e6; %fs
uVars.fSclk             = 100e6;     %update in fGen
uVars.timeToSample      = 0.25;
uVars.extClkDcyc        = 50; % [%]
uVars.IOPort            = 1;
uVars.IOLine            = 4;
uVars.useGPU            = true; %algo, digitizer;
uVars.exportRawData     = false;


% uVars.graphicRequest = acoustoOptics.createGraphicsRequest();
% uVars.graphicRequest.vars.internal = true;
% uVars.graphicRequest.vars.extClk = true;
% 
% uVars.graphics.internal = true;
% uVars.graphics.ch       = 1;
% uVars.graphics.pos      = 1;
% uVars.graphics.mask     = ones(1,9);

% uVars.general.useGPU
% uVars.general.exportRawData = false;


acoustoOptics = acoustoOptics();

acoustoOptics.init();

% acoustoOptics.setMeasVars(uVars);
% 
% % acoustoOptics.prepareGraphics();
% 
% acoustoOptics.configPeripherals();
% 
% tic
% acoustoOptics.measureAndAnlayse();
% toc

acoustoOptics.fullMeasureAndAnalyse(uVars);

timeTable = acoustoOptics.getTimeTable();
% clear acoustoOptics



