close all;
clear all;
clc;
instrreset;

%% Init Objects
addpath(genpath(pwd));
fprintf("Connecting To Stages\n");
stages = stages('Zaber', 'COM3');
stages.connect(); 
fprintf("Connecting to digitizer\n");
daq    = Digitizer();
daq.connect();
fprintf("Creating scan object\n");
cs     = contScan(daq, stages);

%%
aofg = AOFGen();
aofg.init();

uVars = AOFGen.createUserVars();

uVars.fgClk   = 100e6;

uVars.fSig        = 1.25e6;
uVars.fSqnc       = 20e3;
uVars.delay       = 0;
uVars.sigType     = 'Hadamard';
uVars.cycPerPulse = 1;
uVars.sigPower    = 20; %
uVars.fsClk       = 20e6;

aofg.setVars(uVars);
vars = aofg.getVars();
aofg.config();

%% Genreal Parameters
genVars.numOfCh = 1;
genVars.fs = 1e6;

timeToSample = 0.010;
numOfSamples = round(timeToSample * genVars.fs);
genVars.numOfSamples = numOfSamples;

genVars.delay = 0; 
genVars.tVec = genVars.delay  + (0:1:genVars.numOfSamples-1)/genVars.fs;

%% CS Config
uVarsCS = cs.createUserVars();

uVarsCS.numOfCh       = genVars.numOfCh;
uVarsCS.tVec          = genVars.tVec;
uVarsCS.laserRepRate  = 5e3; %[Hz]

uVarsCS.scanDirection = 'bi';   %'bi', 'uni'
uVarsCS.scanType      = 'disc'; % 'cont', 'disc', 'pos';
uVarsCS.stabDist      = 0; %[mm]
uVarsCS.infSinglePos  = false;
uVarsCS.limSinglePos  = 1024;
uVarsCS.allAtOnce     = true;
uVarsCS.stabTime      = 0;

uVarsCS.binSize      = 1;
uVarsCS.binType      = 'mid'; %'mid', 'post';
uVarsCS.adaptBin     = true;

uVarsCS.axScan  = 'Y';  %'X', 'Y', 'Z'
uVarsCS.axDisc1 = 'Z';  %'X', 'Y', 'Z'
uVarsCS.axDisc2 = 'X';  %'X', 'Y', 'Z'
uVarsCS.spanType = 'center'; %'center', 'limits'

uVarsCS.axScanSpan  = 30; %[mm]
uVarsCS.axDisc1Span = 0;  %[mm]
uVarsCS.axDisc2Span = 0;  %[mm]

uVarsCS.axScanStride  = 1;  %[mm]
uVarsCS.axDisc1Stride = 0;  %[mm]
uVarsCS.axDisc2Stride = 0;   %[mm]

uVarsCS.axScanRef  = 69.3;     %[mm] 100
uVarsCS.axDisc1Ref = 108;   %[mm] 71.5
uVarsCS.axDisc2Ref = 0;     %[mm]

% Post Processing
% Internal:
uVarsCS.doPP    = false;
% uVarsCS.filter = @userFilter;
uVarsCS.filter  = 0;
% uVarsCS.filter = 'LPF'; %'LPF', 'BPF', 'HPF'

% External (User functions)
% uVarsCS.postProc    = @MTpostProc;
% uVarsCS.initPlots   = @MTuserInitPlots;
% uVarsCS.plotResults = @MTuserPlotResults;

uVarsCS.userAnalysis = MTUA();

uVarsCS.plotMode  = 'mid'; % 'max', 'mid', 'user'
uVarsCS.chToPlot  = 1;
uVarsCS.posToPlot = 0;

uVarsCS.saveData = false;
uVarsCS.dir = "./Measurements/OpticsPatterns";
uVarsCS.filename = "Transmission-Point";

%Telegram
uVarsCS.tg.text   = false;
uVarsCS.tg.figure = false;
uVarsCS.tg.chatID = '-502071516';
uVarsCS.tg.rep    = 1;

cs.setVars(uVarsCS);
csVars = cs.getVars();

%% DAQ Config
digiVars = Digitizer.uVarsCreate();

digiVars.mode      = 'NPT'; 
digiVars.fs        = genVars.fs;
digiVars.useGPU    = false;
digiVars.channels  = genVars.numOfCh; 

digiVars.triggerDelay = 0;
digiVars.extClk  = false;
digiVars.extTrig = true;

digiVars.timeToSample = timeToSample;
digiVars.avgNum       = 1;
digiVars.numMeas      = csVars.daqAcqNum;

digiVars.coupling = 0;

digiVars.draw = false;

daq.setVars(digiVars);
daq.configure();

%% Scan And Measure 
close all
resCs = cs.scan();

% %%
% size(resCs)
% res = squeeze(max(resCs, [], 5));
% 
% % resAlign = res;
% % resAlign(:, 1:2:end) = circshift(resAlign(:,1:2:end), -2, 1);
% 
% figure();
% imagesc(csVars.scanVecBin, csVars.disc1Vec, res); colorbar
% axis tight equal
