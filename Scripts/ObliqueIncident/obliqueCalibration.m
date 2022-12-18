% close all;
% clear all;
% clc;
% instrreset;

%% Init Objects
addpath(genpath(pwd));
fprintf("Connecting To stages\n");
stages = stages('Zaber', 'COM3');
stages.connect(); 
fprintf("Connecting to digitizer\n");
daq    = Digitizer();
daq.connect();
fprintf("Creating scan object\n");
cs     = contScan(daq, stages);

%% Genreal Parameters
genVars.numOfCh = 16;
genVars.fs = 1e6;

timeToSample = 0.001;
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
uVarsCS.scanType      = 'pos'; % 'cont', 'disc', 'pos';
uVarsCS.stabDist      = 0; %[mm]
uVarsCS.infSinglePos  = true;
uVarsCS.limSinglePos  = 1024;
uVarsCS.allAtOnce     = false;
uVarsCS.stabTime      = 0;

uVarsCS.binSize      = 1;
uVarsCS.binType      = 'mid'; %'mid', 'post';
uVarsCS.adaptBin     = true;

uVarsCS.axScan  = 'Y';  %'X', 'Y', 'Z'
uVarsCS.axDisc1 = 'Z';  %'X', 'Y', 'Z'
uVarsCS.axDisc2 = 'X';  %'X', 'Y', 'Z'
uVarsCS.spanType = 'center'; %'center', 'limits'

uVarsCS.axScanSpan  = 1; %[mm]
uVarsCS.axDisc1Span = 1; %[mm]
uVarsCS.axDisc2Span = 0;  %[mm]

uVarsCS.axScanStride  = 1; %[mm]
uVarsCS.axDisc1Stride = 1;   %[mm]
uVarsCS.axDisc2Stride = 0;   %[mm]

uVarsCS.axScanRef  = 70; %[mm] %94.6
uVarsCS.axDisc1Ref = 81.5; %[mm]
uVarsCS.axDisc2Ref = 10;     %[mm]

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
% uVarsCS.filter      = 0;
% uVarsCS.postProc    = 0;
% uVarsCS.initPlots   = 0;
% uVarsCS.plotResults = 0;

% userAux.sensitivity    = 837;
userAux.idx = 128;
N = numOfSamples;
fs = genVars.fs;
userAux.fVec = (fs/N) *  ( (-N/2) : 1 : (N/2)-1 );
uVarsCS.userAux  = userAux;

uVarsCS.plotMode  = 'mid'; % 'max', 'mid', 'user'
uVarsCS.chToPlot  = 1;
uVarsCS.posToPlot = 0;

uVarsCS.saveData = true;
uVarsCS.dir = ".\";
uVarsCS.filename = "Reflection-NoPDMS";

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
digiVars.avgNum       = 64;
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
