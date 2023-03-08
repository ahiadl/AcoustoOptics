close all;
clear all;
clc;
instrreset;
%% Init Objects
% addpath(genpath(pwd));
fprintf("Connecting To stages\n");
stages = stages('Zaber', 'COM3');
stages.connect(); 
stages.getPosition();
fprintf("Connecting to digitizer\n");
daq    = Digitizer();
daq.connect();
fprintf("Creating scan object\n");
cs     = contScan(daq, stages);

%% Dummy Signal to Measure
fprintf("Connecting to function generator\n");
fg = AOFGen();
fg.init();

fprintf("Connecting to IO\n");
daqIO = IO();
daqIO.connect();

%% Configure the IO
IOVars = daqIO.uVarsCreate;

IOVars.mode = 0; % Output Only
IOVars.port = 1;
IOVars.line = 4;

daqIO.allocPorts(IOVars);
daqIO.open();
daqIO.close();
%% Configure the function generator
updateFGen = true;
% daqIO.close();
if updateFGen
    uVarsFGen = AOFGen.createUserVars();

    uVarsFGen.fgClk   = 100e6;
    
    uVarsFGen.fSig    = 1.25e6;
    uVarsFGen.fSqnc   = 7.5e3;
    uVarsFGen.delay   = 0;
    uVarsFGen.sigType = 'Hadamard'; %'SinePulse'
    uVarsFGen.cycPerPulse = 1;
    
    uVarsFGen.fsClk   = 20e6;
    
    uVarsFGen.clkPower = 50;  %[%]
    uVarsFGen.sigPower = 100;  %[%]
    
    uVarsFGen.shift = 0;

    fg.setVars(uVarsFGen);
    vars = fg.getVars();
    fg.config();

    clkData = fg.vars.clkData;
    sigData = fg.vars.sigData;
end
% daqIO.open();
% daqIO.close();
fs    = vars.sClk.fs;  
fSqnc = vars.sig.fSqnc;

sigTime  = size(sigData,1) / uVarsFGen.fgClk;
sigDepth = sigTime*1500*1e3; %[mm]

figure();
subplot(1,2,1)
plot(sigData);
subplot(1,2,2)
plot(clkData)

delay = 0;
%% Genreal Parameters
genVars.numOfCh = 1;
% genVars.fs = 100e6;
genVars.fs = fs;

timeToSample = sigTime - (144+16)/fs; % this is to make sure 64 clk cycles are available for re-arming the daq.
% timeToSample = sigTime;

numOfSamples = round(timeToSample * genVars.fs);
genVars.numOfSamples = numOfSamples;

genVars.delay = delay; 
genVars.tVec = genVars.delay  + (0:1:genVars.numOfSamples-1)/genVars.fs;

%% CS Config
uVarsCS = cs.createUserVars();

uVarsCS.numOfCh       = genVars.numOfCh;
uVarsCS.tVec          = genVars.tVec;
uVarsCS.laserRepRate  = fSqnc; %[Hz]

uVarsCS.scanDirection = 'bi';   %'bi', 'uni'
uVarsCS.scanType      = 'pos'; % 'cont', 'disc', 'pos';
uVarsCS.stabDist      = 0; %[mm]
uVarsCS.infSinglePos  = false;
uVarsCS.limSinglePos  = 1024;
uVarsCS.allAtOnce     = true;
uVarsCS.stabTime      = 1;

uVarsCS.binSize      = 2048;
uVarsCS.binType      = 'mid'; %'mid', 'post';
uVarsCS.adaptBin     = true;

uVarsCS.axScan  = 'Y';  %'X', 'Y', 'Z'
uVarsCS.axDisc1 = 'Z';  %'X', 'Y', 'Z'
uVarsCS.axDisc2 = 'X';  %'X', 'Y', 'Z'
uVarsCS.spanType = 'center'; %'center', 'limits'

uVarsCS.axScanSpan  = 1; %[mm]
uVarsCS.axDisc1Span = 0;  %[mm]
uVarsCS.axDisc2Span = 0;  %[mm]

uVarsCS.axScanStride  = 1; %[mm]
uVarsCS.axDisc1Stride = 0;   %[mm]
uVarsCS.axDisc2Stride = 0;   %[mm]

uVarsCS.axScanRef  = 122.4;  %[mm]
uVarsCS.axDisc1Ref = 131.1;  %[mm]
uVarsCS.axDisc2Ref = 85.1;   %[mm]

% Post Processing
% Internal:
uVarsCS.doPP        = true;
uVarsCS.filter      = 0;

uVarsCS.userAnalysis  = AOTransducerUA();

uVarsCS.plotMode  = 'mid'; % 'max', 'mid', 'user'
uVarsCS.chToPlot  = 1;
uVarsCS.posToPlot = 0;

uVarsCS.saveData = true;
uVarsCS.dir      = "C:\Users\sahiadl.EED\OneDrive - Technion\Graduate\AcoustoOpticSystem\Measurements\Transducer Pressure Field\Calibrated\";
uVarsCS.filename = "UnFocused-1.25MHz-FocalPoint";

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

digiVars.triggerDelay = delay;
digiVars.extClk  = true;
digiVars.extTrig = true;

digiVars.timeToSample = timeToSample;
digiVars.avgNum       = 1;
digiVars.numMeas      = uVarsCS.binSize;

digiVars.draw = false;

daq.setVars(digiVars);
daq.configure();

data = daq.acquire();
% figure();
% subplot(1,2,1)
% plot(mean(data,3)); 
% subplot(1,2,2)
% imagesc(squeeze(data))

%% Scan And Measure 
close all
clc
daqIO.open();
resCS = cs.scan();
daqIO.close();


%% 
p2p = squeeze(peak2peak(resCS, 5));
p2pLpf = conv(p2p, (1/5)*[1,1,1,1,1], 'same');
env = envelope(p2pLpf, 10, 'peak');
figure(); 
plot(csVars.scanVecBin, p2p); hold on;
plot(csVars.scanVecBin, p2pLpf);
plot(csVars.scanVecBin, env)

%%
dir = "C:/Users/sahiadl.EED/OneDrive - Technion/Graduate/AcoustoOpticSystem\Measurements/Transducer Pressure Field/Calibrated";
file = "15-Jan-2023 11-48-46-Focused-1.25MHz-FocalPoint-HadamardSqnc";

resFocal = load(sprintf("%s/%s.mat", dir, file));
