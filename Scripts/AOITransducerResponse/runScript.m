close all;
clear all;
clc;
instrreset;
%% Init Objects
addpath(genpath(pwd));
fprintf("Connecting To stages\n");
stages = stages('Zaber', 'COM3');
stages.connect(); 
stages = [];
fprintf("Connecting to digitizer\n");
daq    = Digitizer();
daq.connect();
fprintf("Creating scan object\n");
cs     = contScan(daq, stages);
fprintf("Connecting to function generator\n");
fg = fGen();
fg.connect();
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
daqIO.close();
if updateFGen
    fg.reset();
    fGenVars = fg.uVarsCreate();

    sClk = 100e6;
    fClk = 100e6;

    fGenVars.ch{1}.daqFilter     = 50;
    fGenVars.ch{1}.amp           = 2;
    fGenVars.ch{1}.bias          = 0;
    fGenVars.ch{1}.triggerOwner  = true;
    fGenVars.ch{1}.triggerWidth  = 512;
    fGenVars.ch{1}.useExtSclkSrc = false;

    fGenVars.ch{2}.daqFilter     = 50;
    fGenVars.ch{2}.amp           = 1;
    fGenVars.ch{2}.bias          = 0;
    fGenVars.ch{2}.triggerOwner  = false;
    fGenVars.ch{2}.triggerWidth  = 16;
    fGenVars.ch{2}.useExtSclkSrc = false;

    fGenVars.ch{1}.Sclk = sClk;
    fGenVars.ch{2}.Sclk = sClk;

    dataCh1 = fg.getChDataStruct();
    dataCh2 = fg.getChDataStruct();


    %--------------------------------------
    %------- Create the Signal 
    N = 65536;
    sigTime = N/sClk;
    delay = 0;
    fsin = 1.25e6;
    
    %--------------------------------------
    %------- positive-negative Square-Pulse
%     Npulse = (sClk/fsin)/4;
%     sig1 = ones(1,Npulse);
%     sig2 = zeros(1,Npulse);
%     sig = [sig1, -sig1];
%     
%     delaySamples = floor(delay * sClk);
%     padLen = N-length(sig)-delaySamples;
%     sig = [zeros(1,delaySamples), sig, zeros(1,padLen)];

%     %------- Square-Pulse
%     Npulse = (sClk/fsin)/2;
%     sig = ones(1,Npulse);
% 
%     delaySamples = floor(delay * sClk);
%     padLen = N-length(sig)-delaySamples;
%     sig = [zeros(1,delaySamples), sig, zeros(1,padLen)];
    
%     %------- Sin-Pulse
    tVecSin  = (0:1:((sClk/fsin)-1))*(1/sClk);
    Nsin = length(tVecSin);

    sig = sin(2*pi*fsin*tVecSin);

    delaySamples = floor(delay * sClk);
    padLen = N-length(sig)-delaySamples;
    sig = [zeros(1,delaySamples), sig, zeros(1,padLen)];
      
      %--------------------------------------
%     %------- Sin (Pos-Neg)-Pulse
%     tVecSin  = (0:1:((sClk/fsin)-1))*(1/sClk);
%     Nsin = length(tVecSin);
% 
%     sig1 = sin(2*pi*fsin*tVecSin);
% %     sig = [sig1(1:20), -sig1(1:20)];
% %     sig = [sig1(1:40)];
%     sig = [sig1, sig1];
%     delaySamples = floor(delay * sClk);
%     padLen = N-length(sig)-delaySamples;
%     sig = [zeros(1,delaySamples), sig, zeros(1,padLen)];

    %--------------------------------------
    %------- Delta-Pulse
%     pulse = [0, 0.25, 1,  0.25, 0];
%     padLen = N-length(pulse);
%     sig = [pulse, zeros(1, padLen)]; 

%     sig = convSigNorm;
%     sig = convSigFull;
    %--------------------------------------
    
    fVecSig = (sClk/N) *  ( (-N/2) : 1 : (N/2)-1 ); 
    tVecSig = (0:1:N-1)./sClk;
    
    figure();
    subplot(1,2,1)
    plot(tVecSig*1e6, sig);
    xlabel("t[\mu s]");
    subplot(1,2,2)
    plot(fVecSig/1e6, abs(fftshift(fft(fftshift(sig)))).^2)
    xlabel("f[MHz]");
    
    %--------------------------------------
    % Config the Function Generator
    dataCh1.data    = sig;
    dataCh1.dataLen = length(dataCh1.data);
    
    fg.setProperties(fGenVars.ch{1}, fGenVars.ch{2});
    fg.setData(dataCh1, dataCh2);

    fg.configChannel(1); % this will enable the channel

    samplesInClk = floor(sClk/fClk);
    fClkReal     = sClk/samplesInClk;
    if mod(samplesInClk,2)
        samplesInDCHigh  = floor(0.5*samplesInClk)+1;
        samplesInDCLow   = samplesInClk - samplesInDCHigh;
    else
        samplesInDCHigh = 0.5*samplesInClk;
        samplesInDCLow = samplesInDCHigh;
    end
    fprintf("Requested sampling rate %.2f [MHz]\n", fClk/1e6);
    fprintf("Configured sampling rate %.2f [MHz]\n", fClkReal/1e6);

    dataCh2.data = repmat([ones(1,samplesInDCHigh), zeros(1,samplesInDCLow)], 1, 16); % 5Mhz clock
    dataCh2.dataLen = length(dataCh2.data);

    fg.setData(dataCh1, dataCh2);
%     fg.configChannel(2); % this will enable the channel

    updateFGen = true;
end
daqIO.open();
%% Stages Config 
% stages.setVelocity('X', 10);
% stages.setVelocity('Y', 10);
% stages.setVelocity('Z', 10);

%% Genreal Parameters
genVars.numOfCh = 1;
genVars.fs = 100e6;

timeToSample = sigTime;
numOfSamples = round(timeToSample * genVars.fs);
genVars.numOfSamples = numOfSamples;

genVars.delay = delay; 
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
uVarsCS.stabTime      = 1;

uVarsCS.binSize      = 1024;
uVarsCS.binType      = 'mid'; %'mid', 'post';
uVarsCS.adaptBin     = true;

uVarsCS.axScan  = 'Z';  %'X', 'Y', 'Z'
uVarsCS.axDisc1 = 'Y';  %'X', 'Y', 'Z'
uVarsCS.axDisc2 = 'X';  %'X', 'Y', 'Z'
uVarsCS.spanType = 'center'; %'center', 'limits'

uVarsCS.axScanSpan  = 19.3/4; %[mm]
uVarsCS.axDisc1Span = 0;  %[mm]
uVarsCS.axDisc2Span = 0;  %[mm]

uVarsCS.axScanStride  = 0.1; %[mm]
uVarsCS.axDisc1Stride = 0;   %[mm]
uVarsCS.axDisc2Stride = 0;   %[mm]

uVarsCS.axScanRef  = 19.3;    %[mm] %94.6
uVarsCS.axDisc1Ref = 47.4; %[mm]
uVarsCS.axDisc2Ref = 91; %[mm]

% Post Processing
% Internal:
uVarsCS.doPP        = true;
% uVarsCS.filter      = @userFilter;
uVarsCS.filter      = 0;
% uVarsCS.filter      = 'LPF'; %'LPF', 'BPF', 'HPF'

% External (User functions)
% uVarsCS.postProc    = @postProc;
% uVarsCS.initPlots   = @userInitPlots;
% uVarsCS.plotResults = @userPlotResults;

% uVarsCS.filter      = 0;
% uVarsCS.postProc    = 0;
% uVarsCS.initPlots   = 0;
% uVarsCS.plotResults = 0;

% userAux.sensitivity    = 837;
userAux.idx = 128;
N = numOfSamples;
fs = fClkReal;
userAux.fVec = (fs/N) *  ( (-N/2) : 1 : (N/2)-1 );
uVarsCS.userAux  = userAux;

uVarsCS.plotMode  = 'mid'; % 'max', 'mid', 'user'
uVarsCS.chToPlot  = 1;
uVarsCS.posToPlot = 0;

uVarsCS.saveData = true;
uVarsCS.dir = ".\";
uVarsCS.filename = "PA-Flat-ZAxis-FullScan-1MHz";

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
digiVars.extClk  = false;
digiVars.extTrig = true;

digiVars.timeToSample = sigTime;
digiVars.avgNum       = 1;
digiVars.numMeas      = 1;

digiVars.draw = false;

daq.setVars(digiVars);
daq.configure();

data = daq.acquire();
figure(); 
plot(data); 
yyaxis right
plot(data1mhz)

%% Scan And Measure 
close all
daqIO.open();
resCs = cs.scan();
daqIO.close();