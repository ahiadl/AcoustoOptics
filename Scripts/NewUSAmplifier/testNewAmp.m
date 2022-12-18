close all;
clear all;
clc;
instrreset;

%%
addpath(genpath(pwd));
fprintf("Connecting to function generator\n");
fg = fGen();
fg.connect();
fprintf("Connecting to IO\n");
daqIO = IO();
daqIO.connect();

IOVars = daqIO.uVarsCreate;

IOVars.mode = 0; % Output Only
IOVars.port = 1;
IOVars.line = 4;

daqIO.allocPorts(IOVars);
daqIO.open();
daqIO.close();

%% Configure the function generator
daqIO.close();

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
samplesInSinCycle = sClk / fsin;
tVecSin  = (0:1:((samplesInSinCycle)-1)) / sClk;

%--------------------------------------
%     %------- Square-Pulse
%     Npulse = (sClk/fsin)/2;
%     sig = ones(1,Npulse);
% 
%     delaySamples = floor(delay * sClk);
%     padLen = N-length(sig)-delaySamples;
%     sig = [zeros(1,delaySamples), sig, zeros(1,padLen)];
%--------------------------------------
%------- Sin-Pulse
% sigAmp = 0.01;
% Nsin = length(tVecSin);
% sig = sigAmp*sin(2*pi*fsin*tVecSin);
% delaySamples = floor(delay * sClk);
% padLen = N-length(sig)-delaySamples;
% sig = [zeros(1,delaySamples), sig, zeros(1,padLen)];

%--------------------------------------
%------- Delta-Pulse
%     pulse = [0, 0.25, 1,  0.25, 0];
%     padLen = N-length(pulse);
%     sig = [pulse, zeros(1, padLen)]; 

%     sig = convSigNorm;
%     sig = convSigFull;
%--------------------------------------

%--------------------------------------
%------- Hadamard Sequence
sigAmp = 0.11;
hadOrder = 251;
sequence = createSMatrix(hadOrder);
sig = sin(2*pi*fsin*tVecSin);
sigMat = repmat(sig', 1, hadOrder);
sigMat = sigMat .* sequence(1,:);
sigSeq = sigAmp*sigMat(:)';

    NSeq = samplesInSinCycle * hadOrder;
%     padLen = 2*Nseq;
    padLen = 0;
    sig = [sigSeq, zeros(1, padLen)]; 

%     sig = convSigNorm;
%     sig = convSigFull;
%--------------------------------------

fVecSig = (sClk/N) *  ( (-N/2) : 1 : (N/2)-1 ); 
tVecSig = (0:1:N-1)./sClk;

%--------------------------------------
% Config the Function Generator
dataCh1.data    = sig;
dataCh1.dataLen = length(dataCh1.data);

fg.setProperties(fGenVars.ch{1}, fGenVars.ch{2});
fg.setData(dataCh1, dataCh2);

fg.configChannel(1); % this will enable the channel

fg.disableOutput(1);
updateFGen = true;

daqIO.open();