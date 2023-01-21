close all
clear all
% clc
instrreset

fg = fGen();
fg.connect();
%%

fgClk = 100e6;

fGenVars = fg.uVarsCreate();

fGenVars.ch{1}.daqFilter     = 50;
fGenVars.ch{1}.amp           = 1;
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

fGenVars.ch{1}.Sclk = fgClk;
fGenVars.ch{2}.Sclk = fgClk;

dataCh1 = fg.getChDataStruct();
dataCh2 = fg.getChDataStruct();



%% Create Signal
% NOTE: number of samples in the data for each channel should be a
% multiplication of 16.

N = 20080;
% AO Signal
fsin = 1.25e6;
% fSqnc = 5e3;
delay = 0;

sigTime = N/fgClk;
NSin = fgClk/fsin;

% Create a sin
tVecSin  = (1:1:NSin)*(1/fgClk);
Nsin = length(tVecSin);

tVecSinPad = (0:1:N-1)./fgClk;
fVecSig = (fgClk/N) *  ( (-N/2) : 1 : (N/2)-1 ); 
sig = sin(2*pi*fsin*tVecSin);

delaySamples = floor(delay * fgClk);
padLen = N-length(sig)-delaySamples;
sig = [zeros(1,delaySamples), sig, zeros(1,padLen)];

% PreTrigger Test Signal
% sig = linspace(-1,1,N);

dataCh1.data = sig;
% Clock Signal
% dataCh2.data = repmat([ones(1,10), zeros(1,10)], 1, 16); %5MHZ
% dataCh2.data = repmat([ones(1,10), zeros(1,6)], 1, 7);  %20MHZ
% dataCh2.data = repmat([ones(1,3), zeros(1,2)], 1, 16*4);  %20MHZ
%% Calc clk data
fs     = 20e6;

% The data loaded to the AFG must be a multiplication of 16 in manners of data length.
% creating sClk with 16 cycles meets this requirement no matter
% the fFgClk used.
sClkCycles  = 16;
fsNaive     = fs;
sClkDutyCyc = 50/100;

fgClkSamplesPerSClkCyc           = fgClk/fs;

% If there is no round number of fgClk samples in sClk cycle, 
% slow the sClk to the period where fgClk fits in.
if mod(fgClkSamplesPerSClkCyc,1)~=0
    fgClkSamplesPerSClkCyc    =  ceil(fgClkSamplesPerSClkCyc);
    fs                   = fgClk / fgClkSamplesPerSClkCyc;
    fprintf("Notice: The sampling frequence you have chosen cannot be genrated by the AFG. \n the closest sampling frequency is: %d\n", fs);
end

sClkT                = 1/fs;
fgClkSamplesPerFsSig = fgClkSamplesPerSClkCyc * sClkCycles;

%Create sClk Data
cycleData                   = ones(fgClkSamplesPerSClkCyc,1);
dutyCycleIdx                = floor(fgClkSamplesPerSClkCyc*(1-sClkDutyCyc))+1;
cycleData(dutyCycleIdx:end) = 0;
clkData                     = repmat(cycleData, sClkCycles, 1);
dataCh2.data = clkData;
% figure(); plot(clkData)

%%
dataCh1.dataLen = length(dataCh1.data);
dataCh2.dataLen = length(dataCh2.data);

fg.reset();
fg.setProperties(fGenVars.ch{1}, fGenVars.ch{2});
fg.setData(dataCh1, dataCh2);
fg.configChannel(1); % this will enable the channel
fg.configChannel(2); % this will enable the channel