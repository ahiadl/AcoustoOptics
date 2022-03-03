fg = fGen();
fg.connect();
%%

sClk = 100e6;

fGenVars = fg.uVarsCreate();

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



%% Create Signal
% NOTE: number of samples in the data for each channel should be a
% multiplication of 16.

N = 65536;

% AO Signal
fsin = 1.25e6;
delay = 4e-6;

% delay = 0;

sigTime = N/sClk;
% Create a sin
tVecSin  = (0:1:((sClk/fsin)-1))*(1/sClk);
Nsin = length(tVecSin);

tVecSinPad = (0:1:N-1)./sClk;
fVecSig = (sClk/N) *  ( (-N/2) : 1 : (N/2)-1 ); 
sig = sin(2*pi*fsin*tVecSin);

delaySamples = floor(delay * sClk);
padLen = N-length(sig)-delaySamples;
sig = [zeros(1,delaySamples), sig, zeros(1,padLen)];

% PreTrigger Test Signal
% sig = linspace(-1,1,N);

dataCh1.data = sig;
% Clock Signal
dataCh2.data = repmat([ones(1,10), zeros(1,10)], 1, 16); %5MHZ
% dataCh2.data = repmat([ones(1,10), zeros(1,6)], 1, 7);  %20MHZ

%%
dataCh1.dataLen = length(dataCh1.data);
dataCh2.dataLen = length(dataCh2.data);

fg.reset();
fg.setProperties(fGenVars.ch{1}, fGenVars.ch{2});
fg.setData(dataCh1, dataCh2);
fg.configChannel(1); % this will enable the channel
fg.configChannel(2); % this will enable the channel