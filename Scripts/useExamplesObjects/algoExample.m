aoRes = load('D:\Photodiode\PMT-2s\28-Jul-2021 15-15-26-PMT-2s.mat');
rawData = aoRes.resMeas.res(1).rawData;
%%
% Old Algo (for reference)
algo = Algo();

uVars = algo.uVarsCreate();

uVars.fSin       = 1.25e6;              
uVars.fTrain     = 10e3;
uVars.cycInPulse = 1;
uVars.channels   = 1;

uVars.c                 = 1430;  
uVars.phantomDepth      = 4.6e-2;
uVars.distFromPhantom   = 11e-2;
uVars.bufferSizeBytes   = 8*2^20;
uVars.bytesPerSample    = 2;
uVars.preTriggerSamples = 14;
uVars.fExtClk           = 5e6; %fs
uVars.fSclk             = 100e6;
uVars.timeToSample      = 0.25;
uVars.extClkDcyc        = 50; % [%]
uVars.envDC             = 100e3;
uVars.envUS             = 78e3;
uVars.envHar            = 200e3;
uVars.useQuant          = true;
uVars.quantTime         = 0.002;
uVars.useGPU            = true;
uVars.useHadamard       = true;
uVars.highResAO         = false;
uVars.longMeas          = false;

uVars.export.rawData        = false;
uVars.export.netSignal      = false;
uVars.export.deMultiplexed  = false;
uVars.export.reshapedSignal = false;
uVars.export.fftRes         = false;
uVars.export.usCompCmplx    = false;

algoVars = algo.updateAlgoUserVars(uVars);
algo.setRawData(rawData);
res = algo.analyse();

figure();
subplot(1,3,1)
stem(res.phi)
subplot(1,3,2)
stem(aoRes.resMeas.res(1).phi)
subplot(1,3,3)
stem(res.rawPhi - mean(res.rawPhi(1:75)))
%% New Algo
% clc
algoNew = AlgoNew();

uVarsNew = algoNew.createUserVars();

% US Signal
uVars.cycPerPulse       = 1;
uVars.fSin              = 1.25e6;              
uVars.fSqnc             = 10e3;
uVars.frameTime         = 0.002;

% Sampling Clk
uVars.fs                = 5e6; %fs
uVars.sClkDcyc          = 50;% [%]
uVars.fgClk             = 100e6;

% Geometry
uVars.c                 = 1430;
uVars.distFromPhantom   = 11e-2;

% Digitizer
uVars.timeToSample      = 2;
uVars.channels          = 1;
uVars.bufferSizeBytes   = 8*2^20;
uVars.bytesPerSample    = 2;
uVars.preTriggerSamples = 14;

% Frequency
uVars.envDC             = 100e3;
uVars.envUS             = 200e3;

% General Operation
uVars.useFrame            = true;
uVars.useGPU              = false;
uVars.useHadamard         = true;
uVars.contHadamard        = false;
uVars.highResAO           = false;
uVars.analyzeSingleCh     = false;
uVars.contSpeckleAnalysis = false;
 
% Export data
uVars.export.rawData        = false;
uVars.export.netSignal      = false;
uVars.export.deMultiplexed  = false;
uVars.export.reshapedSignal = false;
uVars.export.fftRes         = false;
uVars.export.usCompCmplx    = false;

algoNewVars = algoNew.setVars(uVars);

% simData = zeros(1,  algoNewVars.samples.samplesPerSignalRaw*4);
% preCut  = algoNewVars.samples.preSignalSamples;
% sVec    = algoNewVars.hadamard.sMat(1, :);
% sVec = repmat(sVec, algoNewVars.samples.samplesPerPulse, 1);
% sVec = sVec(:)';
% 
% samplesPerSignal = algoNewVars.samples.samplesPerSignalRaw + algoNewVars.samples.samplesPerSqnc;
% samplesPerSqnc = algoNewVars.samples.samplesPerSqnc;
% simData(preCut+1 : (preCut+2*samplesPerSqnc) ) =...
%         repmat(sVec, 1, 2);

fBar = algoNewVars.freq.fBar;
algoNewVars.freq.fUSIdx;
algoNew.setRawData(rawData);
resNew = algoNew.analyse();

figure();
subplot(1,3,1)
stem(resNew.phi)
subplot(1,3,2)
stem(aoRes.resMeas.res(1).phi)
subplot(1,3,3)
stem(resNew.rawPhi - mean(resNew.rawPhi(50:end)))

figure(); 
subplot(2,2,1)
plot(fBar, resNew.fittedFFTNorm(:,1))
title("bkg")
subplot(2,2,2)
plot(fBar, resNew.fittedFFTNorm(:,31))
title("peak")
subplot(2,2,3)
plot(resNew.fittedFFTNorm(:,1))
subplot(2,2,4)
plot(resNew.fittedFFTNorm(:,31))


disp(resNew.SNR.valRaw)
