function [Sample, digitizer, status] = calcDimensionsNew( vars )
% sine  - the sinus
% pulse - the # of periods of sine
% train - pulse followed by consecutive zeros
% signal - net samples in phantom during all trains
% meas - samples of net signal together with initial propagation time and
% pre trigger samples
% acquisition - all the samples taken in the buffers

% External Clock Params
vars.extClkSamplesPerCyc = vars.fSclk/vars.fExtClk;
vars.extClkCycles = 16;
vars.fExtClkNaive = vars.fExtClk;

% if there is no round number of sclk in on extClk cycle, slow the extClk
% to the period where sclk fits in.
if mod(vars.extClkSamplesPerCyc,1)~=0
    vars.extClkSamplesPerCyc  =  ceil(vars.extClkSamplesPerCyc);
    vars.fExtClk = vars.fSclk / vars.extClkSamplesPerCyc;
    fprintf("Notice: The clock rate you have chosen cannot be genrated by the AFG. the closest rate is: %d\n", vars.extClkTrueFreq);
end

vars.tExtClk = 1/vars.fExtClk; 
vars.extClkSigSamples = vars.extClkSamplesPerCyc*vars.extClkCycles;

% US Excitation params
% Make sure that the train is a multiple of (in manner of sclk): Pulses, extClk cycles and 16.
vars.SclkSamplesInPulse = vars.cycInPulse*(vars.fSclk/vars.fSin);
vars.fTrainNaive = vars.fTrain;
vars.SclkSamplesInTrain = vars.fSclk/vars.fTrain;
factor = (vars.cycInPulse * (vars.fExtClk/vars.fSin)) * vars.extClkSamplesPerCyc;
vars.SclkSamplesInTrain = ceil(vars.SclkSamplesInTrain/vars.factor)*factor;
vars.fTrain = vars.fSclk / vars.SclkSamplesInTrain;

%Samples
vars.samplesPerSin         = vars.fs/vars.fSin;
vars.samplesPerPulse       = vars.samplesPerSin * vars.cycInPulse;
vars.samplesPerTrain       = vars.fSample/vars.fTrain;
vars.numOfTrains           = ceil(vars.timeToSample * vars.fTrain);
vars.prePhantomSamples     = ceil( (vars.distFromPhantom/vars.c) * vars.fs );
vars.samplesPerSignal      = vars.samplesPerTrain * vars.numOfTrains;
vars.samplesPerMeas        = vars.preTriggerSamples + vars.prePhantomSamples + vars.samplesPerSignal;
vars.samplesPerSignalAllCh = vars.samplesPerSignal * vars.channels;
vars.samplesPerMeasAllCh   = vars.samplesPerMeas * vars.channels;

vars.samplesPerPos   = vars.samplesPerPulse * vars.numOfTrains;
vars.numOfPos        = vars.samplesPerTrain / vars.samplesPerPulse;
vars.samplesPerZAxis = vars.numOfPos * vars.samplesPerPulse;

%Time
vars.Tsine  = 1/vars.fsin;            %[s], duration of one period of the sin
vars.Tpulse = vars.cycInPulse*vars.Tsine;     %[s], duration of the pulse
vars.Ttrain = 1/vars.ftrain;      %[s], duration of one period of the train
vars.dts    = 1/vars.fs; %[s], duration of a single sample

vars.tPulseVec  = (0:1:vars.samplesPerPulse  - 1)*dts;
vars.tTrainMeas = (0:1:vars.samplesPerTrain  - 1)*dts;
vars.tSigVec    = (0:1:vars.samplesPerSignal - 1)*dts;
vars.tMeasVec   = (0:1:vars.samplesPerMeas   - 1)*dts;

%Length
vars.dzs      = vars.c*vars.dts;
vars.sinLen   = vars.c*vars.Tsin;
vars.pulseLen = vars.c*vars.Tpulse;
vars.zRes     = vars.pulseLen; 
vars.trainLen = vars.c*vars.tTrain;

vars.zLen =  vars.samplesPerZAxis * vars.dzs; 
vars.zVec = (0:1:(vars.samplesPerZAxis-1))*vars.dzs;

%Frequency
k = 0:vars.samplesPerPos-1;
frequencyBar = Sample.trueSampleFreq * k / samplesPerPos;
frequencyBarShifted = [frequencyBar(samplesPerPos/2+1:end)-Sample.trueSampleFreq, frequencyBar(1:samplesPerPos/2)] ;
USidx = (US.freqSin/Sample.trueSampleFreq)*samplesPerPos+1;

%Digitizer
vars.bufferSize            = 8*2^20; %set buffer size to be 8MB
vars.bytesPerSample        = 2;
vars.samplesPerBufferAllCh = vars.bufferSize / vars.bytesPerSample;
vars.samplesPerBuffer      = vars.samplesPerBufferAllCh / vars.channels;
vars.numOfBuffers          = ceil(vars.samplesPerMeas / vars.samplesPerBuffer);
vars.samplesPerAcq         = vars.samplesPerBuffer * vars.numOfBuffers;
vars.samplesPerAcqAllCh    = vars.samplesPerAcq * vars.channles;
vars.actualSampledTime     = samplesPerAcq * vars.dts;

fprintf(" Samples Per Train: %d\n Trains Per Buffer Per Channel: %d\n Buffer Per Acquisition: %d\n Trains Per Acquisition Per Channel: %d\n Train Size: %d[B]\n Buffer Size: %d[MB]\n Acquisition Size: %d[MB]\n Single Train Time: %d[us]\n Single Buffer Time: %d[ms]\n Acquisition Time: %d[s]\n",...
                SamplingCard.samplesPerTrain, SamplingCard.trainsPerBuffer, SamplingCard.buffersPerAcquisition, SamplingCard.trainsPerAcquisition, SamplingCard.bytesPerTrain, SamplingCard.bytesPerBuffer/2^20, SamplingCard.bytesPerAcquisition/2^20, SamplingCard.timePerTrain*1e6, SamplingCard.timePerBuffer*1e3, SamplingCard.trueTimeSampled);
% [digitizer, vars, status] = calcSampleCardParamsNew(digitizer);

end