function [Sample, digitizer, status] = calcDimensions(US,  Sample, digitizer)
% sine  - the sinus
% pulse - the # of periods of sine
% train - pulse followed by consecutive zeros  

Tsine     = 1/US.freqSin;            %[s], duration of one period of the sin
Tpulse    = US.sinPeriods*Tsine;     %[s], duration of the pulse
Ttrain    = 1/US.freqTrueTrain;      %[s], duration of one period of the train
Tsample   = 1/Sample.trueSampleFreq; %[s], duration of a single sample

Npulse    = ceil(Tpulse/Tsample);    %[#], Amount of samples in one pulse duration
Ntrainraw = Ttrain/Tsample;          %[#], Amount of samples in one train cycle duration, raw because can be not a multiplier of Npulse
Ntrain    = floor(Ntrainraw/16)*16;  % this makes sure that Ntrain divided 16 - sample card demand
Npos      = floor(Ntrain/Npulse);    %[#], Amount of positions we can yield in Z axis

Zres = Npulse*Tsample*US.speed; %[m], possible resolution in Z axis 
Zlen = Npos;
Zbar = Zres*(0:1:Zlen-1);

% Calculating sample card parameters
digitizer.preTriggerSamples  = 0;
digitizer.postTriggerSamples = Ntrain; %NPT
digitizer.samplesPerTrain    = Ntrain; % TS
[digitizer, status] = calcSampleCardParams(digitizer);

% Chopping due to geometry and Cards Bug
trainLength          = Ntrain*Tsample*US.speed;
samplesDelay         = ceil((Sample.transducerHeight - Sample.phantomHeight)/(US.speed*Tsample));
samplesResidue       = Ntrain - mod(samplesDelay, Ntrain);
prePhantomCycles     = ceil(samplesDelay/Ntrain);
preTriggerSamples    = digitizer.cardPreSampling; %BUG in the Card
actualPeriodsSampled = digitizer.actualPeriodsSampled - prePhantomCycles;

%Calculating and Translating frequencies to indexes According to the DFT
%theorm:

samplesPerPos = Npulse*actualPeriodsSampled; %after reshaping
% samplesPerPos = Npulse * samplesPer;
k = 0:samplesPerPos-1;
frequencyBar = Sample.trueSampleFreq * k / samplesPerPos;
frequencyBarShifted = [frequencyBar(samplesPerPos/2+1:end)-Sample.trueSampleFreq, frequencyBar(1:samplesPerPos/2)] ;
USidx = (US.freqSin/Sample.trueSampleFreq)*samplesPerPos+1;

%Export Dimensions:
Sample.Npulse    = Npulse;
Sample.Ntrainraw = Ntrainraw;
Sample.Ntrain    = Ntrain;
Sample.Npos      = Npos;

Sample.Zres = Zres;
Sample.Zbar = Zbar;
Sample.Zlen = Zlen;

Sample.Tsine   = Tsine;     
Sample.Tpulse  = Tpulse;    
Sample.Ttrain  = Ttrain;
Sample.Tsample = Tsample;

Sample.trainLength          = trainLength;
Sample.actualPeriodsSampled = actualPeriodsSampled;
Sample.prePhantomCycles     = prePhantomCycles;
Sample.samplesDelay         = samplesDelay;
Sample.samplesResidue       = samplesResidue;
Sample.preTriggerSamples    = preTriggerSamples;

Sample.frequencyBar = frequencyBar;
Sample.frequencyBarShifted = frequencyBarShifted;
Sample.USidx = USidx;
Sample.samplesPerPos = samplesPerPos;

end