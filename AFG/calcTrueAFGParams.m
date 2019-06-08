function [AFG] = calcTrueAFGParams(AFG)

extClkTrueFreq = AFG.extClkFreq;
Sclk = AFG.Sclk;
freqTrueTrain = AFG.freqTrain;

TExtClk = 1/extClkTrueFreq;
Tsclk   = 1/Sclk;
tTrain  = 1/freqTrueTrain;

tExtClkTrue = TExtClk;

extClkSamples = TExtClk/Tsclk;
extClkSamplesCycle = extClkSamples;

% Make sure that the asked clock can be generated.
if mod(extClkSamples,1)~=0
    extClkSamplesCycle  =  ceil(TExtClk/Tsclk);
    extClkTrueFreq = 1/(extClkSamplesCycle*Tsclk);
    tExtClkTrue       = extClkSamplesCycle*Tsclk; 
    fprintf("Notice: The clock rate you have chosen cannot be genrated by the AFG. the closest rate is: %d\n", extClkTrueFreq);
end

% Make sure that the ext clock is generated using a multiple of 16 sequence
% of data
i=0;
extClkSamples = extClkSamplesCycle;
while (mod(extClkSamples , 16) ~= 0) || (extClkSamples <192)
   extClkSamples = extClkSamples + extClkSamplesCycle;
   i=i+1;
end
extClkRepetitions = extClkSamples/extClkSamplesCycle;

% Make sure the train samples are multiple of the extClk and 16.
% This might generate a train frequency smaller than asked. closed as
% possible.

%calc naive samples in train:
naiveSamplesInTrain = tTrain/Tsclk;
% Make sure that the train is a multiple of (in manner of sclk): Pulses, extClk cycles and 16.
factor = AFG.sinPeriods*(extClkTrueFreq/AFG.freqSin)*extClkSamplesCycle*16;
samplesInTrain = ceil(naiveSamplesInTrain/factor)*factor;

% Make sure that train is multiple of the ext Clock. in manner of Sclk!
rawSamplesInTrain       = ceil(naiveSamplesInTrain/extClkSamplesCycle)*extClkSamplesCycle;
extClkSamplesInTrain = rawSamplesInTrain/extClkSamplesCycle;

% Make sure that in manner of ext clk, will be a multiple of 16 samples in
% one train.
extClkSamplesInTrain = ceil(extClkSamplesInTrain/16)*16;

%number of Sclk samples in one train:
samplesInTrain = extClkSamplesInTrain*extClkSamplesCycle;

%generate The Puls and Train
dt = 1/Sclk;
t = 0:dt:(AFG.sinPeriods/AFG.freqSin)- dt;

pulseData = zeros(samplesInTrain, 1);
pulseData(1:length(t)) = sin(2*pi*AFG.freqSin*t);

if (AFG.debug)
    figure()
    stem(pulseData)
    title('Signal Downloaded to AFG');
end

%Generate The clock

% amp = str2double(AFG.extClkAmp);
offset = str2double(AFG.extClkOff);

clkData = ones(1,extClkSamplesCycle);
clkData(floor(extClkSamplesCycle*(AFG.extClkDcyc/100))+1:end) = 0;
clkData = repmat(clkData, 1, extClkRepetitions)+offset;

if (AFG.debug)
    figure()
    stem(clkData)
    title('Signal Downloaded to AFG');
    
    figure(); 
    plot(linspace(0,AFG.tTrainTrue,samplesInTrain)*1e6, pulseData)
    set(gca, 'FontSize', 20)
    title('Input US pulse')
    xlabel('t[\mus]')
    ylabel('Amplitude')
    xlim([0, 5])
    
end

AFG.trainTrueFreq = 1/(samplesInTrain*Tsclk);
AFG.tTrainTrue = samplesInTrain*Tsclk;
AFG.extClkTrueFreq = extClkTrueFreq;
AFG.tExtClkTrue = tExtClkTrue;

AFG.naiveSamplesInTrain = naiveSamplesInTrain;
AFG.rawSamplesInTrain = rawSamplesInTrain;
AFG.samplesInTrain = samplesInTrain;

AFG.extClkSamples = extClkSamples;
AFG.extClkSamplesCycle = extClkSamplesCycle;
AFG.extClkRepetitions = extClkRepetitions;
AFG.extCLKSamplesInTrain = extClkSamplesInTrain;

AFG.pulseData = pulseData;
AFG.clkData = clkData;






end
    
   