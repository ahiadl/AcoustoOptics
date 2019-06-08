function [digitizer, status] = calcSampleCardParams(digitizer)
AlazarDefs;
status = true;
if ~alazarLoadLibrary(); fprintf('Error: ATSApi library not loaded\n'); status = false; return; end

systemId = int32(1);
boardId = int32(1);

boardHandle = AlazarGetBoardBySystemID(systemId, boardId);
setdatatype(boardHandle, 'voidPtr', 1, 1);
[retCode, ~, ~, bitsPerSample] = AlazarGetChannelInfo(boardHandle, 0, 0); 
if retCode ~= ApiSuccess; fprintf('Error: AlazarGetChannelInfo failed -- %s\n', errorToText(retCode)); status = false; return; end

bytesPerSample = floor((double(bitsPerSample) + 7) / double(8));
maxBufferSize = 8*2^20; % 8[MB]
channels = digitizer.channels;

if digitizer.Mode.NPT
    samplesPerRecord = digitizer.preTriggerSamples + digitizer.postTriggerSamples; %calculated in calcDimensions
    recordsPerAcquisition = 2^ceil(log2(digitizer.periodsToSample));            %user Defined. should be power of 2
    totalSamplingSize =...
        samplesPerRecord * recordsPerAcquisition * digitizer.channels * bytesPerSample;

    bufferSize = 0;
    buffersPerAcquisition = 0.5;

    while bufferSize/(2^20) > 8.5 || bufferSize/(2^20) < 1 
        buffersPerAcquisition = buffersPerAcquisition * 2;
        bufferSize = totalSamplingSize / buffersPerAcquisition;
        recordsPerBuffer = recordsPerAcquisition / buffersPerAcquisition;
    end

    samplesPerBuffer  = samplesPerRecord * recordsPerBuffer * digitizer.channels;
    bytesPerBuffer    = bytesPerSample * samplesPerBuffer;
    
    digitizer.samplesPerTrain                 = samplesPerRecord;
    digitizer.samplesPerBufferPerChannel      = samplesPerBuffer/channels;
    digitizer.samplesPerBuffer                = samplesPerBuffer;
    digitizer.samplesPerAcquisitionPerChannel = (samplesPerBuffer/channels) * buffersPerAcquisition;
    digitizer.samplesPerAcquisition           = samplesPerBuffer  * buffersPerAcquisition;
    
    digitizer.trainsPerAcquisition       = recordsPerBuffer * buffersPerAcquisition; %[per channel]
    digitizer.trainsPerBuffer            = recordsPerBuffer;
    digitizer.buffersPerAcquisition      = buffersPerAcquisition;
    
    digitizer.bytesPerTrain              = samplesPerRecord * bytesPerSample; 
    digitizer.bytesPerBuffer             = bytesPerBuffer;
    digitizer.bytesPerAcquisition        = bytesPerBuffer * buffersPerAcquisition;

    digitizer.totalSamplingSize     = totalSamplingSize;
    digitizer.recordsPerAcquisition = recordsPerAcquisition;
    
    % Due To BUG in the sampling card:
    digitizer.actualPeriodsSampled = recordsPerAcquisition - buffersPerAcquisition;
    
elseif digitizer.Mode.TS
    digitizer.recordsPerAcquisition      = 1;

    samplesPerTrain = digitizer.samplesPerTrain;
    rawTrainsPerAcquisition = ceil(digitizer.timeToSample * digitizer.freqSample/samplesPerTrain); % make sure a complete number of train cycles is in the acquisition
    maxTrainPerBuffer = floor(maxBufferSize / (bytesPerSample*samplesPerTrain*channels)); % how many trains can be fitted into a buffer of 8MB
    buffersPerAcquisition = ceil(rawTrainsPerAcquisition/maxTrainPerBuffer); % how many buffers will be needed
    trainsPerAcquisition = buffersPerAcquisition * maxTrainPerBuffer; % make sure that the number of trains will be fir into complete bumber of 8MB buffers
    
    samplesPerBufferPerChannel = maxTrainPerBuffer * samplesPerTrain;
    bytesPerBuffer = samplesPerBufferPerChannel * channels * bytesPerSample;
    totalAcquisitionSize = buffersPerAcquisition * bytesPerBuffer;
    
    digitizer.samplesPerTrain                 = samplesPerTrain;
    digitizer.samplesPerBufferPerChannel      = samplesPerBufferPerChannel;
    digitizer.samplesPerBuffer                = samplesPerBufferPerChannel * channels;
    digitizer.samplesPerAcquisitionPerChannel = samplesPerBufferPerChannel * buffersPerAcquisition;
    digitizer.samplesPerAcquisition           = samplesPerBufferPerChannel * buffersPerAcquisition * channels;
    
    digitizer.rawTrainsPerAcquisition    = rawTrainsPerAcquisition;
    digitizer.trainsPerAcquisition       = trainsPerAcquisition; %[per channel]
    digitizer.trainsPerBuffer            = maxTrainPerBuffer;
    digitizer.buffersPerAcquisition      = buffersPerAcquisition;
    
    digitizer.bytesPerTrain              = samplesPerTrain * bytesPerSample; 
    digitizer.bytesPerBuffer             = bytesPerBuffer;
    digitizer.bytesPerAcquisition        = totalAcquisitionSize;
    
    digitizer.timePerTrain               = samplesPerTrain / digitizer.freqSample;
    digitizer.timePerBuffer              = samplesPerBufferPerChannel / digitizer.freqSample;
    digitizer.trueTimeSampled            = samplesPerBufferPerChannel * buffersPerAcquisition / digitizer.freqSample ;
    % Due To BUG in the sampling card:
    digitizer.actualPeriodsSampled = trainsPerAcquisition - 1; %buffersPerAcquisition;
    
    digitizer.netTrains =  ceil(digitizer.timeToSample / digitizer.timePerTrain);
    digitizer.netSamples = digitizer.netTrains*samplesPerTrain;
    
end

digitizer.bufferCount = 4;
% Due To BUG in the sampling card:
% SamplingCard.actualPeriodsSampled = recordsPerAcquisition - buffersPerAcquisition;
digitizer.cardPreSampling = 15; %14 to 4 ch, 24 to 1 ch.
end

