function [boardHandle, SamplingCard, status] = sampleCardInit(SamplingCard)

%--------------------------------
%--- Init Variables and board ---
%--------------------------------

AlazarDefs;

if ~alazarLoadLibrary(); fprintf('Error: ATSApi library not loaded\n'); status = false; return; end

systemId = int32(1);
boardId = int32(1);

boardHandle = AlazarGetBoardBySystemID(systemId, boardId);
setdatatype(boardHandle, 'voidPtr', 1, 1);

%----------------------------
%--- Configure Clk & Trig ---
%----------------------------

status = confExtClk(boardHandle, SamplingCard.freqSample);
if(~status); return; end

[SamplingCard.channelMask, status] = confChnls(boardHandle, SamplingCard.channels);
if(~status); return; end

status = confExtTrig(boardHandle, SamplingCard.freqSample);
if(~status); return; end

%----------------------------
%--- Configure Dimensions ---
%----------------------------
bufferCount = SamplingCard.bufferCount;

retCode = AlazarSetRecordSize(boardHandle, SamplingCard.preTriggerSamples, SamplingCard.postTriggerSamples);
if retCode ~= ApiSuccess; fprintf('Error: AlazarSetRecordSize failed -- %s\n', errorToText(retCode)); status = false; return; end

if (SamplingCard.Mode.TS)
    fprintf(" Samples Per Train: %d\n Trains Per Buffer Per Channel: %d\n Buffer Per Acquisition: %d\n Trains Per Acquisition Per Channel: %d\n Train Size: %d[B]\n Buffer Size: %d[MB]\n Acquisition Size: %d[MB]\n Single Train Time: %d[us]\n Single Buffer Time: %d[ms]\n Acquisition Time: %d[s]\n",...
    SamplingCard.samplesPerTrain, SamplingCard.trainsPerBuffer, SamplingCard.buffersPerAcquisition, SamplingCard.trainsPerAcquisition, SamplingCard.bytesPerTrain, SamplingCard.bytesPerBuffer/2^20, SamplingCard.bytesPerAcquisition/2^20, SamplingCard.timePerTrain*1e6, SamplingCard.timePerBuffer*1e3, SamplingCard.trueTimeSampled);
    admaFlags = ADMA_EXTERNAL_STARTCAPTURE + ADMA_TRIGGERED_STREAMING;
    retCode = AlazarBeforeAsyncRead(boardHandle, SamplingCard.channelMask, 0, SamplingCard.samplesPerBufferPerChannel, 1, hex2dec('7FFFFFFF'), admaFlags);
    
else
     fprintf("Number of Buffers: %d\nBuffer's size: %d[MB]\nRecords Per Buffer: %d\nSamples Per Record: %d\nTotal Sampling Size: %d[MB]\n",...
            SamplingCard.buffersPerAcquisition, SamplingCard.bytesPerBuffer/2^20, SamplingCard.recordsPerBuffer, SamplingCard.samplesPerRecord, SamplingCard.totalSamplingSize/2^20);
        
     admaFlags = ADMA_EXTERNAL_STARTCAPTURE + ADMA_NPT; % Select AutoDMA flags as required
     retCode = AlazarBeforeAsyncRead(boardHandle, SamplingCard.channelMask, -int32(SamplingCard.preTriggerSamples), SamplingCard.samplesPerRecord, SamplingCard.recordsPerBuffer, hex2dec('7FFFFFFF'), admaFlags);
     
end

if retCode ~= ApiSuccess; fprintf('Error: AlazarBeforeAsyncRead failed -- %s\n', errorToText(retCode)); status = false; return; end

% ---------------------------------------------
% --- Create pointers to CPU memory Buffers ---
% ---------------------------------------------

% Create an array of DMA buffers on the Card's memory.
SamplingCard.buffers = cell(1, bufferCount);
for j = 1 : bufferCount
    pbuffer = AlazarAllocBuffer(boardHandle, SamplingCard.bytesPerBuffer);
    if pbuffer == 0; fprintf('Error: AlazarAllocBuffer %u samples failed\n', SamplingCard.samplesPerBuffer); status = false; return; end
    SamplingCard.buffers(1, j) = { pbuffer };
end

    
end

