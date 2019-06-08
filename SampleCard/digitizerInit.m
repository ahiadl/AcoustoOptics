function [digitizer, status] = digitizerInit(digitizer)

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

status = confExtClk(boardHandle, digitizer.freqSample);
if(~status); return; end

[digitizer.channelMask, status] = confChnls(boardHandle, digitizer.channels);
if(~status); return; end

status = confExtTrig(boardHandle, digitizer.freqSample);
if(~status); return; end

%----------------------------
%--- Configure Dimensions ---
%----------------------------
bufferCount = digitizer.bufferCount;

retCode = AlazarSetRecordSize(boardHandle, digitizer.preTriggerSamples, digitizer.postTriggerSamples);
if retCode ~= ApiSuccess; fprintf('Error: AlazarSetRecordSize failed -- %s\n', errorToText(retCode)); status = false; return; end

if (digitizer.Mode.TS)
    fprintf(" Samples Per Train: %d\n Trains Per Buffer Per Channel: %d\n Buffer Per Acquisition: %d\n Trains Per Acquisition Per Channel: %d\n Train Size: %d[B]\n Buffer Size: %d[MB]\n Acquisition Size: %d[MB]\n Single Train Time: %d[us]\n Single Buffer Time: %d[ms]\n Acquisition Time: %d[s]\n",...
    digitizer.samplesPerTrain, digitizer.trainsPerBuffer, digitizer.buffersPerAcquisition, digitizer.trainsPerAcquisition, digitizer.bytesPerTrain, digitizer.bytesPerBuffer/2^20, digitizer.bytesPerAcquisition/2^20, digitizer.timePerTrain*1e6, digitizer.timePerBuffer*1e3, digitizer.trueTimeSampled);
    admaFlags = ADMA_EXTERNAL_STARTCAPTURE + ADMA_TRIGGERED_STREAMING;
    retCode = AlazarBeforeAsyncRead(boardHandle, digitizer.channelMask, 0, digitizer.samplesPerBufferPerChannel, 1, hex2dec('7FFFFFFF'), admaFlags);
    
else
     fprintf("Number of Buffers: %d\nBuffer's size: %d[MB]\nRecords Per Buffer: %d\nSamples Per Record: %d\nTotal Sampling Size: %d[MB]\n",...
            digitizer.buffersPerAcquisition, digitizer.bytesPerBuffer/2^20, digitizer.recordsPerBuffer, digitizer.samplesPerRecord, digitizer.totalSamplingSize/2^20);
        
     admaFlags = ADMA_EXTERNAL_STARTCAPTURE + ADMA_NPT; % Select AutoDMA flags as required
     retCode = AlazarBeforeAsyncRead(boardHandle, digitizer.channelMask, -int32(digitizer.preTriggerSamples), digitizer.samplesPerRecord, digitizer.recordsPerBuffer, hex2dec('7FFFFFFF'), admaFlags);
     
end

if retCode ~= ApiSuccess; fprintf('Error: AlazarBeforeAsyncRead failed -- %s\n', errorToText(retCode)); status = false; return; end

% ---------------------------------------------
% --- Create pointers to CPU memory Buffers ---
% ---------------------------------------------

% Create an array of DMA buffers on the Card's memory.
digitizer.buffers = cell(1, bufferCount);
for j = 1 : bufferCount
    pbuffer = AlazarAllocBuffer(boardHandle, digitizer.bytesPerBuffer);
    if pbuffer == 0; fprintf('Error: AlazarAllocBuffer %u samples failed\n', digitizer.samplesPerBuffer); status = false; return; end
    digitizer.buffers(1, j) = { pbuffer };
end

digitizer.boardHandle = boardHandle;
    
end

