function [RawData] = acquireDCLevel(app)
AlazarDefs

acquisitionLength_sec = 10e-3;

% TODO: Select the number of samples in each DMA buffer
samplesPerBufferPerChannel = 1e6;
samplesPerSec = 10e6;

samplesPerAcquisition = uint32(floor((samplesPerSec * acquisitionLength_sec + 0.5)));
buffersPerAcquisition = uint32(floor((samplesPerAcquisition + samplesPerBufferPerChannel - 1) / samplesPerBufferPerChannel));

% TODO: Select which channels to capture (A, B, or both)
channelMask = CHANNEL_A + CHANNEL_B + CHANNEL_C + CHANNEL_D;
channelCount = 4;
bufferDataOut = uint16(zeros(1, samplesPerBufferPerChannel));

% Get the sample and memory size
[retCode, app.boardHandle, maxSamplesPerRecord, bitsPerSample] = AlazarGetChannelInfo(app.boardHandle, 0, 0);
if retCode ~= ApiSuccess
    fprintf('Error: AlazarGetChannelInfo failed -- %s\n', errorToText(retCode));
    return
end

bufferCount = uint32(4);

% Calculate the size of each buffer in bytes
bytesPerSample = floor((double(bitsPerSample) + 7) / double(8));
samplesPerBuffer = samplesPerBufferPerChannel * channelCount;
bytesPerBuffer = bytesPerSample * samplesPerBuffer;

retCode = AlazarBeforeAsyncRead(app.boardHandle, channelMask, 0, samplesPerBufferPerChannel, 1, hex2dec('7FFFFFFF'), ADMA_EXTERNAL_STARTCAPTURE + ADMA_CONTINUOUS_MODE);
if retCode ~= ApiSuccess
    fprintf('Error: AlazarBeforeAsyncRead failed -- %s\n', errorToText(retCode));
    return
end

buffers = cell(1, bufferCount);
for j = 1 : bufferCount
    pbuffer = AlazarAllocBuffer(app.boardHandle, bytesPerBuffer);
    if pbuffer == 0
        fprintf('Error: AlazarAllocBuffer %u samples failed\n', samplesPerBuffer);
        return
    end
    buffers(1, j) = { pbuffer };
end

% Post the buffers to the board
for bufferIndex = 1 : bufferCount
    pbuffer = buffers{1, bufferIndex};
    retCode = AlazarPostAsyncBuffer(app.boardHandle, pbuffer, bytesPerBuffer);
    if retCode ~= ApiSuccess
        fprintf('Error: AlazarPostAsyncBuffer failed -- %s\n', errorToText(retCode));
        return
    end
end

retCode = AlazarStartCapture(app.boardHandle);
if retCode ~= ApiSuccess
    fprintf('Error: AlazarStartCapture failed -- %s\n', errorToText(retCode));
    return
end

buffersCompleted = 0;
captureDone = false;
success = false;

while ~captureDone

    bufferIndex = mod(buffersCompleted, bufferCount) + 1;
    pbuffer = buffers{1, bufferIndex};

    % Wait for the first available buffer to be filled by the board
    [retCode, app.boardHandle, bufferOut] = ...
        AlazarWaitAsyncBufferComplete(app.boardHandle, pbuffer, 5000);
    if retCode == ApiSuccess
        % This buffer is full
        bufferFull = true;
        captureDone = false;
    elseif retCode == ApiWaitTimeout
        % The wait timeout expired before this buffer was filled.
        % The board may not be triggering, or the timeout period may be too short.
        fprintf('Error: AlazarWaitAsyncBufferComplete timeout -- Verify trigger!\n');
        bufferFull = false;
        captureDone = true;
    else
        % The acquisition failed
        fprintf('Error: AlazarWaitAsyncBufferComplete failed -- %s\n', errorToText(retCode));
        bufferFull = false;
        captureDone = true;
    end

    if bufferFull
        % TODO: Process sample data in this buffer.
          
        setdatatype(bufferOut, 'uint16Ptr', 1, samplesPerBuffer);
        
        bufferRawData = bufferOut.Value;
        
        % Make the buffer available to be filled again by the board
        retCode = AlazarPostAsyncBuffer(app.boardHandle, pbuffer, bytesPerBuffer);
        if retCode ~= ApiSuccess
            fprintf('Error: AlazarPostAsyncBuffer failed -- %s\n', errorToText(retCode));
            captureDone = true;
        end

        % Update progress
        buffersCompleted = buffersCompleted + 1;
        if buffersCompleted >= buffersPerAcquisition
            captureDone = true;
            success = true;
        end

    end % if bufferFull

end % while ~captureDone

retCode = AlazarAbortAsyncRead(app.boardHandle);
if retCode ~= ApiSuccess
    fprintf('Error: AlazarAbortAsyncRead failed -- %s\n', errorToText(retCode));
end

for bufferIndex = 1:bufferCount
    pbuffer = buffers{1, bufferIndex};
    retCode = AlazarFreeBuffer(app.boardHandle, pbuffer);
    if retCode ~= ApiSuccess
        fprintf('Error: AlazarFreeBuffer failed -- %s\n', errorToText(retCode));
    end
    clear pbuffer;
end



RawData = uint16([]);
for i = 1:channelCount
    RawData(:,i) = bufferRawData(i:channelCount:end);
end




end

