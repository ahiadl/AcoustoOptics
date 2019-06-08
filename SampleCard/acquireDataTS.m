function [bufferDataOut, status] = acquireDataTS(boardHandle, SamplingCard)
AlazarDefs

bufferCount = SamplingCard.bufferCount;

status = true;
channels = SamplingCard.channels; %parfor issues;
bufferDataOut = uint16(zeros(SamplingCard.buffersPerAcquisition, SamplingCard.samplesPerBuffer));

preparationsTime = tic;
%-----------------------------
%----- Configure AutoDMA -----
%-----------------------------
admaFlags =  ADMA_EXTERNAL_STARTCAPTURE + ADMA_TRIGGERED_STREAMING; % Select AutoDMA flags as required
retCode = AlazarBeforeAsyncRead(boardHandle, SamplingCard.channelMask, 0, SamplingCard.samplesPerBufferPerChannel, 1, hex2dec('7FFFFFFF'), admaFlags);
if retCode ~= ApiSuccess; fprintf('Error: AlazarBeforeAsyncRead failed -- %s\n', errorToText(retCode)); status = false; return; end

% Post the buffers to the board
for bufferIndex = 1 : bufferCount
    pbuffer = SamplingCard.buffers{1, bufferIndex};
    retCode = AlazarPostAsyncBuffer(boardHandle, pbuffer, SamplingCard.bytesPerBuffer);
    if retCode ~= ApiSuccess; fprintf('Error: AlazarPostAsyncBuffer failed -- %s\n', errorToText(retCode)); status = false; return; end
end

%------ Should Start Aquisition From This Point -----%
% Arm the board system to wait for triggers
% startTime = tic;
retCode = AlazarStartCapture(boardHandle);
if retCode ~= ApiSuccess; fprintf('Error: AlazarStartCapture failed -- %s\n', errorToText(retCode)); status = false; return; end
% startTime = toc(startTime);
i=1;
buffersCompleted = 0;
captureDone = false;
buf = 1;
% Wait for sufficient data to arrive to fill a buffer, process the buffer,
% and repeat until the acquisition is complete
captureTime = tic;
while ~captureDone
    oneBuffer = tic;
    bufferIndex = mod(buffersCompleted, bufferCount) + 1;
    pbuffer = SamplingCard.buffers{1, bufferIndex};
    
    % Wait for the first available buffer to be filled by the board
%     waitBufferTime = tic;
    [retCode, boardHandle, bufferOut] = AlazarWaitAsyncBufferComplete(boardHandle, pbuffer, 5000);
%     waitBufferTime = toc( waitBufferTime)
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
%     fprintf("Buffer Status: %d", bufferFull);
    if bufferFull
%         dataTypeTime = tic;
        setdatatype(bufferOut, 'uint16Ptr', 1, SamplingCard.samplesPerBuffer); 
%         dataTypeTime = toc(dataTypeTime)
        % Save Results
        %------------------------------------------------------------------
        % Separate Cahnnels
%         copyDataTime = tic;
        bufferDataOut(buf,:) = bufferOut.Value;
%         copyDataTime = toc(copyDataTime)
        buf = buf +1;
        %------------------------------------------------------------------
        % Make the buffer available to be filled again by the board
%         resetBufferTime = tic;
        retCode = AlazarPostAsyncBuffer(boardHandle, pbuffer, SamplingCard.bytesPerBuffer);
%         resetBufferTimeVec(i) = toc(resetBufferTime);
%         i=i+1;
        if retCode ~= ApiSuccess
            fprintf('Error: AlazarPostAsyncBuffer failed -- %s\n', errorToText(retCode));
            captureDone = true;
        end
        % Update progress
        buffersCompleted = buffersCompleted + 1;
        if buffersCompleted >= SamplingCard.buffersPerAcquisition
            captureDone = true;
        end
    end
    timing(buf) = toc(oneBuffer);
end
captureTime = toc(captureTime);
% disp(timing)
% disp(resetBufferTimeVec)

% abortTime = tic;
% Abort the acquisition
% retCode = AlazarAbortAsyncRead(boardHandle);
% if retCode ~= ApiSuccess; fprintf('Error: AlazarAbortAsyncRead failed -- %s\n', errorToText(retCode)); status = false; return; end
% status = true;

end


