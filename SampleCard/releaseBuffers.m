function [] = releaseBuffers(Digitizer)
% Release the buffers
AlazarDefs
 
for bufferIndex = 1:Digitizer.bufferCount
    pbuffer = Digitizer.buffers{1, bufferIndex};
    retCode = AlazarFreeBuffer(Digitizer.boardHandle, pbuffer);
    if retCode ~= ApiSuccess
        fprintf('Error: AlazarFreeBuffer failed -- %s\n', errorToText(retCode));
    end
    clear pbuffer;
end
end

