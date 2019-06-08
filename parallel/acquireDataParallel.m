function [bufferDataOut, status, TimeStamp] = acquireDataParallel(boardHandle,SamplingCard)
%     if ~Debug.dontOpenIO; openIO(IOdev); end
%     if ~libisloaded('ATSApi.lib')
%         loadlibrary('ATSApi.lib')
%     end
    fullAcqTime = tic;
    [bufferDataOut, status] = acquireData(boardHandle, SamplingCard);
    TimeStamp = toc(fullAcqTime);

%     if ~Debug.dontOpenIO; closeIO(IOdev); end
end

