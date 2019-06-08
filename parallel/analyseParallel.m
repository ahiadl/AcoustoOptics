function [rawData, U, I, fftRes, TimeStamp] = analyseParallel(bufferDataOut, SamplingCard, Sample, Plotting)
    analyseTime = tic;
    
    if isempty(bufferDataOut)
        rawData = []; U = []; I = [];fftRes = [];
        TimeStamp(3) = 0;
        TimeStamp(4) = toc(analyseTime);
        return
    end
    
    postProcessingTime = tic;
    [rawData] = rawDataPostProcessingwithchopping(bufferDataOut,SamplingCard, Sample);
    TimeStamp(1) = toc(postProcessingTime);

    reshapeTime = tic;
    U = reshapeRawData(rawData, Sample);      %reshape data
    TimeStamp(2) = toc(reshapeTime);
    
    fftTime = tic;
    [I, fftRes]  = signalProcessing(U, Sample, Plotting);   %fft for each Z pos on this (X,Y) pos
    TimeStamp(3) = toc(fftTime);
    TimeStamp(4) = toc(analyseTime);
end

