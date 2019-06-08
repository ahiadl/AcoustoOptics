function [finalRawData] = rawDataPostProcessing(bufferDataOut,SamplingCard, Sample)

rawData = [];

% Separate Multiplexed Raw Data to Cahnnels
for buf = 1: SamplingCard.buffersPerAcquisition
%     oneBufferChannelSeperate = tic;
    tmpRawData = uint16([]);
    
    %seperate channels
    for i = 1:SamplingCard.channels
        tmpRawData(:,:,i) = bufferDataOut(buf,i:SamplingCard.channels:end);
    end
    
    %chop margins according to card bug
    tmpRawData(:,(end-(SamplingCard.samplesPerTrain-(SamplingCard.cardPreSampling+1))):end,:) = [];
    tmpRawData(:,1:SamplingCard.cardPreSampling,:) = [];
    
    %reshape to channels dims and concat
    if buf == 1
        tmpRawData(:,(end - Sample.samplesResidue):end, :) = [];
        tmpRawData(:,1:Sample.samplesDelay,:) = [];
        tmpRawData = reshape(tmpRawData, SamplingCard.samplesPerTrain, SamplingCard.trainsPerBuffer-1-prePhantomCycles, SamplingCard.channels);
    else
        tmpRawData = reshape(tmpRawData, SamplingCard.samplesPerTrain, SamplingCard.trainsPerBuffer-1, SamplingCard.channels);
    end
        rawData = cat(2,rawData,tmpRawData);
%     oneBufferChannelSeperate = toc(oneBufferChannelSeperate)
end

finalRawData = zeros(size(rawData));
% Convert rawData from Samples Code to Volts
% convertFormatTime = tic;
for ch =1:SamplingCard.channels
    finalRawData(:,:,ch) = convertUnsignedSampleToVolts(rawData(:,:,ch), SamplingCard.voltsRange);
end
% convertFormatTime = toc(convertFormatTime)

end

