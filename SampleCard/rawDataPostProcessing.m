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
    rawData = cat(2,rawData, tmpRawData);
%     oneBufferChannelSeperate = toc(oneBufferChannelSeperate)
end

% chop margins according to card bug
% rawData(:,(end-(SamplingCard.samplesPerTrain-(SamplingCard.cardPreSampling+1))):end,:) = [];
% rawData(:,1:SamplingCard.cardPreSampling,:) = [];

%Chopping the first cycles before the first pulse entered the phantom
rawData(:,(end - Sample.samplesResidue+1):end,:) = [];
rawData(:,1:Sample.samplesDelay,:) = [];
rawData = reshape(rawData, SamplingCard.samplesPerTrain, Sample.actualPeriodsSampled, SamplingCard.channels);

finalRawData = zeros(size(rawData));
% Convert rawData from Samples Code to Volts
% convertFormatTime = tic;
for ch =1:SamplingCard.channels
    finalRawData(:,:,ch) = convertUnsignedSampleToVolts(rawData(:,:,ch), SamplingCard.voltsRange);
end
% convertFormatTime = toc(convertFormatTime)

end

