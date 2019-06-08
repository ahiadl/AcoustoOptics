function [rawData] = rawDataPostProcessingNew(bufferDataOut,digitizer, Sample)

rawData = [];
useGPU = false;
% Separate Multiplexed Raw Data to Cahnnels
% allBufferChannelSeperate = tic;
if (useGPU)
    tmp0 = convertUnsignedSampleToVolts(bufferDataOut, digitizer.voltsRange);
    rawData = reshape(tmp0(:), digitizer.channels,  digitizer.samplesPerAcquisitionPerChannel);
%     rawData = [ch * samplesPerAcquisitionPerChannel]

    %chop unneccesarry cycles
    rawData(:,1:digitizer.cardPreSampling) = []; % Card Bug
    rawData(:,1:Sample.samplesDelay) = [];       % Spread untill first pulse hitting the phantom
    rawData(:,(digitizer.netSamples+1):end) = []; % extra trains to fill a buffer
else
    for buf = 1: digitizer.buffersPerAcquisition
    %     oneBufferChannelSeperate = tic;
        tmpRawData = uint16([]);

        tmpRawData = reshape(bufferDataOut(buf,:), 1, digitizer.channels, digitizer.samplesPerBufferPerChannel);

        if (buf == 1)
            tmpRawData(:,:,1:digitizer.cardPreSampling) = [];
            tmpRawData(:,:,1:Sample.samplesDelay) = [];
        elseif (buf == digitizer.buffersPerAcquisition)
            tmpRawData(:,:,(end-(digitizer.samplesPerTrain-(digitizer.cardPreSampling+1))):end) = [];
            tmpRawData(:,:,(end - Sample.samplesResidue+1):end) = [];
        end
        rawData = cat(3,rawData, tmpRawData);
    %     oneBufferChannelSeperate = toc(oneBufferChannelSeperate)
    end
    
    % permuteTime = tic;
    rawData = permute(rawData, [1 3 2]);
    % permuteTime = toc(permuteTime)

    % finalReshape = tic;
    rawData = reshape(rawData, digitizer.samplesPerTrain, Sample.actualPeriodsSampled, digitizer.channels);
    % finalReshape = toc(finalReshape)
    rawData = convertUnsignedSampleToVolts(rawData, digitizer.voltsRange);
end
% allBufferChannelSeperate = toc(allBufferChannelSeperate)
% Convert rawData from Samples Code to Volts
% convertFormatTime = tic;

% convertFormatTime = toc(convertFormatTime)

end

