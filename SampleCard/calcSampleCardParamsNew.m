function [vars, status] = calcSampleCardParamsNew(digitizer, vars)

vars.bufferSize            = 8*2^20; %set buffer size to be 8MB
vars.bytesPerSample        =  2;
vars.samplesPerBufferAllCh = vars.bufferSize / vars.bytesPerSample;
vars.samplesPerBuffer      = vars.samplesPerBufferAllCh / vars.channels;
vars.numOfBuffers          = ceil(vars.samplesPerMeas / vars.samplesPerBuffer);
vars.samplesPerAcq         = vars.samplesPerBuffer * vars.numOfBuffers;
vars.samplesPerAcqAllCh    = vars.samplesPerAcq * vars.channles;
vars.actualSampledTime     = samplesPerAcq * vars.dts;

digitizer.bufferCount = 4;
digitizer.bufferPerAcquisition = vars.numOfBuffers;

% Due To BUG in the sampling card:
% SamplingCard.actualPeriodsSampled = recordsPerAcquisition - buffersPerAcquisition;
% digitizer.cardPreSampling = 15; %14 to 4 ch, 24 to 1 ch.
end

