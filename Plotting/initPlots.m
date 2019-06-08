function [Plotting] = initPlots(Sample,Digitizer, Debug, Stages, Plotting)

Plotting.channels = Sample.channels;

% if (Plotting.plotFFT)
%  [Plotting] = initFFTplot(Digitizer,Plotting, Sample);
% end
% 
% if(Plotting.plotSpatialResult)
%     Plotting.hPS  = stem(Plotting.spatialPowerSpectrumHandle, Stages.Ypos, zeros(1,length(Stages.Ypos)));
% end
% 
% if (Plotting.plotPlanesImage)
%     cla(Plotting.planeImageHandle);
%     Plotting.hISC = imagesc(Plotting.planeImageHandle,...
%                                  'XData', Stages.Ypos,...
%                                  'YData', Sample.Zbar/1e3,...
%                                  'CData', zeros(length(Sample.Zbar), length(Stages.Ypos)));
%     colorbar(Plotting.planeImageHandle)
%     axis(Plotting.planeImageHandle,'fill');
% end

if sum(Debug.Acquisition + Debug.postProcessing + Debug.Reshape + Debug.FFT) > 0 
    if (Debug.Acquisition)
        figAcq = figure();
        set(figAcq, 'Position', [0, 575 900, 425]); %[left bottom width height]
        Plotting.hACQ = stem(zeros(1,Digitizer.samplesPerAcquisitionPerChannel));
        title("Multiplexed Raw Data")
    end

    if (Debug.postProcessing)
        figPP = figure();
        set(figPP, 'Position', [950, 575, 900, 425]);
        Plotting.hPP = stem(zeros(1,Digitizer.samplesPerTrain));
        title('Post Processing - 1 Record Raw Data')
    end

    if(Debug.Reshape)
        figRS = figure();
        set(figRS, 'Position', [0, 50, 900, 425]);
        Plotting.hRS = stem(zeros(1,Sample.samplesPerPos));
        title("Signal Reshape")
    end

    if(Debug.FFT)
        figFFT = figure();
        set(figFFT, 'Position', [950, 50, 900, 425]);
        Plotting.hFFTdebug = stem(Sample.frequencyBarShifted, zeros(1,Sample.samplesPerPos));
        title("Signal FFT")
    end
end

end

