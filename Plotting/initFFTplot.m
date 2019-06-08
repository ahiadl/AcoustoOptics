function [hFFT, hFFTMarkers, hFFTt] = initFFTplot(channels, fftAxes, fBar, f)

    close all; 
    cla(fftAxes,'reset')
%     hFFT = plot(fftPlotHandle, Sample.frequencyBarShifted, zeros(1,length(Sample.frequencyBarShifted)));
%     hold(Plotting.fftPlotHandle, 'on');
%     hFFTMarkers = plot(fftPlotHandle, f, 0, 'g+');
%     fftPlotHandle.XLim = [f-0.0001*f, f+0.0001*f];
%     hold(fftPlotHandle, 'off');
    
    close all; 
    cla(fftAxes,'reset')
    for i = 1:channels
        hFFT(i) = plot(fftAxes, fBar, zeros(1,length(fBar)));
        legendstr{i} = ['Ch ', num2str(i)];
        hold(fftAxes, 'on');
    end

    for i = 1:channels
        hFFTMarkers(i) = plot(fftAxes, f, 0, 'g+');
    end
%     fftAxes.XLim = [f-0.0002*f, f+0.0002*f];
    hold(fftAxes, 'off');
    
    hFFTt = title(fftAxes, 'Current FFT');
    legend(fftAxes, legendstr, 'Location', 'northeast')

end

