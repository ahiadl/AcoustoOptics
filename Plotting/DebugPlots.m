function [TimeStamp] = DebugPlots(Debug, Plotting, currentRawDataOut, rawData, U, I, fftRes)
    plottingTime = tic;
    if(Debug.Acquisition)
        set(Plotting.hACQ, 'Ydata',currentRawDataOut(1,Plotting.channelToPlot:Plotting.channels:end))
    end

    if(Debug.postProcessing)% && ~Debug.Parallel)
        set(Plotting.hPP, 'Ydata',rawData(:,1,Plotting.channelToPlot))
    end

    if(Debug.Reshape)% && ~Debug.Parallel)
        set(Plotting.hRS, 'Ydata', U(Plotting.Zcoor,:,Plotting.channelToPlot)); 
    end

    if(Debug.FFT)
        set(Plotting.hFFTdebug, 'Ydata', fftshift(abs(fftRes))); 
    end
    
%     if (Plotting.plotFFT)
% %         set(Plotting.hFFT,'Ydata',fftshift(abs(fftRes)))
% %         set(Plotting.hFFTMarkers, 'Xdata', Plotting.fBar(Plotting.USidx), 'Ydata', abs(fftRes(Plotting.USidx)), 'Marker', '+', 'MarkerEdgeColor', 'g') 
%         
%         for i = 1:Plotting.channels
%             set(Plotting.hFFT(i),'Ydata',fftshift(abs(fftRes(:,:,i))))
%         end
%         
%         legend(Plotting.fftPlotHandle,{'Ch 1','Ch 2','Ch 3','Ch 4'}, 'Location', 'northeast')
%         
%         for i = 1:Plotting.channels
%             set(Plotting.hFFTMarkers(i), 'Xdata', Plotting.fBar(Plotting.USidx), 'Ydata', abs(fftRes(:,Plotting.USidx,i)), 'Marker', '+', 'MarkerEdgeColor', 'g') 
%         end
%     end
%     
%     if (Plotting.plotSpatialResult)
%         set(Plotting.hPS,'Ydata',I)
%     end
    
    TimeStamp = toc(plottingTime);
    
end

