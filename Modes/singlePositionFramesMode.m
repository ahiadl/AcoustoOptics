function [I, params, status] = singlePositionFramesMode(userParams, app)
%-------------------------------
% Init the Instruments
%-------------------------------
status = true;
close all;
% releaseing all connections to innstruments from outside matlab
instrreset 
pause(1)
instrreset

% Taking structs out of userPparams
General = userParams.General;
Plotting = userParams.Plotting;
Debug = userParams.Debug;
US = userParams.US;
Timing = initTiming(userParams.Timing);

[measDir, rawDataDir, ~] = fileSystemInitFlow(General);
diary([measDir,'\log.txt']); diary('on');
[AFG, Stages, IOdev, digitizer, Sample] = peripheralsInitFlow(userParams);

US.freqTrueTrain = AFG.trainTrueFreq;

%-------------------------------
% Calculating and Initializing
%-------------------------------

%initializing
% power spectrum of vs z, timeFrame, ch and repetitions 
V = zeros(Sample.Zlen, digitizer.channels, Timing.numOfTimeFrames, Timing.repPerTimeFrame, );
%total power spectrum vs z, timeFrame, and repetition
I = zeros(Sample.Zlen, Timing.numOfTimeFrames, Timing.repPerTimeFrame);

% Initializing Plots
% initPlotsSinglePoint(app, Sample.frequencyBarShifted)
app.initPlotsSinglePoint(Sample.frequencyBar, Timing)
% Plotting = initPlotsSinglePosition(Plotting, Timing, Digitizer, Sample);
% Plotting.fBar = Sample.frequencyBar;
% Plotting.USidx = Sample.USidx;
% Plotting.channels = digitizer.channels;
% Plotting.timeFrameSamples = zeros(1,Timing.numOfTimeFrames);

%-------------------------------
% Scan and Process
%-------------------------------

fprintf("%s - Starting Scanning and Sampling\n", datestr(datetime('now')));
% app.updateDimensions(Sample, Stages);
% app.transferData([], 1, 0);
Action = {'Full Acquisition'; 'Post Processing'; 'Reshape'; 'FFT'; 'Full Analysing'; 'Plotting';'Single Point'};
TimeStamp = zeros(length(Action),1);
fullScanTime = tic;


fftWidth = 200;
fftIdxShort = fftWidth+1;

for k = 1:scanNum
    for i = 1:Timing.numOfTimeFrames
    fprintf("------------------Starting New Time Frame----------------------\n")
    fprintf('start time: %s\n', datestr(datetime('now')))
    
    digitizer.timeToSample = Timing.timeFrames(i);
    [Sample, digitizer, ~] = calcDimensions(US, Sample, digitizer);
    [digitizer, ~]         = digitizerInit(digitizer);
    
    Plotting.fBar = Sample.frequencyBar;
    Plotting.USidx = Sample.USidx;
    Plotting.timeFrameSamples(i) = Sample.samplesPerPos;
    fftIdxToPlot = Sample.USidx-fftWidth:Sample.USidx+fftWidth;
    fBarReduced = Sample.frequencyBar(fftIdxToPlot);
    for j=1:Timing.repPerTimeFrame
        fprintf("------------------Starting New Measurement----------------------\n")
        fprintf('start time: %s\n', datestr(datetime('now')))
        singlePoint = tic;
        
        if ~Debug.dontOpenIO; openIO(IOdev); end

        fullAcqTime = tic;
        if (digitizer.Mode.TS)
            [bufferDataOut, status] = acquireDataTSdig(digitizer);
        else
            [bufferDataOut, status] = acquireDataNPTdig(digitizer);
        end
        if(~status); return; end
        TimeStamp(1) = toc(fullAcqTime);
        
        if ~Debug.dontOpenIO; closeIO(IOdev); end
        
        [~, ~, I(:,i,j,k), V(:,:,i,j,k), fftRes, TimeStamp(2:5)] = analyse(bufferDataOut, digitizer, Sample, Plotting);

%         t = (0:1:Sample.samplesPerPos-1)*(1/5e6);
%         fftRes = repmat((2/Sample.samplesPerPos)*fft(tissue(:)*sin(2*pi*(1.25e6)*t),[],2), 1,1,4) + 0.1*rand(Sample.Zlen, Sample.samplesPerPos, 4);            
%         V(:,:,i,j) = abs(fftRes(:,Sample.USidx,:));
%         I(:,i,j) = mean(V(:,:,i,j),2);
            
             
%         TimeStamp(6) = DebugPlots(Debug, Plotting, bufferDataOut, rawData, U, I(:,Xidx, Plotting.Zcoor), fftRes);
        fftResReduced = fftRes(:,fftIdxToPlot,:);
        app.transferSinglePositionDataAndPlot( V, I, fftResReduced, fftIdxShort, i,fBarReduced )
        rawDataName = sprintf('T=%.2f-r=%.2f.mat', Timing.timeFrames(i), j);
        if General.saveRawData
            save([rawDataDir,'\', rawDataName],  'fftResReduced', 'fBarReduced', '-v7.3');
        end
       
        TimeStamp(7) = toc(singlePoint);  
        fprintf("Done Scan for (T,r) = (%8.3f,%8.3f)\n", Timing.timeFrames(i), j);
        pause(0.1);

    end
%     updatePlots
    releaseBuffers(digitizer)
    end
end


fullScanTime = toc(fullScanTime);
fprintf("Scan Completed Successfully in %8.3f seconds.\n", fullScanTime);

params.Stages       = Stages;
params.SamplingCard = digitizer;
params.General      = General;
params.Plotting     = Plotting;
params.AFG          = AFG;
params.Sample       = Sample;

Results.I = I;
Results.V = V;
Results.fftResReduced = fftResReduced;

Results.params = params;

save([measDir,'\paramsAndTotalResults.mat'],'-struct', 'Results');

end