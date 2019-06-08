function [I, params, status] = consistency2(userParams, app)
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

% Initializing results arrays
data.phiCh       = zeros(Sample.Zlen, Timing.numOfTimeFrames, Timing.numOfSets, Timing.numOfQuant(end), digitizer.channels);
data.phiQuant    = zeros(Sample.Zlen, Timing.numOfTimeFrames, Timing.numOfSets, Timing.numOfQuant(end));
data.phiSets     = zeros(Sample.Zlen, Timing.numOfTimeFrames, Timing.numOfSets);
data.phiFrame    = zeros(Sample.Zlen, Timing.numOfTimeFrames);
data.phiSetsStd  = zeros(Sample.Zlen, Timing.numOfTimeFrames, Timing.numOfSets);
data.phiFrameStd = zeros(Sample.Zlen, Timing.numOfTimeFrames);

% Initializing Plots
app.initPlotsSinglePoint(Sample.frequencyBar, Timing)

%-------------------------------
% Scan and Process
%-------------------------------

fprintf("%s - Starting Scanning and Sampling\n", datestr(datetime('now')));
Action = {'Full Acquisition'; 'Post Processing'; 'Reshape'; 'FFT'; 'Full Analysing'; 'Plotting';'Single Point'};
TimeStamp = zeros(length(Action),1);
fullScanTime = tic;

fftWidth = 200;
fftIdxShort = fftWidth+1;
data.fftIdx = fftIdxShort;


for i = 1:Timing.numOfTimeFrames
    fprintf("------------------Starting a New Time Frame----------------------\n")
    fprintf('start time: %s\n', datestr(datetime('now')))
    
    digitizer.timeToSample = Timing.timeToSample(i);
    [Sample, digitizer, ~] = calcDimensions(US, Sample, digitizer);
    [digitizer, ~]         = digitizerInit(digitizer);
    
    
    Plotting.fBar = Sample.frequencyBar;
    Plotting.USidx = Sample.USidx;
    Plotting.timeFrameSamples(i) = Sample.samplesPerPos;
    fftIdxToPlot = Sample.USidx-fftWidth:Sample.USidx+fftWidth;
    data.fBar = Sample.frequencyBar(fftIdxToPlot);
    for j=1:Timing.numOfSets
        fprintf("------------------Starting a New Set----------------------\n")
        for k=1:Timing.numOfQuant(i)
            fprintf("set number: %d", k)

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

            [~, ~, data.phiQuant(:,i,j,k), data.phiCh(:,i,j,k,:), fftRes, TimeStamp(2:5)] = analyse(bufferDataOut, digitizer, Sample, Plotting);
            [data.phiSets(:,i,j), data.phiFrame(:,i), data.phiSetsStd(:,i,j), data.phiFrameStd(:,i)] = ...
                                                                                averageResults(data.phiCh(:,i,j,k,:));
            data.fftRes = fftRes(:,fftIdxToPlot,:);
            
            app.transferSinglePositionDataAndPlot( data, [i,j,k] )
            rawDataName = sprintf('T=%.2f-r=%.2f.mat', Timing.timeFrames(i), j);
            if General.saveRawData
                save([rawDataDir,'\', rawDataName],  'fftResReduced', 'fBarReduced', '-v7.3');
            end

            TimeStamp(7) = toc(singlePoint);  
            fprintf("Done Scan for (T,r) = (%8.3f,%8.3f)\n", Timing.timeFrames(i), j);
%             rotate diffuser
            pause(0.1);
        end
    end
%     updatePlots
    releaseBuffers(digitizer)
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