function [I, params, status] = fullMeasMode(userParams, app)
%-------------------------------
% Init the Instruments
%-------------------------------
close all;
% releaseing all connections to innstruments from outside matlab
instrreset 
pause(1)
instrreset

% Taking structs out of userPparams
General = userParams.General;
Plotting = userParams.Plotting;
Debug = userParams.Debug;


[measDir, rawDataDir, ~] = fileSystemInitFlow(General);
diary([measDir,'\log.txt']); diary('on');
[AFG, Stages, IOdev, digitizer, Sample] = peripheralsInitFlow(userParams);

%-------------------------------
% Calculating and Initializing
%-------------------------------

%initializing
I = zeros(Stages.Ylen, Stages.Xlen, Sample.Zlen);

% Initializing Plots

Plotting = initPlots(Sample, digitizer, Debug, Stages, Plotting);
% Plotting.fBar = Sample.frequencyBar;
% Plotting.USidx = Sample.USidx;
% Plotting.channels = digitizer.channels;

fftWidth = 100;
fftIdxToPlot = [Sample.USidx-fftWidth:Sample.USidx+fftWidth];
fftIdxShort = fftWidth+1; 
fBar = Sample.frequencyBar(fftIdxToPlot);

% app.updateDimensions(Sample.Zbar, fBar, fftIdxShort, Stages);
app.initScanPlots(Sample.Zbar, fBar, fftIdxShort, Stages)

%-------------------------------
% Scan and Process
%-------------------------------

fprintf("%s - Starting Scanning and Sampling\n", datestr(datetime('now')));

Action = {'Full Acquisition'; 'Post Processing'; 'Reshape'; 'FFT'; 'Full Analysing'; 'Plotting';'Single Point'};
TimeStamp = zeros(length(Action),1);
fullScanTime = tic;


timeOfSimulation = (2.25*digitizer.timeToSample+1) * (Stages.Xlen*Stages.Ylen);

fprintf("Measurement will take approx. %d [s]", timeOfSimulation);
pause(1);

for Xidx = 1:Stages.Xlen%tmp
    if ~Debug.dontMoveStages && (length(Stages.Zaber_Devices) == 2)
         Zaber_MoveAbs(Stages.Zaber_Devices(1), Stages.Xpos(Xidx));
    end
    for Yidx = 1:Stages.Ylen
        fprintf("------------------Starting New Position on Y axis----------------------\n")
        fprintf('start time: %s\n', datestr(datetime('now')))
        singlePoint = tic;
        if ~Debug.dontMoveStages
            Zaber_MoveAbs(Stages.Zaber_Devices(2), Stages.Ypos(Yidx)); 
        end
           
        % Sampling
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

        % Data Analysis
        [rawData, U, I(Yidx,Xidx,:), ~, fftRes, TimeStamp(2:5)] = analyse(bufferDataOut, digitizer, Sample, Plotting);
        TimeStamp(6) = DebugPlots(Debug, Plotting, bufferDataOut, rawData, U, I(:,Xidx, Plotting.Zcoor), fftRes);   
        app.transferScanDataAndPlot(I, fftRes(:,fftIdxToPlot,:), Xidx);
        
        TimeStamp(7) = toc(singlePoint);
        timingTable  = table(Action, TimeStamp);
        disp(timingTable);
        
        fprintf("Done Scan for (x,y) = (%8.3f,%8.3f)\n",Stages.Xpos(Xidx),Stages.Ypos(Yidx));
        pause(1);
        rawDataName = sprintf('X=%.2f-Y=%.2f.mat',Stages.Xpos(Xidx), Stages.Ypos(Yidx));
        if General.saveRawData
            save([rawDataDir,'\', rawDataName], 'U');
        end
        
    end
%     if (Plotting.plotPlanesImage)
%         set(Plotting.hISC,  'CData', flip(permute(I(:,Xidx,:), [3,1,2]),1));
%     end
%     app.transferData(I, 0, Xidx);
    fprintf("Done Scan for Complete Y");
end

releaseBuffers(digitizer)
fullScanTime = toc(fullScanTime);
fprintf("Scan Completed Successfully in %8.3f seconds.\n", fullScanTime);

params.Stages       = Stages;
params.SamplingCard = digitizer;
params.General      = General;
params.Plotting     = Plotting;
params.AFG          = AFG;
params.Sample       = Sample;

Results.I = I;
Results.params = params;

save([measDir,'\paramsAndTotalResults'],'-struct', 'Results');

end