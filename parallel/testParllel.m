function [I, params, status] = testParllel(userParams, app)
%-------------------------------
% Init the Instruments
%-------------------------------
status = true;
close all;
% releaseing all connections to innstruments from outside matlab
instrreset 
pause(1)
if isempty(gcp('nocreate'))
%      parpool;
end

% Taking structs out of userPparams
General = userParams.General;
Plotting = userParams.Plotting;
Debug = userParams.Debug;

%------ Arbitrary Function generator -------
fprintf("Initiating the Arbitrary Function Generator\n");

[AFG] = calcTrueClockParams(userParams.AFG);
if (~Debug.dontActivateAFG)
    [status] = AFGinit(AFG);
end
US = userParams.US;
US.freqTrueTrain = AFG.trainTrueFreq; 
userParams.Sample.trueSampleFreq = AFG.extClkTrueFreq;
if (~status); printf("Arbitrary Function Generator is no valid\n"); return; end
fprintf("----------Done INIT AFG------------\n");

%-------------- Zaber Stages ----------------
fprintf("Starting The Stages Control\n");
[Stages, status] = StagesInit(userParams.Stages);
if ~status; printf("Stages init is not valid\n"); return; end
fprintf("----------Done INIT Stages------------\n");

%------------- I/O controller ---------------
fprintf("Starting The I/O controller\n");
IOdev = initIOControl();
fprintf("Done INIT I/O controller\n");

%------------- Dimensions ---------------
fprintf("Calculating Problem Dimensions according to User Params\n");
[Sample, SamplingCard, status] = calcDimensions(US, userParams.Sample, userParams.SamplingCard);
if(~status); return; end
fprintf("----------Done Calculate Dimensions------------\n");

%------------- Sampling Card ----------------
fprintf("Starting The Sampling Card\n");
[boardHandle, SamplingCard, status] = sampleCardInit(SamplingCard);
if(~status); return; end
fprintf("Done INIT The Sampling Card\n");

%-------------------------------
% Calculating and Initializing
%-------------------------------

%initializing
Zlen = Sample.Npos;
Zbar = Sample.Zres*(0:1:Zlen-1);
I = zeros(Stages.Ylen, Stages.Xlen, Zlen);
Sample.Zbar = Zbar;

%Initializing Plots
if (Plotting.plotFFT)
    close all; 
    Plotting.hFFT = plot(Plotting.fftPlotHandle, Sample.frequencyBarShifted, zeros(1,length(Sample.frequencyBarShifted)));
end
if(Plotting.plotSpatialResult)
    Plotting.hPS  = stem(Plotting.spatialPowerSpectrumHandle, Stages.Ypos, zeros(1,length(Stages.Ypos)));
end

if (Plotting.plotPlanesImage)
    cla(Plotting.planeImageHandle);
    Plotting.hISC = imagesc(Plotting.planeImageHandle,...
                                 'XData', Stages.Ypos,...
                                 'YData', Zbar/1e3,...
                                 'CData', zeros(length(Zbar), length(Stages.Ypos)));
    colorbar(Plotting.planeImageHandle)
end

if sum(Debug.Acquisition + Debug.postProcessing + Debug.Reshape + Debug.FFT) > 0 
    if (Debug.Acquisition)
        figAcq = figure();
        set(figAcq, 'Position', [0, 575 900, 425]); %[left bottom width height]
        Plotting.hACQ = stem(zeros(1,SamplingCard.periodsToSample * SamplingCard.samplesPerRecord));
        title("Multiplexed Raw Data")
    end

    if (Debug.postProcessing)
        figPP = figure();
        set(figPP, 'Position', [950, 575, 900, 425]);
        Plotting.hPP = stem(zeros(1,SamplingCard.samplesPerRecord));
        title('1 Record Raw Data')
    end

    if(Debug.Reshape)
        figRS = figure();
        set(figRS, 'Position', [0, 50, 900, 425]);
        Plotting.hRS = stem(zeros(1,Sample.Npulse*SamplingCard.actualPeriodsSampled));
        title("Signal Reshape")
    end

    if(Debug.FFT)
        figFFT = figure();
        set(figFFT, 'Position', [950, 50, 900, 425]);
        Plotting.hFFTdebug = stem(Sample.frequencyBarShifted, zeros(1,Sample.Npulse*SamplingCard.actualPeriodsSampled));
%         hFFTdebug = stem(zeros(1,Sample.Npulse*SamplingCard.periodsToSample));
        title("Signal FFT")
    end
end
%-------------------------------
% Scan and Process
%-------------------------------

% previouseDataOut = [];
% currentDataOut = [];
% 
% analysisCluster = parcluster('analyseCluster');
% sampleCluster = parcluster('sampleCluster');
% 
% jobSample  = createJob(sampleCluster);
% jobAnalyze = createJob(analysisCluster);
% 
% sampleTask   = createTask(jobSample, @acquireDataParallel, 3, {boardHandle,SamplingCard}); % job, handle, output arguments, input arguments 
% analysisTask = createTask(jobAnalyze, @analyse,    5, {previouseDataOut, Sample, SamplingCard, Plotting});

spmd
    loadlibrary C:\Windows\System32\ATSApi.dll
end

fprintf("Starting Scanning and Sampling\n");
app.updateDimensions(Sample, Stages);
app.transferData([], 1, 0);
Action = {'Full Acquisition'; 'Post Processing'; 'Reshape'; 'FFT'; 'Full Analysing'; 'Plotting';'Single Point'};
TimeStamp = zeros(length(Action),1);
fullScanTime = tic;

for Xidx = 1:Stages.Xlen%tmp
    if ~Debug.dontMoveStages && (length(Stages.Zaber_Devices) == 2)
         Zaber_MoveAbs(Stages.Zaber_Devices(1), Stages.Xpos(Xidx));
    end
    for Yidx = 1:Stages.Ylen
        fprintf("------------------Starting New Position on Y axis----------------------\n")
        singlePoint = tic;
        if ~Debug.dontMoveStages
            Zaber_MoveAbs(Stages.Zaber_Devices(2), Stages.Ypos(Yidx)); 
        end
         if ~Debug.dontOpenIO; openIO(IOdev); end
        sampleAndAnalyse = tic;
        parfor i = 1:2
           if i==1
              [currentDataOut{i}, ~, acquireTime{i}] = acquireDataParallel(boardHandle,SamplingCard);
           else
              [rawDataTmp{i}, Utmp{i}, Itmp{i}, fftResTmp{i}, processTime{i}] = analyse(previouseDataOut, SamplingCard, Sample, Plotting)
           end
        end
        if ~Debug.dontOpenIO; closeIO(IOdev); end
        previouseDataOut = currentDataOut{1};
        TimeStamp(1)     = acquireTime{1};
        
        TimeStamp(2:5) = processTime{2}; 
        rawData        = rawDataTmp{2};
        U              = Utmp{2};
        I(Yidx,Xidx,:) = Itmp{2};
        fftRes         = fftResTmp{2};

%         submit(jobSample);
%         submit(jobAnalyze);
% 
%         wait(jobSample);
%         wait(jobAnalyze);
%         
%         sampleOutput = fetchOutputs(jobSample);
%         analyseOutput = fetchOutputs(jobAnalyze);
%         
%         previouseDataOut = sampleOutput{1};
%         TimeStamp(1)     = sampleOutput{3};
%         
%         rawData        = analyseOutput{1};
%         U              = analyseOutput{2};
%         I(Yidx,Xidx,:) = analyseOutput{3};
%         fftRes         = analyseOutput{4};
%         TimeStamp(2:5) = analyseOutput{5};
        
        sampleAndAnalyse = toc(sampleAndAnalyse);
        disp(sampleAndAnalyse);
        
        TimeStamp(6) = DebugPlots(Debug, currentDataOut{1}, rawData, U, fftRes, I(:,Xidx, Plotting.Zcoor));
       
        TimeStamp(7) = toc(singlePoint);
        timingTable  = table(Action, TimeStamp);
        disp(timingTable);
        
        fprintf("Done Scan for (x,y) = (%8.3f,%8.3f)\n",Stages.Xpos(Xidx),Stages.Ypos(Yidx));
        pause(0.1)
    end

    set(Plotting.hISC, 'CData', reshape(I(:,Xidx,:),Stages.Ylen, Sample.Npos)');
    app.transferData(I, 0, Xidx);
    fprintf("Done Scan for Complete Y");
end
releaseBuffers(boardHandle,SamplingCard)
fullScanTime = toc(fullScanTime);
fprintf("Scan Completed Successfully in %8.3f seconds.\n", fullScanTime);

params.Stages       = Stages;
params.SamplingCard = SamplingCard;
params.General      = General;
params.Plotting     = Plotting;
params.AFG          = AFG;
params.Sample       = Sample;

end