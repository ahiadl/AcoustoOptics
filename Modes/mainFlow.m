function [I, params, status] = mainFlow(userParams, app)
%-------------------------------
% Init the Instruments
%-------------------------------
status = true;
close all;
% releaseing all connections to innstruments from outside matlab
instrreset 
pause(1)
instrreset
% curPool = gcp('nocreate');
% if isempty(curPool)
%      curPool = parpool;
% end

% Taking structs out of userPparams
General = userParams.General;
Plotting = userParams.Plotting;
Debug = userParams.Debug;
General.Parallel = false;
Debug.Parallel = General.Parallel;
%------ Arbitrary Function generator -------
fprintf("Initiating the Arbitrary Function Generator\n");

[AFG] = calcTrueAFGParams(userParams.AFG);
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
if ~Debug.dontOpenIO; IOdev = initIOControl(); closeIO(IOdev); end
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
I = zeros(Stages.Ylen, Stages.Xlen, Sample.Zlen);
previouseDataOut = [];
% currentDataOut = [];
% rawData = [];
% U = [];
%init parallel job of analysis
% analysisCluster = parcluster('analyseCluster');
% curPool = gcp();
% if (isempty(curPool) ||(isempty(curPool) && strcmp(curPool.cluster.Profile, 'analysisCluster')))
%     poolAnalysis = parpool(analysisCluster);
% else
%     poolAnalysis = curPool;
% end
    
% jobAnalyze = createJob(analysisCluster);
% analysisTask = createTask(jobAnalyze, @analyse,    5, {previouseDataOut, Sample, SamplingCard, Plotting});

% Initializing Plots
Plotting = initPlots(Sample, SamplingCard, Debug, Stages, Plotting);
Plotting.fBar = Sample.frequencyBar;
Plotting.USidx = Sample.USidx;
Plotting.channels = SamplingCard.channels;

% % Initializing Parallel Mechanism
% funcList = {@acquireDataTS, @analyseParallel};
% results{1}  = cell(1,nargout(@acquireDataTS));
% results{2}  = cell(1,nargout(@analyseParallel));
% dataList{1} = cell(1,nargin(@acquireDataTS));
% dataList{2} = cell(1,nargin(@analyseParallel));
% dataList{1} = {boardHandle, SamplingCard};
% dataList{2} = {previouseDataOut, SamplingCard, Sample, Plotting};

%-------------------------------
% Scan and Process
%-------------------------------

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
       
        if General.Parallel
            
            sampleAndAnalyse = tic();
            parfor idx = 1 : length(funcList)
                 [results{idx}{:}] = funcList{idx}(dataList{idx}{:});
            end
            sampleAndAnalyse = toc(sampleAndAnalyse);
            disp(sampleAndAnalyse);
            
            TimeStamp(2:5) = results{2}{5};
            dataList{2}{1} = results{1}{1};
            TimeStamp(6) = DebugPlots(Debug, Plotting, results{1}{1}, results{2}{1:4});
        else
            % Sampling
            if ~Debug.dontOpenIO; openIO(IOdev); end

            fullAcqTime = tic;
            if (SamplingCard.Mode.TS)
                [bufferDataOut, status] = acquireDataTS(boardHandle, SamplingCard);
            else
                [bufferDataOut, status] = acquireDataNPT(boardHandle, SamplingCard);
            end
            if(~status); return; end
            TimeStamp(1) = toc(fullAcqTime);

            if ~Debug.dontOpenIO; closeIO(IOdev); end
                
            % Data Analysis
            [rawData, U, I(Yidx,Xidx,:), fftRes, TimeStamp(2:5)] = analyse(bufferDataOut, SamplingCard, Sample, Plotting);
                TimeStamp(6) = DebugPlots(Debug, Plotting, bufferDataOut, rawData, U, I(:,Xidx, Plotting.Zcoor), fftRes);   
        end
 
        TimeStamp(7) = toc(singlePoint);
        timingTable  = table(Action, TimeStamp);
        disp(timingTable);
        
        fprintf("Done Scan for (x,y) = (%8.3f,%8.3f)\n",Stages.Xpos(Xidx),Stages.Ypos(Yidx));
        pause(0.1);
        save('D:\Results\', 'bufferDataOut', 'U')
        
    end
    if (Plotting.plotPlanesImage)
        set(Plotting.hISC,  'CData', flip(permute(I(:,Xidx,:), [3,1,2]),1));
    end
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