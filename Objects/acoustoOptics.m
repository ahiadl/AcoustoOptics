classdef acoustoOptics < handle
    %CONSIST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        owner;
        
        % Objects:
        fGen;
        digitizer;
        IO;
        algo;
        fileSystem;
        graphics;

        % Data:
        rawData;
        result;
        timeTable;

        % Variables:
        uVars;    % Reduced from user to acoustoOptics
        extVars;  % Extended from acoustoOptics to submodules
        measVars; % Measurement vars (after calculation) given from submodules to acoustoOptics
        oldExtVars;
        limits
        % Control Vars:
        changeLog;
        connected;
        periAvail;
        runLive;
        graphicsNames;
        owned;
    end
    
    methods (Static) 
        function aoVars = aoVarsCreate()
            % Algo
            aoVars.c                 = [];
            aoVars.fSin              = [];              
            aoVars.fTrain            = [];
            aoVars.cycInPulse        = []; 
            aoVars.phantomDepth      = [];
            aoVars.distFromPhantom   = [];
            aoVars.fExtClk           = []; %fs
            aoVars.timeToSample      = [];
            aoVars.extClkDcyc        = []; % [%]
            aoVars.fastAnalysis      = false;
            aoVars.useQuant          = false;
            aoVars.quantTime         = [];
            aoVars.useHadamard       = false;
            aoVars.envDC             = [];
            aoVars.envUS             = [];
            aoVars.envHar            = [];  
            % algo & fGen
            aoVars.fSclk             = [];     %update in fGen
            
            % algo & digitizer
            aoVars.channels          = []; %update in digitizer
            
            % IO
            aoVars.IOPort            = [];
            aoVars.IOLine            = [];
            
            % General: algo & digitizer
            aoVars.useGPU            = false; %algo, digitizer;
            
            aoVars.useVirtualData    = false;
            aoVars.limitByN          = [];
            aoVars.N                 = [];
            aoVars.dispTimeTable     = true;
            aoVars.skipParamsCheck   = false;
        end
        
        function uVars  = uVarsCreate()
            % Algo and Peripherals
            uVars.ao         = acoustoOptics.aoVarsCreate();
            % Figures
            uVars.figs       = AOGraphics.createUserVars();
            % File system     
            uVars.fileSystem = fileSystemAO.uVarsCreate();
        end
        
        function ouVars = ownerObjUserVarsCreate()
           ouVars.ao         =  acoustoOptics.aoVarsCreate();
           ouVars.figs       =  AOGraphics.createGraphicsUserVars();
           ouVars.fileSystem =  AOFileSystem.aoFSVarsCreate();
        end
        
        function newExtVars = extVarsCreate()
            newExtVars.algo      = Algo.uVarsCreate();
            newExtVars.fGen      = fGen.uVarsCreate();
            newExtVars.digitizer = Digitizer.uVarsCreate();
            newExtVars.IO        = IO.uVarsCreate();
            newExtVars.figs      = AOGraphics.createGraphicsVars();
        end

        function log = checkParamsChanged(oldVars, newVars)
            fields = {'fGen', 'algo', 'digitizer', 'IO'};%, 'graphics'};
            changesFound = false;
            for i = 1:length(fields)
                names = fieldnames(newVars.(fields{i}));
                log.(fields{i}) = false;
                for j = 1:numel(names)
                    if strcmp(names{j}, 'figs')
                        continue;
                    end
                    log.(fields{i}) = log.(fields{i}) ||... 
                        (isempty(oldVars.(fields{i}).(names{j}))) ||...
                        (~isequal(oldVars.(fields{i}).(names{j}), newVars.(fields{i}).(names{j})));
                        
                    if log.(fields{i}) && strcmp(fields{i}, 'algo') && ...
                             (strcmp(names{j}, 'fSin')       || strcmp(names{j}, 'fTrain') ||...
                              strcmp(names{j}, 'cycInPulse') || strcmp(names{j}, 'channels') || strcmp(names{j}, 'useHadamard'))
                        log.fGen = 1;
                    end
                end
                changesFound = changesFound || log.(fields{i});
            end
            if changesFound
                fprintf("Found changes in variables - should re-config peripheral\n") 
            else
                fprintf("No changes in user variables\n")
            end
        end
    end
    
    methods
        
        function this = acoustoOptics()
            fprintf("AOI: ------- Creating AcoustoOptic ----------\n");
            fprintf("AOI: 1. Creating an Arbitrary Function Generator Object\n");
            this.fGen = fGen();
            fprintf("AOI: 2. Creating a Digitizer Object\n");
            this.digitizer = Digitizer();
            fprintf("AOI: 3. Creating an IO Object\n");
            this.IO = IO();
            fprintf("AOI: 4. Creating an Algorithm Object\n");
            this.algo = Algo();
            fprintf("AOI: 5. Creating a File System Object\n");
            this.fileSystem = fileSystemAO(this);
            fprintf("AOI: 6. Creating a Graphics Object\n");
            this.graphics = AOGraphics();
            
            this.uVars   = this.uVarsCreate();
            this.extVars = this.extVarsCreate();
            
            this.limits.measLimit = 10;
            this.algo.setMeasLimit(this.limits.measLimit);
            
            this.graphicsNames = this.graphics.getGraphicsNames();
            
            this.connected = false;
            this.runLive = false;
        end 

        function init(this, owner)
            fprintf("AOI: ------- Initiating AcoustoOptic ----------\n");
            if this.connected   
                fprintf("AOI: Acousto Optics system is already connected, no need to reconnect.\n");
            else
                fprintf("AOI: 1. Resetting Instruments\n")
                instrreset;
                fprintf("AOI: 2. Connecting to fGen\n")
                % Configure and Activate fGen
                this.fGen.connect()
                this.periAvail.fGen = this.fGen.hardwareAvailable;
                fprintf("AOI: 3. Connecting to Digitizer\n")
                % Start the digitizer
                this.digitizer.connect();
                this.periAvail.digitizer = this.digitizer.system.hardwareAvailable; 
                fprintf("AOI: 4. Connecting to IO\n")
                % Start the IO
                this.IO.connect();
                this.periAvail.IO = this.IO.hardwareAvailable;

                this.connected = true;
                
                if nargin >1
                    this.owner = owner;
                    this.owned = true;
                end
            end
            this.periAvail.completeHardwareAvail =  this.periAvail.digitizer && this.periAvail.fGen && this.periAvail.IO;
        end
        
        % Variables handling
        function setMeasVars(this, uVars)
%             user gives reduced variables. this function convert these 
%             variables to extended variables suitable for each
%             submodule\peripheral. this function extend the variables by
%             adding default acousto optics values to the user reduced
%             variables
            fprintf("AOI: ------- Sets Extended Vars From Reduced Vars ----------\n")
            this.uVars = uVars;
            
            fprintf("AOI:  1. Sets Algorithm reduced Vars\n")
            % Algo 
            algoVars      = this.algo.uVarsCreate();

            algoVars.c               = uVars.ao.c;
            algoVars.bufferSizeBytes = 8*2^20;
            algoVars.bytesPerSample  = 2;
            algoVars.fSin            = uVars.ao.fSin;
            algoVars.fTrain          = uVars.ao.fTrain;
            algoVars.cycInPulse      = uVars.ao.cycInPulse;
            algoVars.phantomDepth    = uVars.ao.phantomDepth;
            algoVars.distFromPhantom = uVars.ao.distFromPhantom;
            algoVars.fExtClk         = uVars.ao.fExtClk;

            algoVars.extClkDcyc      = uVars.ao.extClkDcyc;
            algoVars.fastAnalysis    = uVars.ao.fastAnalysis;
            algoVars.useQuant        = uVars.ao.useQuant;
            algoVars.quantTime       = uVars.ao.quantTime;
            algoVars.useHadamard     = uVars.ao.useHadamard;
            algoVars.envDC           = uVars.ao.envDC;
            algoVars.envUS           = uVars.ao.envUS;
            algoVars.envHar          = uVars.ao.envHar;
            
            if uVars.ao.timeToSample > this.limits.measLimit
                AO.splitMeas = true;
                algoVars.longMeas = true;
                [algoVars.timeToSample, AO.splitNum] = this.calcSplitTimeToSample(uVars.ao.timeToSample);
            else 
                AO.splitMeas = false;
                AO.splitNum = 1;
                algoVars.longMeas = false;
                algoVars.timeToSample    = uVars.ao.timeToSample;
            end
            
            algoVars.export.netSignal      = uVars.figs.validStruct.netSignal || uVars.fileSystem.saveNetSignal;
            algoVars.export.deMultiplexed  = uVars.figs.validStruct.deMul || uVars.fileSystem.saveDemultiplexed;
            algoVars.export.reshapedSignal = uVars.figs.validStruct.reshapedSignal || uVars.fileSystem.saveReshapedSignal;
            algoVars.export.fftRes         = uVars.fileSystem.saveFFT;
            algoVars.export.usCompCmplx    = uVars.fileSystem.saveFFT;

            fprintf("AOI:  2. Sets Function Generator reduced Vars\n")
            %fGen - static configurations
            fGenVars      = this.fGen.uVarsCreate();

            fGenVars.ch{1}.daqFilter     = 50;
            fGenVars.ch{1}.amp           = 1;
            fGenVars.ch{1}.bias          = 0;
            fGenVars.ch{1}.triggerOwner  = true;
            fGenVars.ch{1}.triggerWidth  = 512;
            fGenVars.ch{1}.useExtSclkSrc = false;

            fGenVars.ch{2}.daqFilter     = 50;
            fGenVars.ch{2}.amp           = 1;
            fGenVars.ch{2}.bias          = 0;
            fGenVars.ch{2}.triggerOwner  = false;
            fGenVars.ch{2}.triggerWidth  = 16;
            fGenVars.ch{2}.useExtSclkSrc = false;

            % Algo & fGen
            algoVars.fSclk      = uVars.ao.fSclk; %update in fGen
            fGenVars.ch{1}.Sclk = algoVars.fSclk;
            fGenVars.ch{2}.Sclk = algoVars.fSclk; 

            fprintf("AOI:  3. Sets Digitizer reduced Vars\n")
            % Digitizer 
            digitizerVars = this.digitizer.uVarsCreate();

            digitizerVars.mode = 'TS'; 
            digitizerVars.fs   = algoVars.fExtClk;
            
            % Algo & Digitizer
            algoVars.preTriggerSamples  = digitizerVars.preTriggerSamples;
            algoVars.useGPU             = uVars.ao.useGPU;
            algoVars.channels           = uVars.ao.channels;

            digitizerVars.useGPU        = algoVars.useGPU;
            digitizerVars.channels      = algoVars.channels;           
            
            % IO
            IOVars = this.IO.uVarsCreate;

            IOVars.mode = 0; % Output Only
            IOVars.port = uVars.ao.IOPort;
            IOVars.line = uVars.ao.IOLine;

            % Calculate measurement variables and dimensions
            fprintf("AOI:  4. Calculating Algorithm Variables\n");
            this.measVars.algo = this.algo.updateAlgoUserVars(algoVars);

            % Complete digitizer variables from calculated algorithm vars
            fprintf("AOI:  5. Extracting Digitizer Variables from Algorithm\n");
            digitizerVars.numOfBuffers          = this.measVars.algo.digitizer.numOfBuffers;
            digitizerVars.bufferSizeBytes       = this.measVars.algo.digitizer.bufferSizeBytes;
            digitizerVars.samplesPerAcq         = this.measVars.algo.digitizer.samplesPerAcq;
            digitizerVars.samplesPerBuffer      = this.measVars.algo.digitizer.samplesPerBuffer;
            digitizerVars.samplesPerAcqAllCh    = this.measVars.algo.digitizer.samplesPerAcqAllCh;
            digitizerVars.samplesPerBufferAllCh = this.measVars.algo.digitizer.samplesPerBufferAllCh;
           
            % FileSystem Vars
            fprintf("AOI:  6. Extracting File system variables.\n");
            if AO.splitMeas
                fileSystemVars.saveNetSignal        = false;
                fileSystemVars.saveDemultiplexed    = false;
                fileSystemVars.saveReshapedSignal   = false;
                fileSystemVars.saveFFT              = false;
            else
                fileSystemVars.saveNetSignal        = uVars.fileSystem.saveNetSignal;
                fileSystemVars.saveDemultiplexed    = uVars.fileSystem.saveDemultiplexed && uVars.useHadamard;
                fileSystemVars.saveReshapedSignal   = uVars.fileSystem.saveReshapedSignal;
                fileSystemVars.saveFFT              = uVars.fileSystem.saveFFT;
            end
            fileSystemVars.saveRawData          = uVars.fileSystem.saveRawData;
            
            fileSystemVars.splitMeas = AO.splitMeas;
            fileSystemVars.splitNum  = AO.splitNum;
            
            fileSystemVars.saveResults          = uVars.fileSystem.saveResults;
            fileSystemVars.saveFigs             = uVars.fileSystem.saveFigs;
            fileSystemVars.dontSaveVars         = uVars.fileSystem.dontSaveVars;
            
            fileSystemVars.dirPath              = uVars.fileSystem.dirPath;
            fileSystemVars.projName             = uVars.fileSystem.projName;
            fileSystemVars.resDirName           = uVars.fileSystem.resDirName;

            fileSystemVars.extProject            = uVars.fileSystem.extProject;

            this.fileSystem.setUserVars(fileSystemVars);

            % Set Graphics variables
            fprintf("AOI:  7. Sets Graphic Variables.\n")

            figsVars = AOGraphics.createGraphicsVars();
            
            if uVars.figs.zIdx > this.measVars.algo.samples.numOfPos
                figsVars.zIdx = 1;
            else
                figsVars.zIdx = uVars.figs.zIdx;
            end
            
            if uVars.figs.ch > this.measVars.algo.digitizer.channels
                figsVars.ch = 1;
            else
                figsVars.ch = uVars.figs.ch;
            end
            
            if uVars.figs.quant > this.measVars.algo.samples.numOfQuant
                figsVars.quant = 1;
            else
                figsVars.quant = uVars.figs.quant;
            end
            figsVars.reopenFigures   = uVars.figs.reopenFigures;
            figsVars.plotPhiInd      = uVars.figs.plotPhiInd;
            figsVars.displayFullFFT  = uVars.figs.displayFullFFT;
            figsVars.FFTenv          = uVars.figs.FFTenv;
            figsVars.intExt          = uVars.figs.intExt;
            figsVars.validStruct     = uVars.figs.validStruct;
            
            if uVars.ao.fastAnalysis
%                 figsVars.validStruct.FFT     = false;
                figsVars.validStruct.qAvgChFFT = false;
            end
            
            figsVars.extH = uVars.figs.extH;
            
            figsVars.numOfChannels = uVars.ao.channels;
            figsVars.df            = this.measVars.algo.freq.df;
            figsVars.fBar          = this.measVars.algo.freq.frequencyBarShifted;
            figsVars.fIdx          = this.measVars.algo.freq.fIdxPosShift;
            figsVars.zVec          = this.measVars.algo.len.zVecUSRes;
            
            figsVars.tUS      = this.measVars.algo.timing.tUS;
            figsVars.tExtClk  = this.measVars.algo.timing.tExtClk;
            
            if ~AO.splitMeas
                figsVars.tRawData = this.measVars.algo.timing.tAcqVec;
                figsVars.tMeasSig = this.measVars.algo.timing.tMeasVec;
                figsVars.tNetSig  = this.measVars.algo.timing.tSigVec;
                figsVars.tPos     = this.measVars.algo.timing.tPosVec; 
            else
                figsVars.validStruct.fullSignal     = false;
                figsVars.validStruct.measSamples    = false;
                figsVars.validStruct.netSignal      = false;
                figsVars.validStruct.deMul          = false;
                figsVars.validStruct.reshapedSignal = false;
                
                figsVars.tRawData = [];
                figsVars.tMeasSig = [];
                figsVars.tNetSig  = [];
                figsVars.tPos     = [];
            end
                
            figsVars.zAxis    = this.measVars.algo.len.zVecUSRes;
            
            this.measVars.figs = figsVars;
            this.graphics.setUserVars(figsVars);
            
            % Writing Extended Vars to Acousto Optics Object
            fprintf("AOI:  8. Set AO Object Variables.\n");
            % Set AO vars
            AO.useVirtualData   = uVars.ao.useVirtualData;
            AO.limitByN         = uVars.ao.limitByN;
            AO.N                = uVars.ao.N;
            AO.dispTimeTable    = uVars.ao.dispTimeTable;
            AO.skipParamsCheck  = uVars.ao.skipParamsCheck;
            
            fprintf("AOI:  9. Set Extended Vars\n")
            newExtVars.AO         = AO;
            newExtVars.algo       = algoVars;
            newExtVars.fGen       = fGenVars;
            newExtVars.digitizer  = digitizerVars;
            newExtVars.IO         = IOVars;
            newExtVars.fileSystem = fileSystemVars;
            newExtVars.figs       = figsVars;

            % Set extended
            this.oldExtVars = this.extVars;
            this.extVars    = newExtVars;
            
            if AO.skipParamsCheck
                this.changeLog.fGen      = true;
                this.changeLog.algo      = true;
                this.changeLog.digitizer = true;
                this.changeLog.IO        = true;
            else
                this.changeLog = this.checkParamsChanged(this.oldExtVars, this.extVars);
            end
                    
            
            this.measVars.AO = this.extVars.AO;
        end
        
        function aoVars = getVars(this)
%             fprintf("------- Downloading AO Vars ----------\n")
            aoVars.measVars  = this.measVars;
            aoVars.extVars   = this.extVars;
            aoVars.uVars     = this.uVars; 
            
            % We put an empty graphic request in order to decrease the size
            % of the vars mat file.
            aoVars.measVars.figs = AOGraphics.createGraphicsVars();
            aoVars.extVars.figs  = AOGraphics.createOwnerVars();
            aoVars.uVars.figs    = AOGraphics.createUserVars();
            
            aoVars.measVars.algo.timing.tSigVec  = [];
            aoVars.measVars.algo.timing.tMeasVec = [];
            aoVars.measVars.algo.timing.tAcqVec  = [];
        end
        
        function [splitTime, splitNum] = calcSplitTimeToSample(this, time)
            splitTime = this.limits.measLimit+1;
            splitNum = 0;
            while splitTime > this.limits.measLimit
                splitNum = splitNum+1;
                splitTime = round(time/splitNum, 3);
            end
        end
        
        % Peripherals Configurations
        function configPeripherals(this)
            fprintf("AOI: ------- Configuring Peripheral ----------\n")

            % Config the FS
            fprintf("AOI:  ** Configuring FileSystem\n")
            this.fileSystem.configFileSystem();
            
            % Reconstruct figures
            % Should be before fGen configuration
            fprintf("AOI:  ** Configuring Graphics\n")
            this.graphics.setGraphicsDynamicVars();
            
            % Download signal to fGen and Activate
            if this.changeLog.fGen
                fprintf("AOI:  ** Configuring fGen\n")
                if this.periAvail.fGen
                    this.configfGen();
                end
            end
            
            % Config digitizer
            if this.changeLog.digitizer
                fprintf("AOI:  ** Configuring Digitizer\n")
                if this.periAvail.digitizer
                    this.configDigitizer();
                end
            end
            
            % Config the IO
            if this.changeLog.IO
                fprintf("AOI:  ** Configuring IO\n")
                if this.periAvail.IO
                    this.configIO()
                end
            end

            fprintf("AOI:  ** Done configuring Peripherals\n")
        end
        
        function configfGen(this)
            [sigData, clkData] = this.algo.createSignalsForfGen();
            
            dataCh1 = this.fGen.getChDataStruct();
            dataCh2 = this.fGen.getChDataStruct();
            
            dataCh1.data = sigData;
            dataCh2.data = clkData;
            dataCh1.dataLen = length(dataCh1.data);
            dataCh2.dataLen = length(dataCh2.data);

            this.measVars.fGen = this.extVars.fGen; 
            
            this.fGen.reset();
            this.fGen.setProperties(this.extVars.fGen.ch{1}, this.extVars.fGen.ch{2});
            this.fGen.setData(dataCh1, dataCh2);
            this.fGen.configChannel(1);
            this.fGen.configChannel(2);
            
            this.graphics.setData([], sigData', clkData')
            
            if this.measVars.figs.validStruct.extClk
                this.graphics.plotAFGSignals('extClk');
            end
            
            if this.measVars.figs.validStruct.usSignal
                this.graphics.plotAFGSignals('usSignal');
            end
        end
    
        function configDigitizer(this)
            this.digitizer.setDigitizerVars(this.extVars.digitizer);
            this.digitizer.configure();
            this.measVars.digitizer = this.digitizer.getVars;
        end
        
        function configIO(this)
           this.IO.allocPorts(this.extVars.IO);
           this.measVars.IO = this.extVars.IO;
        end

        function resetDigitizerMem(this)
         	this.digitizer.releaseBuffers();  
        end

        % Measurements algorithms types     
        function res = runAcoustoOptics(this)
            % Save AO variables (if needed, decided by fileSystem)
            this.fileSystem.saveVarsToDisk();            
            % Measure data and analyse it with AO Algorithm
            if this.measVars.AO.splitMeas
                resArr(this.measVars.AO.splitNum) = struct('qAvgChFFT', [], 'unFittedFFT', [], 'unFittedFFTShift', [], 'fitModel', [], 'fittedFFT', [], 'phi', []);
%                 resArr(this.measVars.AO.splitNum) = struct('qAvgChFFT', []);
                for i = 1:this.measVars.AO.splitNum
                    fprintf("AOI: Split %d/%d: \n", i, this.measVars.AO.splitNum);
                    this.fileSystem.updateSplitInd(i);
                    resArr(i) = this.measureAndAnalyse();
                    this.graphics.setData(resArr(i));
                    this.plotAll();
                end
                this.result = this.algo.analyseSplittedData( resArr, this.measVars.AO.splitNum);
                this.result.splitedRes = resArr;
            else 
                this.result = this.measureAndAnalyse();
            end
            
            % Set data to graphics object
            this.timeTable.setResultsToGraphics = tic;
            this.graphics.setData(this.result);
            this.timeTable.setResultsToGraphics = toc(this.timeTable.setResultsToGraphics);
            
            % Plot Results (if needed, decided by AO)
            this.timeTable.plotAll = tic;
            this.plotAll();
            this.timeTable.plotAll = toc(this.timeTable.plotAll);
            
            % Save results (if needed, decided by fileSystem)
            % fprintf("AOI: Saving Results To Disk\n");
            this.timeTable.saveData = tic;
            this.fileSystem.saveResultsToDisk(this.result)
            this.fileSystem.closeFileSystem();
            this.timeTable.saveData = toc(this.timeTable.saveData);
            
            % Collect time statistics
            if this.measVars.AO.dispTimeTable
                this.updateTimeTable();
            end
            res = this.result;
            fprintf ("AOI: Done AO\n") 
        end
        
        function res = liveAcoustoOptics(this) 
            if this.measVars.AO.limitByN
                this.fileSystem.initLiveAOFS(this.getVars());
                this.runLive = true;
                if this.measVars.AO.splitMeas
                    res(this.measVars.AO.N) = struct('splitedRes', [], 'qAvgChFFT', [], 'unFittedFFT', [], 'unFittedFFTShift', [], 'fitModel', [], 'fittedFFT', [], 'phi', []);
                else
                    res(this.measVars.AO.N) = struct('qAvgChFFT', [], 'unFittedFFT', [], 'unFittedFFTShift', [], 'fitModel', [], 'fittedFFT', [], 'phi', []);
                end
                for i=1:this.measVars.AO.N
                    this.fileSystem.setLiveAOInd(i);
                    res(i) = this.runAcoustoOptics();
                    pause(0.1)
                    if ~this.runLive
                        fprintf("AOI: Live AO was STOPPED.\n");
                        break;
                    end
                end
                
                this.fileSystem.closeLiveAoFileSystem();
            else %infinite number of measurements
                %saving is not allowed in unlimited live AO
                this.runLive = true;
                while(this.runLive)
                   this.runAcoustoOptics();
                   pause(0.1);
                end
                fprintf("AOI: Live AO was STOPPED.\n");
                res = this.result;
            end           
        end
        
        function res = measureAndAnalyse(this)
           if this.measVars.AO.useVirtualData
               % Generate virtual data (DEBUG)
                bufferDataOut = this.createVirtualData();
            elseif this.periAvail.completeHardwareAvail
                % Get data from digitizer (Measure)
                fprintf ("AOI: Acquiring...\n")
                this.timeTable.acq = tic;
                this.IO.open();
                bufferDataOut      = this.digitizer.acquireDataTS();
                this.IO.close();
                this.timeTable.acq = toc(this.timeTable.acq);

                this.timeTable.moveData = tic;
                this.algo.setRawData(bufferDataOut);
                bufferDataOut           = gather(bufferDataOut);
                this.timeTable.moveData = toc(this.timeTable.moveData);
                
                if this.measVars.AO.splitMeas
                    this.fileSystem.saveRawDataToDisk(bufferDataOut)
                    bufferDataOut = [];
                end
            else
                 % No source of data (Error)
                 fprintf("AOI: No source of data found.\n");  
                 res = [];
                 return;
            end
            
            fprintf ("AOI: Analyzing!\n")
            this.timeTable.analyse = tic;
            res                    = this.algo.analyse();
            this.result            = res;
            if ~this.measVars.AO.splitMeas
                this.result.rawData    = bufferDataOut;
            end
            this.timeTable.analyse = toc(this.timeTable.analyse);
        end
        
        % Misc
        function bufferDataOut = createVirtualData(this)
            this.timeTable.acq = tic;
            fus              = this.measVars.algo.usSignal.fSin;
            samplesPerPulse  = this.measVars.algo.samples.samplesPerPulse;
            samplesPerTrain  = this.measVars.algo.samples.samplesPerTrain;
            numOfTrains      = this.measVars.algo.samples.numOfTrains;
            samplesPerSignal = this.measVars.algo.samples.samplesPerSignal;
            samplesPerAcq    = this.measVars.algo.digitizer.samplesPerAcq;
            
            tPulse   = this.measVars.algo.timing.tUS(1:samplesPerPulse);
            sig      = sin(2*pi*fus*tPulse);
            sigTrain = zeros(1, samplesPerTrain);
            
            sigTrain(1:samplesPerPulse) = sig;
            
            sig1 = repmat(sigTrain, 1, numOfTrains);
            sig2 = repmat(circshift(sigTrain,   samplesPerPulse), 1, numOfTrains);
            sig3 = repmat(circshift(sigTrain, 2*samplesPerPulse), 1, numOfTrains);
            sig4 = repmat(circshift(sigTrain, 3*samplesPerPulse), 1, numOfTrains);
            
            prePhantom  = this.measVars.algo.samples.prePhantomSamples;
            preTrigger  = this.measVars.algo.digitizer.preTriggerSamples;
            postSamples = samplesPerAcq - (prePhantom + preTrigger + samplesPerSignal);
            
            preSig  = zeros(1, prePhantom+preTrigger);
            postSig = zeros(1, postSamples);
            
            totSig1 = [preSig, sig1, postSig];
            totSig2 = [preSig, sig2, postSig];
            totSig3 = [preSig, sig3, postSig];
            totSig4 = [preSig, sig4, postSig];

            bufferDataOut = [totSig1; totSig2; totSig3; totSig4];
            noise         = 2*rand(4, samplesPerAcq);
            
            bufferDataOut = bufferDataOut + noise;
            this.algo.setRawData(bufferDataOut);
            this.timeTable.acq = toc(this.timeTable.acq);
            this.timeTable.moveData = tic;
        end
        
        function setSamplingTime(this, timeToSample)
            this.algo.uVars.timeToSample = timeToSample;
            this.vars.algo = this.algo.updateAlgoUserVars(this.vars.algo.uVars);
        	this.configDigitizer();
        end
     
        function timeTable = getTimeTable(this)
           timeTable = this.timeTable; 
        end
        
        function timeTable = getInnerTimeTables(this)
            timeTable.algo      = this.algo.getTimeTable();
            timeTable.digitizer = this.digitizer.getTimeTable();
        end
        
        function updateTimeTable(this)
            inTimeTable              = this.getInnerTimeTables();
            this.timeTable.algo      = inTimeTable.algo;
            this.timeTable.digitizer = inTimeTable.digitizer;
            
            if this.owned
                this.owner.displayAOTimeTable();
            end
        end
        
        function result = getResult(this)
            result = this.result; 
        end
        
        %Graphics
        function plotAll(this)
            if this.measVars.figs.validStruct.fullSignal
                if ~isfield(this.result, 'rawData')
                    fprintf("AOI: You asked to display Raw Data, but rawData isn't available in reults.\n");
                else
                    this.graphics.plotSignal('fullSignal');
                end
            end
            
            if this.measVars.figs.validStruct.measSamples
                if ~isfield(this.result, 'rawData')
                    fprintf("AOI: You asked to display Measured Samples, but rawData isn't available in reults.\n");
                else
                    this.graphics.plotSignal('measSamples');
                end
            end

            if this.measVars.figs.validStruct.netSignal
                if ~isfield(this.result, 'netSignal')
                    fprintf("AOI: You asked to display Net Signal, but netSignal isn't available in reults.\n");
                else
                    this.graphics.plotSignal('netSignal');
                end
            end
            
            if this.measVars.figs.validStruct.deMul
                if ~isfield(this.result, 'deMultiplexed')
                    fprintf("AOI: You asked to display Demultiplexed Signal, but deMul isn't available in reults.\n");
                else
                    this.graphics.plotSignal('deMul');
                end
            end
            
            if this.measVars.figs.validStruct.reshapedSignal
                if ~isfield(this.result, 'reshapedSignal')
                    fprintf("AOI: You asked to display Reshaped Singal, but reshapedSignal isn't available in reults.\n");
                else
                    this.graphics.plotReshapedSignal();
                end
            end  
      
            if this.measVars.figs.validStruct.qAvgChFFT
                 if ~isfield(this.result, 'qAvgChFFT')
                    fprintf("AOI: You asked to display quant-averaged FFT, but qAvgFFT isn't available in reults.\n");
                 else
                    this.graphics.plotFFT('qAvgChFFT');
                 end
            end

            if this.measVars.figs.validStruct.unFittedFFT
                 if ~isfield(this.result, 'unFittedFFTShift')
                    fprintf("AOI: You asked to display unfitted FFT, but unFittedFFTShift isn't available in reults.\n");
                 else
                    this.graphics.plotFFT('unFittedFFT');
                 end
            end
            
            if this.measVars.figs.validStruct.fittedFFT
                 if ~isfield(this.result, 'fittedFFT')
                    fprintf("AOI: You asked to display fitted FFT, but fittedData isn't available in reults.\n");
                 else
                    this.graphics.plotFFT('fittedFFT');
                 end
            end
            
            if this.measVars.figs.validStruct.phi
                if ~isfield(this.result, 'phi')
                    fprintf("AOI: You asked to display Phi, but phi isn't available in reults.\n");
                else
                    this.graphics.plotPhi();
                end
            end  
        end
        
        function gH = getGraphicsHandle(this)
            gH = this.graphics;
        end
        
        function setFigsVars(this, figs)
            this.graphics.setUserVars(figs);
            this.graphics.setGraphicsDynamicVars()
            this.plotAll();
        end
        
        % Saved data Operation
        
        % The following 3 functions should be called consecutively in order
        % to load data and analyse/plot it.
        function vars = loadVarsToAO(this, arg, varargin)
            % Importing vars and sets them into the object but not fully
            % configuring the object.
            vars = [];
            if isstring(arg) || ischar(arg)
                %argument is file path
                vars = load(arg);
            elseif isstruct(arg) 
                %argument is struct
                vars = arg;
            else
                fprintf("AOI: vars type not supported. Input argument should be filepath or variables struct.\n");
                return;
            end
            
            % this functionality allows to load variables without faulting
            % on graphics. if user mentions 'noFigs' or simply suppliend 
            % only a filename as input, graphics will be automatically
            % disabled. If user is interested in plotting results later on,
            % he may supply the AOGraphics uVars struct in this point and
            % spare him/herself the need in calling another 
            % ao.setMeasVars() later on.
            if isempty(varargin) || strcmp(varargin{1}, 'noFigs')
                vars.figs = AOGraphics.createGraphicsUserVars();
            elseif strcmp(varargin{2}, 'createFigs')
                if length(varargin)<2 || isempty(varargin{2})
                    fprintf("AOI: You asked to create load vars with figures suppoet but didn't supplied figs struct.\n");
                    return; 
                end
                vars.figs = varargin{2};
            end
                        
            if isfield(vars, 'uVars')
                uVarsLoaded = vars.uVars;
                % Deactivate file system - all flags are false by default
                % this is to make sure analysis will not overrun the same file
                %TODO: apply functionality similar to graphics.
                uVarsLoaded.fileSystem = fileSystemAO.uVarsCreate(); 
                this.setMeasVars(uVarsLoaded); %algo vars are calculated in here.
                fprintf("AOI: New rawData was loaded to AOI.\n");
            else
                fprintf("AOI: Error while loading variables: Your file has no uVars field.\n");
            end
        end
        
        function rawData = loadRawDataToAO(this, arg)
            rawData = [];
            if isstring(arg) || ischar(arg)
                %argument is file path
                results = load(arg);
            elseif isstruct(arg)
                %argument is struct
                results = arg;
            else
                fprintf("AOI: rawData type not supported, input argument should be filepath or variables struct/\n");
                return;
            end
               
            if isfield(results, 'rawData')
                this.rawData = results.rawData;
                rawData = results.rawData;
                fprintf("AOI: New rawData was loaded to AOI.\n");
            else
                fprintf("AOI: Error while loading raw data: Your file has no rawData field\n");
            end
        end
        
        function res = analyseLoadedData(this)
            %this function assums your data was already loaded and is
            %available in this.rawData.
            
            % Setup graphics
            this.graphics.setGraphicsDynamicVars();
            this.fileSystem.configFileSystem();
            this.fileSystem.saveVarsToDisk();

            % Analyse
            fprintf ("AOI: Analyzing!\n")
            this.algo.setRawData(this.rawData);
            res                    = this.algo.analyse();
            this.result            = res;
            this.result.rawData    = this.rawData;
            
            % Plot reseults
            this.plotAll();
            
            % Save results
            fprintf("AOI: Saving Results To Disk\n");
            this.fileSystem.saveResultsToDisk(this.result) 
            
            fprintf ("AOI: Done AO\n")
        end
        
        function saveVarsToDisk(this, path)
            % path is relative to projPath
            this.fileSystem.saveVarsToDisk(this.getAOVars(), path);
        end
        
    end
end