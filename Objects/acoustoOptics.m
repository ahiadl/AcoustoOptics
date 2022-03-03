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
            aoVars.cycPerPulse       = 1;      %[#] 
            aoVars.fSin              = 1.25e6; %[Hz]              
            aoVars.fSqnc             = 10e3; %[Hz]
            aoVars.frameTime         = 0.002; %[s]

            % Sampling Clk
            aoVars.fs                = 5e6;    %[Hz]
            aoVars.sClkDcyc          = 50;     %[%]
            aoVars.fgClk             = 100e6;  %[S/s]

            % Digitizer
            aoVars.timeToSample      = 0.1; %[s]
            aoVars.channels          = 2; %[#]

            % Frequency
            aoVars.envDC             = 100e3;
            aoVars.envUS             = 78e3;

            % Geometry & Length
            aoVars.c                 = 1500;   %[m/s]
            aoVars.distFromPhantom   = 0; %[m]

            % General Operation
            aoVars.useFrame            = true;
            aoVars.useGPU              = true;
            aoVars.useHadamard         = true;
            aoVars.contHadamard        = false;
            aoVars.highResAO           = false;
            aoVars.analyzeSingleCh     = true;
            aoVars.contSpeckleAnalysis = false;

            aoVars.useVirtualData   = false;
            aoVars.limitByN         = true;
            aoVars.N                = 10;
            aoVars.dispTimeTable    = true; % relevant only with GUI
            aoVars.skipParamsCheck  = false;
            aoVars.extCropDef       = false;
            aoVars.extCropSamples   = 1000;
            
            % IO
            aoVars.IOPort = 1;
            aoVars.IOLine = 4;

            % Export
            aoVars.exportData.meas        = false;
            aoVars.exportData.signal      = false;
            aoVars.exportData.deMul       = false;
            aoVars.exportData.reshaped    = false;
            aoVars.exportData.fft         = false;
            aoVars.exportData.usCompCmplx = false;
        end
        
        function uVars  = createUserVars()
            % Algo and Peripherals
            uVars.ao         = acoustoOptics.aoVarsCreate();
            % Figures
            uVars.figs       = AOGraphics.createUserVars();
            % File system     
            uVars.fileSystem = fileSystemAO.uVarsCreate();
        end

        function newExtVars = extVarsCreate()
            newExtVars.algo      = Algo.createUserVars();
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
                             (strcmp(names{j}, 'fSin')       || strcmp(names{j}, 'fSqnc') ||...
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
            
            this.uVars   = this.createUserVars();
            this.extVars = this.extVarsCreate();
            
%             this.limits.measLimit = 10;
            this.algo.setMeasLimit(10);
            
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
        function setVars(this, uVars)
%             user gives reduced variables. this function convert these 
%             variables to extended variables suitable for each
%             submodule\peripheral. this function extend the variables by
%             adding default acousto optics values to the user reduced
%             variables
            fprintf("AOI: ------- Sets Extended Vars From Reduced Vars ----------\n")
            this.uVars = uVars;

            %-----------------
            % Algo 
            %-----------------
            fprintf("AOI:  1. Set & Calculate Algorithm Vars\n")
            algoVars      = this.algo.createUserVars();

            % US Signal
            algoVars.cycPerPulse     = uVars.ao.cycPerPulse;
            algoVars.fSin            = uVars.ao.fSin;
            algoVars.fSqnc           = uVars.ao.fSqnc;
            algoVars.frameTime       = uVars.ao.frameTime;
            
            % Sampling Clk
            algoVars.fs              = uVars.ao.fs;
            algoVars.sClkDcyc        = uVars.ao.sClkDcyc;
            algoVars.fgClk           = uVars.ao.fgClk;
            
            % Geometry
            algoVars.c               = uVars.ao.c;
            algoVars.distFromPhantom = uVars.ao.distFromPhantom;
            
            % Digitizer
            algoVars.timeToSample       = uVars.ao.timeToSample;
            algoVars.channels           = uVars.ao.channels;

            % Frequency
            algoVars.envDC           = uVars.ao.envDC;
            algoVars.envUS           = uVars.ao.envUS;
            
            % General Operation
            algoVars.useFrame            = uVars.ao.useFrame;
            algoVars.useGPU              = uVars.ao.useGPU;
            algoVars.useHadamard         = uVars.ao.useHadamard;
            algoVars.contHadamard        = uVars.ao.contHadamard;
            algoVars.highResAO           = uVars.ao.highResAO;
            algoVars.analyzeSingleCh     = uVars.ao.analyzeSingleCh;
            algoVars.contSpeckleAnalysis = uVars.ao.contSpeckleAnalysis;
            
            algoVars.cutArtfct           = uVars.ao.cutArtfct;
            algoVars.artfctIdxVec        = uVars.ao.artfctIdxVec;
            
            % Export
            algoVars.export.meas           = uVars.ao.exportData.meas        || uVars.figs.validStruct.meas     || uVars.fileSystem.saveMeas;
            algoVars.export.signal         = uVars.ao.exportData.signal      || uVars.figs.validStruct.signal   || uVars.fileSystem.saveSignal;
            algoVars.export.deMul          = (uVars.ao.exportData.deMul      || uVars.figs.validStruct.deMul    || uVars.fileSystem.saveDeMul) && algoVars.useHadamard;
            algoVars.export.reshaped       = uVars.ao.exportData.reshaped    || uVars.figs.validStruct.reshaped || uVars.fileSystem.saveReshaped;
            algoVars.export.fft            = uVars.ao.exportData.fft         || uVars.fileSystem.saveFFT;
            algoVars.export.usCompCmplx    = uVars.ao.exportData.usCompCmplx || uVars.fileSystem.saveFFT;

            this.measVars.algo = this.algo.setVars(algoVars);
            AO.splitMeas = this.measVars.algo.general.splitMeas;
            AO.splitNum  = this.measVars.algo.general.splitNum;
            
            %-----------------
            %fGen
            %-----------------
            fprintf("AOI:  2. Set Function Generator Vars\n")
            fGenVars      = this.fGen.uVarsCreate();

            fGenVars.ch{1}.daqFilter     = 50;
            fGenVars.ch{1}.amp           = 2;
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

            fGenVars.ch{1}.Sclk = algoVars.fgClk;
            fGenVars.ch{2}.Sclk = algoVars.fgClk; 

            %-----------------
            % Digitizer
            %-----------------
            fprintf("AOI:  3. Sets Digitizer Vars\n")
            digiVars = this.digitizer.uVarsCreate();

            digiVars.mode     = 'TS'; 
            digiVars.fs       = algoVars.fs;
            digiVars.useGPU   = algoVars.useGPU;
            digiVars.channels = algoVars.channels;           
            
            digiVars.triggerDelay = 0;
            digiVars.extClk  = true;
            digiVars.extTrig = true;
            digiVars.draw    = false;
            
            digiVars.bufferSizeBytes = 4*(2^20);
            digiVars.samplesPerMeas  = this.measVars.algo.samples.samplesPerMeas;
            
            digiVars.exportCropped = uVars.figs.validStruct.cropped && ~AO.splitMeas;
            if this.uVars.ao.extCropDef
                digiVars.croppedSamples = this.uVars.ao.extCropSamples;
            else
                digiVars.croppedSamples = 2*this.measVars.algo.samples.samplesPerSqnc;
            end
            %-----------------
            % IO
            %-----------------
            IOVars = this.IO.uVarsCreate;

            IOVars.mode = 0; % Output Only
            IOVars.port = uVars.ao.IOPort;
            IOVars.line = uVars.ao.IOLine;

            %-----------------
            % FileSystem Vars
            %-----------------
            fprintf("AOI:  4. Extracting File system variables.\n");
            if AO.splitMeas
                fileSystemVars.saveSignal   = false;
                fileSystemVars.saveDeMul    = false;
                fileSystemVars.saveReshaped = false;
                fileSystemVars.saveFFT      = false;
            else
                fileSystemVars.saveSignal   = uVars.fileSystem.saveSignal;
                fileSystemVars.saveDeMul    = uVars.fileSystem.saveDeMul && uVars.useHadamard;
                fileSystemVars.saveReshaped = uVars.fileSystem.saveReshaped;
                fileSystemVars.saveFFT      = uVars.fileSystem.saveFFT;
            end
            fileSystemVars.saveMeas = uVars.fileSystem.saveMeas;
            
            fileSystemVars.splitMeas = AO.splitMeas;
            fileSystemVars.splitNum  = AO.splitNum;
            fileSystemVars.useFrame  = algoVars.useFrame;
            fileSystemVars.chToAnalyze = this.measVars.algo.general.analysisReps;
            
            fileSystemVars.saveResults  = uVars.fileSystem.saveResults;
            fileSystemVars.saveFigs     = uVars.fileSystem.saveFigs;
            fileSystemVars.dontSaveVars = uVars.fileSystem.dontSaveVars;
            
            fileSystemVars.dirPath    = uVars.fileSystem.dirPath;
            fileSystemVars.projName   = uVars.fileSystem.projName;
            fileSystemVars.resDirName = uVars.fileSystem.resDirName;

            fileSystemVars.extProject = uVars.fileSystem.extProject;

            this.fileSystem.setUserVars(fileSystemVars);

            %-----------------
            % Graphics
            %-----------------
            fprintf("AOI:  5. Sets Graphic Variables.\n")

            figsVars = this.graphics.createGraphicsVars();
            
            if uVars.figs.depthIdx > this.measVars.algo.samples.numOfPos
                figsVars.depthIdx = 1;
            else
                figsVars.depthIdx = uVars.figs.depthIdx;
            end
            
            if uVars.figs.ch > this.measVars.algo.general.channels
                figsVars.ch = 1;
            else
                figsVars.ch = uVars.figs.ch;
            end
            
            if uVars.figs.frame > this.measVars.algo.samples.framesPerSignal
                figsVars.frame = 1;
            else
                figsVars.frame = uVars.figs.frame;
            end
            figsVars.reopenFigures   = uVars.figs.reopenFigures;
            figsVars.displayFullFFT  = uVars.figs.displayFullFFT;
            figsVars.FFTenv          = uVars.figs.FFTenv;
            figsVars.intExt          = uVars.figs.intExt;
            
            figsVars.validStruct.usSignal = uVars.figs.validStruct.usSignal;
            figsVars.validStruct.extClk   = uVars.figs.validStruct.extClk;
            
            figsVars.validStruct.cropped  = uVars.figs.validStruct.cropped  && ~AO.splitMeas;
            figsVars.validStruct.meas     = uVars.figs.validStruct.meas     && this.measVars.algo.general.export.meas;
            figsVars.validStruct.signal   = uVars.figs.validStruct.signal   && this.measVars.algo.general.export.signal;
            figsVars.validStruct.deMul    = uVars.figs.validStruct.deMul    && this.measVars.algo.general.export.deMul && algoVars.useHadamard;
            figsVars.validStruct.reshaped = uVars.figs.validStruct.reshaped && this.measVars.algo.general.export.reshaped;
            
            figsVars.validStruct.unFittedFFTAllCh    = uVars.figs.validStruct.unFittedFFTAllCh;
            figsVars.validStruct.unFittedFFTSingleCh = uVars.figs.validStruct.unFittedFFTSingleCh;
            figsVars.validStruct.fittedFFTSingleCh   = uVars.figs.validStruct.fittedFFTSingleCh;
            figsVars.validStruct.fittedFFTAllCh      = uVars.figs.validStruct.fittedFFTAllCh;
            figsVars.validStruct.avgFFT              = uVars.figs.validStruct.avgFFT;
            figsVars.validStruct.normFFT             = uVars.figs.validStruct.normFFT;
            
            figsVars.validStruct.phi    = uVars.figs.validStruct.phi;
            figsVars.validStruct.rawPhi = uVars.figs.validStruct.rawPhi;
            figsVars.validStruct.phiLog = uVars.figs.validStruct.phiLog;
            
            figsVars.depthAxType     = uVars.figs.depthAxType;
            
            if ~uVars.ao.useFrame
                figsVars.validStruct.unFittedFFT = false;
            end
            
            figsVars.extH = uVars.figs.extH;
            
            figsVars.numOfChannels = uVars.ao.channels;
            figsVars.df            = this.measVars.algo.freq.df;
            figsVars.fBar          = gather(this.measVars.algo.freq.fBar);
            figsVars.fIdx          = this.measVars.algo.freq.fUSIdx;
            figsVars.depthVecNorm  = this.measVars.algo.len.depthNorm;
            figsVars.depthVecCntr  = this.measVars.algo.len.depthCntr;
            figsVars.depthVecZero  = this.measVars.algo.len.depthZero;
            figsVars.depthVecIdx   = this.measVars.algo.len.depthIdx;

            figsVars.depthVecNormRaw = this.measVars.algo.len.depthVecAlgo;
            
            figsVars.tVecUS      = this.measVars.algo.timing.tVecUS;
            figsVars.tVecSampleClk  = this.measVars.algo.timing.tVecSampleClk;
            
            if ~AO.splitMeas
                figsVars.tVecMeas = this.measVars.algo.timing.tVecMeas;
                figsVars.tVecSig  = this.measVars.algo.timing.tVecSig;
                figsVars.tVecPos  = this.measVars.algo.timing.tVecPosPerFrame; 
            else
%                 figsVars.validStruct.measSamples    = false;
%                 figsVars.validStruct.netSignal      = false;
%                 figsVars.validStruct.deMul          = false;
%                 figsVars.validStruct.reshaped = false;
                
                figsVars.tVecMeas = [];
                figsVars.tVecSig  = [];
                figsVars.tVecPos  = [];
            end
            
            figsVars.analyzeSingleCh = uVars.ao.analyzeSingleCh;
            figsVars.singleChIdx     = uVars.figs.singleChIdx;
            
            this.graphics.setUserVars(figsVars);
            this.measVars.figs = figsVars;

            %-----------------
            % AO
            %-----------------
            
            % Writing Extended Vars to Acousto Optics Object
            fprintf("AOI:  8. Set AO Object Variables.\n");
            % Set AO vars
            AO.useVirtualData   = uVars.ao.useVirtualData;
            AO.limitByN         = uVars.ao.limitByN;
            AO.N                = uVars.ao.N;
            AO.dispTimeTable    = uVars.ao.dispTimeTable;
            AO.skipParamsCheck  = uVars.ao.skipParamsCheck;
            
            newExtVars.AO         = AO;
            newExtVars.algo       = algoVars;
            newExtVars.fGen       = fGenVars;
            newExtVars.digitizer  = digiVars;
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
            
            % We put an empty graphic request in order to shrink the size
            % of the vars mat file.
            aoVars.measVars.figs = AOGraphics.createGraphicsVars();
            aoVars.extVars.figs  = AOGraphics.createOwnerVars();
%             aoVars.uVars.figs    = AOGraphics.createUserVars();
            aoVars.uVars.figs.extH = [];

            aoVars.measVars.algo.timing.tSigVec  = [];
            aoVars.measVars.algo.timing.tMeasVec = [];
            aoVars.measVars.algo.timing.tAcqVec  = [];
        end
        
        function setMeasLimit(this, time)
            this.algo.setMeasLimit(time)            
        end
        
        function [res] = getResultStruct(this)
            [~, resStruct, ~, resSplit] = this.algo.getResultsStruct();
            
            if this.measVars.AO.splitMeas
                res = resSplit;
            else
                res = resStruct;
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
            sigData = this.measVars.algo.usSignal.sigData;
            clkData = this.measVars.algo.sClk.clkData;
            
            dataCh1 = this.fGen.getChDataStruct();
            dataCh2 = this.fGen.getChDataStruct();
            
            dataCh1.data = sigData;
            dataCh2.data = clkData;
            dataCh1.dataLen = length(dataCh1.data);
            dataCh2.dataLen = length(dataCh2.data);

            this.measVars.fGen = this.extVars.fGen; 
            this.measVars.fGen.sigData = sigData;
            this.measVars.fGen.clkData = clkData;
            
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
            this.digitizer.setVars(this.extVars.digitizer);
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
            [~, ~, dataInSplit, ~] = this.algo.getResultsStruct();
            % Save AO variables (if needed, decided by fileSystem)
            this.fileSystem.saveVarsToDisk();            
            % Measure data and analyse it with AO Algorithm
            if this.measVars.AO.splitMeas
                for i = 1:this.measVars.AO.splitNum
                    fprintf("AOI: Split %d/%d: ", i, this.measVars.AO.splitNum);
                    dataInSplit(:,i) = this.measureAndAnalyse();
                end
                this.result = this.algo.analyseSplittedData(dataInSplit);
            else 
                this.result = this.measureAndAnalyse();
            end
            
            % Set data to graphics object
            this.timeTable.setResultsToGraphics = tic;
            this.graphics.setData(gather(this.result));
            this.timeTable.setResultsToGraphics = toc(this.timeTable.setResultsToGraphics);
            
            % Plot Results (if needed, decided by AO)
            this.timeTable.plotAll = tic;
            this.plotAll();
            this.timeTable.plotAll = toc(this.timeTable.plotAll);
            
            % Save results (if needed, decided by fileSystem)
            this.timeTable.saveData = tic;
            this.fileSystem.saveResultsToDisk(this.result)
            this.fileSystem.closeFileSystem();
            this.timeTable.saveData = toc(this.timeTable.saveData);
            
            % Collect time statistics
            if this.measVars.AO.dispTimeTable
                this.updateTimeTable();
            end
            res = this.result;
        end
        
        function res = liveAcoustoOptics(this) 
            if this.measVars.AO.limitByN
                this.fileSystem.initLiveAOFS(this.getVars());
                this.runLive = true;
 
                for i=1:this.measVars.AO.N
                    fprintf("AOI: Live AO: %d/%d.\n", i, this.measVars.AO.N);
                    this.fileSystem.setLiveAOInd(i);
                    if this.measVars.algo.general.useFrame
                        res(i) = this.runAcoustoOptics();
                    else
                        res = this.runAcoustoOptics();
                    end
                    pause(1)
                    if ~this.runLive
                        fprintf("AOI: Live AO was STOPPED.\n");
                        break;
                    end
                end
                
                this.fileSystem.closeLiveAoFileSystem();
                this.result = res;
            else %infinite number of measurements
                %saving is not allowed in unlimited live AO
                this.runLive = true;
                while(this.runLive)
                   this.runAcoustoOptics();
                   pause(1);
                end
                fprintf("AOI: Live AO was STOPPED.\n");
                res = this.result;
            end           
        end
        
        function res = measureAndAnalyse(this)
           if this.measVars.AO.useVirtualData
               % Generate virtual data (DEBUG)
                fprintf ("AOI: Generating synthetic data...")
                bufferDataOut = this.createVirtualData();
           elseif this.periAvail.completeHardwareAvail
                % Get data from digitizer (Measure)
                fprintf ("AOI: Acquiring...")
                this.timeTable.acq = tic;
                this.IO.open();
                bufferDataOut      = this.digitizer.acquire();
                this.IO.close();
                this.timeTable.acq = toc(this.timeTable.acq);
                this.displayCroppedData();
            else
                 % No source of data (Error)
                 fprintf("AOI: No source of data found.\n");  
                 res = [];
                 return;
            end
            
            this.timeTable.moveData = tic;
            this.algo.setRawData(bufferDataOut);
            bufferDataOut           = gather(bufferDataOut);
            this.timeTable.moveData = toc(this.timeTable.moveData);

            if this.measVars.AO.splitMeas
                this.fileSystem.saveMeasToDisk(bufferDataOut)
                bufferDataOut = [];
            end
            
            fprintf ("Analyzing...")
            this.timeTable.analyse = tic;
            res                    = this.algo.analyse();

            fprintf ("Done!\n")
            this.timeTable.analyse = toc(this.timeTable.analyse);
        end
        
        % Misc
        function bufferDataOut = createVirtualData(this)
            this.timeTable.acq = tic;
            fus              = this.measVars.algo.usSignal.fSin;
            samplesPerPulse  = this.measVars.algo.samples.samplesPerPulse;
            samplesPerSqnc   = this.measVars.algo.samples.samplesPerSqnc;
            sqncPerSignal    = this.measVars.algo.samples.sqncPerSignal;
            samplesPerSignal = this.measVars.algo.samples.samplesPerSignal;
            samplesPerAcq    = this.measVars.algo.digitizer.samplesPerAcq;
            
            tPulse   = this.measVars.algo.timing.tSigVec(1:samplesPerPulse);
            sig      = sin(2*pi*fus*tPulse);
            sigSqnc = zeros(1, samplesPerSqnc);
            
            sigSqnc(1:samplesPerPulse) = sig;
            
            prePhantom  = this.measVars.algo.samples.prePhantomSamples;
            preTrigger  = this.measVars.algo.digitizer.preTriggerSamples;
            postSamples = samplesPerAcq - (prePhantom + preTrigger + samplesPerSignal);
            
            preSig  = zeros(1, prePhantom+preTrigger);
            postSig = zeros(1, postSamples);
            
            bufferDataOut = zeros(this.measVars.digitizer.vars.channels, samplesPerAcq);
            
            for i=1:this.measVars.digitizer.vars.channels
                tmpSig = repmat(circshift(sigSqnc,   (i-1)*samplesPerPulse), 1, sqncPerSignal);
                bufferDataOut(i,:) = [preSig, tmpSig, postSig]; %+  2*rand(1, samplesPerAcq);
            end

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
            this.graphics.plotSignal('meas');
            this.graphics.plotSignal('signal');
            this.graphics.plotSignal('deMul');
            this.graphics.plotReshapedSignal();
            this.graphics.plotFFT('unFittedFFTAllCh');
            this.graphics.plotFFT('unFittedFFTSingleCh')
            this.graphics.plotFFT('fittedFFTSingleCh');
            this.graphics.plotFFT('fittedFFTAllCh');
            this.graphics.plotFFT('avgFFT');
            this.graphics.plotFFT('normFFT');
            this.graphics.plotPhi('phi');
            this.graphics.plotPhi('rawPhi');
            this.graphics.plotPhi('phiLog');
        end
        
        function gH = getGraphicsHandle(this)
            gH = this.graphics;
        end
        
        function setFigsVars(this, figs)
            this.graphics.setUserVars(figs);
            this.graphics.setGraphicsDynamicVars()
            this.graphics.setSingleChIdx(figs.singleChIdx)
            this.plotAll();
        end
        
        function displayCroppedData(this)
            this.graphics.setCroppedData(this.digitizer.getCroppedData)
            this.graphics.plotSignal('cropped');
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