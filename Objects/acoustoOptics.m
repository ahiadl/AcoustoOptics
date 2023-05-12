classdef acoustoOptics < handle
    %CONSIST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ownerGUI;
        
        % Objects:
        fGen;
        digitizer;
        IO;
        algo;
        fileSystem;
        graphics;

        % Data:
        resultNew;
        rawData;
        result;
        splits;

        % Variables:
        uVars;    % Reduced from user to acoustoOptics
        extVars;  % Extended from acoustoOptics to submodules
        measVars; % Measurement vars (after calculation) given from submodules to acoustoOptics
        oldExtVars;
        
        % Control Vars:
        changeLog;
        connected;
        periAvail;
        stopMeas;
        pauseMeas;
        contLive;
        graphicsNames;
        ownedByGUI;

        % Run Statistics
        timeTable;
    end
    
    methods (Static) 
        function aoVars = aoVarsCreate()
            % Algo
            aoVars.cycPerPulse       = 1;      %[#] 
            aoVars.fSin              = 1.25e6; %[Hz]              
            aoVars.fSqnc             = 10e3; %[Hz]
            aoVars.frameTime         = 0.002; %[s]
            
            % Sampling Clk
            aoVars.fs       = 5e6;    %[Hz]
            aoVars.sClkDcyc = 50;     %[%]
            aoVars.fgClk    = 100e6;  %[S/s]
 
            % Digitizer
            aoVars.timeToSample   = 0.1; %[s]
            aoVars.channels       = 2; %[#]
            aoVars.extCropDef     = false;
            aoVars.extCropSamples = 1000;
            aoVars.measTimeLimit   = 1;

            % Frequency
            aoVars.envDC = 100e3;
            aoVars.envUS = 50e3;
 
            % Geometry & Length
            aoVars.c               = 1500;   %[m/s]
            aoVars.distFromPhantom = 0; %[m]

            % General Operation
            aoVars.useFrame        = true;
            aoVars.useGPU          = true;
            aoVars.useHadamard     = true;
            aoVars.analyzeSingleCh = true;
            aoVars.acCoupling      = true;
            aoVars.keepSplits      = false;

            aoVars.displayBuildUp  = false;
            aoVars.displayBUEvery  = 10;

            aoVars.useCalibration  = false;           
            aoVars.timeToSampleCal = 2;
            
            aoVars.cutArtfct    = false;
            aoVars.artfctIdxVec = [];
            
            aoVars.calcMuEff  = false;
            aoVars.muEffModel = 'Uniform';
            aoVars.muEffIdxs  = [];
            
            aoVars.extSNRDef   = false;
            aoVars.extNoiseIdx = [];
            aoVars.extPeakIdx  = [];

            aoVars.noiseStd = [];

            aoVars.useVirtualData      = false;
            aoVars.virtualDataNoiseSTD = 0;
            aoVars.limitByN            = true;
            aoVars.N                   = 10;
            aoVars.dispTimeTable       = true; % relevant only with GUI

            aoVars.skipParamsCheck = false;

            aoVars.uploadToTelegram = false;
            aoVars.telegramChatID   = '-512325870';            
            
            % fGen
            aoVars.usPower = 100;
            
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

            this.graphicsNames = this.graphics.getGraphicsNames();
            
            this.connected = false;
            this.stopMeas  = false;
            this.pauseMeas = false;
            this.contLive  = false;
        end 

        function init(this, ownerGUI)
            fprintf("AOI: ------- Initiating AcoustoOptic ----------\n");
            if this.connected   
                fprintf("AOI: Acousto Optics system is already connected, no need to reconnect.\n");
            else
                fprintf("AOI: 1. Resetting Instruments\n")
                instrreset; %#ok<INSTRR> 
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
                    this.ownerGUI = ownerGUI;
                    this.ownedByGUI = true;
                end
            end
            this.periAvail.completeHardwareAvail =  this.periAvail.digitizer && this.periAvail.fGen && this.periAvail.IO;
        end
        
        % Complete runs
        function res = AOI(this, uVars)
            this.setVars(uVars);
            this.configurePeripherals();
            res = this.runAcoustoOptics();
        end
        
        function res = calibrate(this, uVars)
            this.setVars(uVars);
            this.configPeripherals();
            res = this.runCalibration();
        end
        
        % Variables
        function algoVars = setAlgoVars(this)
            uVarsAO = this.uVars;
            
            %-----------------
            % Algo 
            %-----------------
            fprintf("AOI:  1. Set & Calculate Algorithm Vars\n")
            algoVars      = this.algo.createUserVars();

            % US Signal
            algoVars.cycPerPulse     = uVarsAO.ao.cycPerPulse;
            algoVars.fSin            = uVarsAO.ao.fSin;
            algoVars.fSqnc           = uVarsAO.ao.fSqnc;
            algoVars.frameTime       = uVarsAO.ao.frameTime;
            
            % Sampling Clk
            algoVars.fs              = uVarsAO.ao.fs;
            algoVars.sClkDcyc        = uVarsAO.ao.sClkDcyc;
            algoVars.fgClk           = uVarsAO.ao.fgClk;
            
            % Geometry
            algoVars.c               = uVarsAO.ao.c;
            algoVars.distFromPhantom = uVarsAO.ao.distFromPhantom;
            
            % Digitizer
            algoVars.timeToSample       = uVarsAO.ao.timeToSample;
            algoVars.channels           = uVarsAO.ao.channels;

            % Frequency
            algoVars.envDC           = uVarsAO.ao.envDC;
            algoVars.envUS           = uVarsAO.ao.envUS;
            
            % General Operation
            algoVars.useFrame            = uVarsAO.ao.useFrame;
            algoVars.useGPU              = uVarsAO.ao.useGPU;
            algoVars.useHadamard         = uVarsAO.ao.useHadamard;
            algoVars.analyzeSingleCh     = uVarsAO.ao.analyzeSingleCh;
            algoVars.calibrate           = uVarsAO.ao.useCalibration;
            algoVars.acCoupling          = uVarsAO.ao.acCoupling;
            algoVars.measTimeLimit       = uVarsAO.ao.measTimeLimit;

            algoVars.calcMuEff           = uVarsAO.ao.calcMuEff;
            algoVars.muEffModel          = uVarsAO.ao.muEffModel;
            algoVars.muEffIdxs           = uVarsAO.ao.muEffIdxs;
            
            algoVars.extSNRDef   = uVarsAO.ao.extSNRDef;  
            algoVars.extNoiseIdx = uVarsAO.ao.extNoiseIdx;
            algoVars.extPeakIdx  = uVarsAO.ao.extPeakIdx ;

            algoVars.noiseStd            = uVarsAO.ao.noiseStd;

            algoVars.cutArtfct           = uVarsAO.ao.cutArtfct;
            algoVars.artfctIdxVec        = uVarsAO.ao.artfctIdxVec;
            
            % Export
            algoVars.export.meas           = uVarsAO.ao.exportData.meas        || uVarsAO.figs.validStruct.meas     || uVarsAO.fileSystem.saveMeas;
            algoVars.export.signal         = uVarsAO.ao.exportData.signal      || uVarsAO.figs.validStruct.signal   || uVarsAO.fileSystem.saveSignal;
            algoVars.export.deMul          = (uVarsAO.ao.exportData.deMul      || uVarsAO.figs.validStruct.deMul    || uVarsAO.fileSystem.saveDeMul) && algoVars.useHadamard;
            algoVars.export.reshaped       = uVarsAO.ao.exportData.reshaped    || uVarsAO.figs.validStruct.reshaped || uVarsAO.fileSystem.saveReshaped;
            algoVars.export.fft            = uVarsAO.ao.exportData.fft         || uVarsAO.fileSystem.saveFFT;
            algoVars.export.usCompCmplx    = uVarsAO.ao.exportData.usCompCmplx || uVarsAO.fileSystem.saveFFT;

            this.measVars.algo = this.algo.setVars(algoVars);
        end
        
        function fGenVars = setFGenVars(this)
            uVarsAO = this.uVars;
            %-----------------
            %fGen
            %-----------------
            fprintf("AOI:  2. Set Function Generator Vars\n")
            fGenVars      = this.fGen.uVarsCreate();

            fGenVars.ch{1}.daqFilter     = 50;
            fGenVars.ch{1}.amp           = uVarsAO.ao.usPower /100 *2;
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

            fGenVars.ch{1}.Sclk = this.measVars.algo.sClk.fgClk;
            fGenVars.ch{2}.Sclk = this.measVars.algo.sClk.fgClk; 
        end
        
        function digiVars = setDigitizerVars(this)
            uVarsAO = this.uVars;
            %-----------------
            % Digitizer
            %-----------------
            fprintf("AOI:  3. Sets Digitizer Vars\n")
            digiVars = this.digitizer.uVarsCreate();

            if this.measVars.algo.general.splitMeas
                digiVars.mode         = 'ContNPT';
                digiVars.timeToSample = this.measVars.algo.timing.frameT;
                digiVars.numAvg  = 1;
                digiVars.numMeas = this.measVars.algo.samples.framesPerSignal;
            else
                digiVars.mode     = 'TS';
                digiVars.samplesPerMeas = this.measVars.algo.samples.samplesPerMeas;
            end

            digiVars.fs       = this.measVars.algo.sClk.fs;
            digiVars.useGPU   = this.measVars.algo.general.useGPU;
            digiVars.channels = this.measVars.algo.general.channels;           
            
            digiVars.triggerDelay = 0;
            digiVars.extClk  = true;
            digiVars.extTrig = true;
            digiVars.draw    = false;
            
            digiVars.bufferSizeBytes = 16*(2^20);
            
            
            digiVars.exportCropped = uVarsAO.figs.validStruct.cropped && ~AO.splitMeas;
            if this.uVars.ao.extCropDef
                digiVars.croppedSamples = this.uVars.ao.extCropSamples;
            else
                digiVars.croppedSamples = 2*this.measVars.algo.samples.samplesPerSqnc;
            end
        end
        
        function IOVars = setIOVars(this)
            uVarsAO = this.uVars;
            %-----------------
            % IO
            %-----------------
            IOVars = this.IO.uVarsCreate;

            IOVars.mode = 0; % Output Only
            IOVars.port = uVarsAO.ao.IOPort;
            IOVars.line = uVarsAO.ao.IOLine;
        end
        
        function fileSystemVars = setFileSystemVars(this, uVarsFS)
            %-----------------
            % FileSystem Vars
            %-----------------
            fprintf("AOI:  4. Extracting File system variables.\n");
            if this.measVars.algo.general.splitMeas
                fileSystemVars.saveSignal   = false;
                fileSystemVars.saveDeMul    = false;
                fileSystemVars.saveReshaped = false;
                fileSystemVars.saveFFT      = false;
            else
                fileSystemVars.saveSignal   = uVarsFS.saveSignal;
                fileSystemVars.saveDeMul    = uVarsFS.saveDeMul && this.measVars.algo.general.useHadamard;
                fileSystemVars.saveReshaped = uVarsFS.saveReshaped;
                fileSystemVars.saveFFT      = uVarsFS.saveFFT;
            end

            fileSystemVars.saveMeas = uVarsFS.saveMeas;
            
            fileSystemVars.splitMeas   = this.measVars.algo.general.splitMeas;
            fileSystemVars.splitNum    = this.measVars.algo.general.splitNum;
            fileSystemVars.useFrame    = this.measVars.algo.general.useFrame;
            fileSystemVars.chToAnalyze = this.measVars.algo.general.analysisReps;
            
            fileSystemVars.saveResults  = uVarsFS.saveResults;
            fileSystemVars.saveFigs     = uVarsFS.saveFigs;
            fileSystemVars.dontSaveVars = uVarsFS.dontSaveVars;
            
            fileSystemVars.dirPath    = uVarsFS.dirPath;
            fileSystemVars.projName   = uVarsFS.projName;
            fileSystemVars.resDirName = uVarsFS.resDirName;

            fileSystemVars.extProject = uVarsFS.extProject;

            this.fileSystem.setUserVars(fileSystemVars);
        end
        
        function figsVars = setGraphicsVars(this, uVarsFigs)
            %-----------------
            % Graphics
            %-----------------
            fprintf("AOI:  5. Sets Graphic Variables.\n")

            figsVars = this.graphics.createGraphicsVars();
            
            if uVarsFigs.depthIdx > this.measVars.algo.samples.numOfPos
                figsVars.depthIdx = 1;
            else
                figsVars.depthIdx = uVarsFigs.depthIdx;
            end
            
            if uVarsFigs.ch > this.measVars.algo.general.channels
                figsVars.ch = 1;
            else
                figsVars.ch = uVarsFigs.ch;
            end
            
            if uVarsFigs.frame > this.measVars.algo.samples.framesPerSignal
                figsVars.frame = 1;
            else
                figsVars.frame = uVarsFigs.frame;
            end
            figsVars.reopenFigures   = uVarsFigs.reopenFigures;
            figsVars.displayFullFFT  = uVarsFigs.displayFullFFT;
            figsVars.FFTenv          = uVarsFigs.FFTenv;
            figsVars.intExt          = uVarsFigs.intExt;
            
            figsVars.validStruct.usSignal = uVarsFigs.validStruct.usSignal;
            figsVars.validStruct.extClk   = uVarsFigs.validStruct.extClk;
            
            figsVars.validStruct.cropped  = uVarsFigs.validStruct.cropped  && ~this.measVars.algo.general.splitMeas;
            figsVars.validStruct.meas     = uVarsFigs.validStruct.meas     &&  this.measVars.algo.general.export.meas;
            figsVars.validStruct.signal   = uVarsFigs.validStruct.signal   &&  this.measVars.algo.general.export.signal;
            figsVars.validStruct.deMul    = uVarsFigs.validStruct.deMul    &&  this.measVars.algo.general.export.deMul && algoVars.useHadamard;
            figsVars.validStruct.reshaped = uVarsFigs.validStruct.reshaped &&  this.measVars.algo.general.export.reshaped;
            
            figsVars.validStruct.rawFFT         = uVarsFigs.validStruct.rawFFT;
            figsVars.validStruct.fittingModel   = uVarsFigs.validStruct.fittingModel;
            figsVars.validStruct.fittedPowerFFT = uVarsFigs.validStruct.fittedPowerFFT;
            figsVars.validStruct.finalFFT       = uVarsFigs.validStruct.finalFFT;
            figsVars.validStruct.calibration    = uVarsFigs.validStruct.calibration;
            
            figsVars.validStruct.phi    = uVarsFigs.validStruct.phi;
            figsVars.validStruct.rawPhi = uVarsFigs.validStruct.rawPhi;
            figsVars.validStruct.phiLog = uVarsFigs.validStruct.phiLog;
            
            figsVars.depthAxType     = uVarsFigs.depthAxType;
            
            if ~this.measVars.algo.general.useFrame
                figsVars.validStruct.unFittedFFT = false;
            end
            
            figsVars.extH = uVarsFigs.extH;
            
            figsVars.numOfChannels = this.measVars.algo.general.channels;
            figsVars.df            = this.measVars.algo.freq.df;
            figsVars.fBar          = gather(this.measVars.algo.freq.fBar);
            figsVars.fIdx          = this.measVars.algo.freq.fUSIdx;
            figsVars.depthVecNorm  = this.measVars.algo.len.depthNorm;
            figsVars.depthVecCntr  = this.measVars.algo.len.depthCntr;
            figsVars.depthVecZero  = this.measVars.algo.len.depthZero;
            figsVars.depthVecIdx   = this.measVars.algo.len.depthIdx;

            figsVars.depthVecNormRaw = this.measVars.algo.len.depthVecAlgo;
            
            figsVars.tVecUS          = this.measVars.algo.timing.tVecUS;
            figsVars.tVecSampleClk   = this.measVars.algo.timing.tVecSampleClk;
            
            if ~this.measVars.algo.general.splitMeas
                figsVars.tVecMeas = this.measVars.algo.timing.tVecMeas;
                figsVars.tVecSig  = this.measVars.algo.timing.tVecSig;
                figsVars.tVecPos  = this.measVars.algo.timing.tVecPosPerFrame; 
            else                
                figsVars.tVecMeas = [];
                figsVars.tVecSig  = [];
                figsVars.tVecPos  = [];
            end
            
            figsVars.analyzeSingleCh = this.measVars.algo.general.analyzeSingleCh;
            figsVars.singleChIdx     = uVarsFigs.singleChIdx;
            figsVars.dispMuEff       = this.measVars.algo.general.calcMuEff;
            
            this.graphics.setUserVars(figsVars);
            this.measVars.figs = figsVars;
        end
        
        function setVars(this, uVars)
            this.uVars = uVars;

            uVarsAO = this.uVars;
            
            fprintf("AOI: ------- Sets Extended Vars From Reduced Vars ----------\n")
            
            %Algorithm:
            algoVars = this.setAlgoVars();
            newExtVars.algo = algoVars;
            
            AO.splitMeas = this.measVars.algo.general.splitMeas;
            AO.splitNum  = this.measVars.algo.general.splitNum;
            
            % US Excitation:
            fGenVars = setFGenVars(this);
            
            % Sampling:
            digiVars = this.setDigitizerVars();
            
            %IO:
            IOVars = this.setIOVars();

            %File System:
            fileSystemVars = this.setFileSystemVars(uVarsAO.fileSystem);

            % Graphics:
            figsVars = this.setGraphicsVars(uVarsAO.figs);

            %-----------------
            % AO
            %-----------------
            
            % Writing Extended Vars to Acousto Optics Object
            fprintf("AOI:  8. Set AO Object Variables.\n");
            % Set AO vars
            AO.useVirtualData      = uVarsAO.ao.useVirtualData;
            AO.virtualDataNoiseSTD = uVarsAO.ao.virtualDataNoiseSTD;
            AO.limitByN            = uVarsAO.ao.limitByN;
            AO.N                   = uVarsAO.ao.N;
            AO.dispTimeTable       = uVarsAO.ao.dispTimeTable;
            AO.skipParamsCheck     = uVarsAO.ao.skipParamsCheck;
            
            AO.displayBuildUp      = uVarsAO.ao.displayBuildUp;
            AO.displayBUEvery      = uVarsAO.ao.displayBUEvery;
            
            AO.uploadToTelegram = uVarsAO.ao.uploadToTelegram;
            AO.telegramChatID   = uVarsAO.ao.telegramChatID;
            
            AO.topUpMode = false;
            AO.topUpIdx  = 0;
            
            AO.useGPU     = uVarsAO.ao.useGPU;
            AO.keepSplits = uVarsAO.ao.keepSplits;
            
            AO.longMeas = this.uVars.ao.timeToSample > 16;

            newExtVars.AO         = AO;
            
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

        % Get Functions
        function aoVars = getVars(this)
%             fprintf("------- Downloading AO Vars ----------\n")
            aoVars.measVars  = this.measVars;
            aoVars.extVars   = this.extVars;
            aoVars.uVars     = this.uVars; 
            
            % We put an empty graphic request in order to shrink the size
            % of the vars mat file.
            aoVars.measVars.figs = AOGraphics.createGraphicsVars();
            aoVars.extVars.figs  = AOGraphics.createOwnerVars();
            aoVars.uVars.figs.extH = [];

            aoVars.measVars.algo.timing.tVecSig  = [];
            aoVars.measVars.algo.timing.tVecMeas = [];
            
            aoVars.measVars.algo.hadamard.sMatInvSingleSqnc = [];
            aoVars.measVars.algo.hadamard.sMatInvSqnc       = [];
        end

        function [res] = getResultStruct(this)
            [~, resStruct, ~, resSplit] = this.algo.getResultsStruct();
            
            if this.measVars.AO.splitMeas
                res = resSplit;
            else
                res = resStruct;
            end  
            
        end
        
        function result = getResult(this)
            result = this.result; 
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
        
        function resetDigitizer(this)
            this.digitizer.resetDigitizer();
            this.configDigitizer();
        end

        % Measurements algorithms types
        function res = runAcoustoOptics(this)
            % reset results memory:
            this.result = [];
            this.rawData = [];
            this.resultNew = [];
            
            res = this.runMeasureAndAnalyze();
        end
        
        function res = runCalibration(this)
            uiwait(msgbox("Please block US beam","Calibration"));
            res = this.runMeasureAndAnalyze();
            this.algo.setCalibrationData(res);
            uiwait(msgbox("Calibration Completed. Please remove US Block","Calibration"));
        end
        
        function res = runTopUp(this)
            this.measVars.AO.topUpMode = true;
            this.measVars.AO.topUpIdx = this.measVars.AO.topUpIdx +1;
            this.resultNew = [];
            res = this.runMeasureAndAnalyze();
        end
        
        function res = runMeasureAndAnalyze(this)
            % Save AO variables (if needed, decided by fileSystem)
            if ~this.measVars.AO.topUpMode
                this.fileSystem.saveVarsToDisk(); 
            end
            keepSplits = this.measVars.AO.keepSplits;
            % Measure data and analyse it with AO Algorithm
            if this.measVars.AO.splitMeas % Split-mode long measurement
                this.stopMeas = false;
                T(this.measVars.AO.splitNum) = this.algo.timeTable;

                if this.measVars.AO.topUpMode
                    offset = this.measVars.AO.splitNum*this.measVars.AO.topUpIdx;
                else
                    this.algo.resetCurSplit(); % if starting a new measurement reset the remains of the older one
                    models = this.algo.getResultsModels();
                    this.splits = models.baseChSplitArr;
                    offset = 0;
                end

                splitNum = offset+this.measVars.AO.splitNum;
                splitStartIdx = offset+1;
                if ~keepSplits
                    this.splits(end, splitNum) = this.splits(end,end);
                end

                this.IO.open();
                for i = splitStartIdx:splitNum
                    fprintf("AOI:");
                    if this.measVars.AO.longMeas
                        
                    end
                    fprintf("Split %d/%d: ", i, splitNum);
                    this.fileSystem.updateSplitInd(i);
                    % Check if stop button was pushed
                    if this.stopMeas
                        this.IO.close();
                        if this.pauseMeas && this.ownedByGUI
                            fprintf("AOI: Measurement Paused.\n");
                            uiwait();
                            this.IO.open();
                            this.pauseMeas = false;
                            this.stopMeas  = false;
                        else
                            fprintf("AOI: Aborting Measurement.\n");
                            break;
                        end
                    end
                    this.algo.initTimeTable();

                    % Bring current Split Datta
                    curSplit = this.measureAndAnalyse();
                    if keepSplits
                        this.splits(:,i) = curSplit;
                    else
                        this.splits = curSplit;
                    end
                    
                    % Full so-far analysis and plotting
                    if  this.measVars.AO.displayBuildUp
                        this.algo.addCurrentSplit(curSplit, i);
                        if i==2 || ~mod(i, this.measVars.AO.displayBUEvery) || i==splitNum
                            this.result = this.algo.reconCurSplit();
                            this.graphics.setData(gather(this.result));
                            this.plotAll();
                            if this.measVars.AO.uploadToTelegram
                                tgprint(this.measVars.AO.telegramChatID, this.graphics.figs.phi.handles.cur.ax, 'photo');
                                tgprintf(this.measVars.AO.telegramChatID, sprintf("AOI: Split %d/%d: ", i, this.measVars.AO.splitNum));
                            end
                        end
                    end

                    T(i) = this.algo.getTimeTable();
                    
                    % No display by default, because no full reconstruction is
                    % performed in split-mode

                    pause(0.001); % For GUI Responsivity
                end
                this.IO.close();

                if ~this.stopMeas
                    this.stopMeas = true;
                    if keepSplits
                        this.algo.initTimeTable();
                        this.result = this.algo.reconSplittedData(this.splits);
                        T(i+1) = this.algo.getTimeTable();
                    end
                else
                    res = [];
                    return
                end
            else % No split
                % Regular single-meas operation:
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
            % In case of split mode, exclude raw Data from saved file as it
            % saved separately during measurement.
            if this.measVars.AO.splitMeas
                saveRes = this.result;
                for i=1:splitNum
                    saveRes.splitRes(i).export.meas = [];
                end
            else
                saveRes = this.result;
            end

            this.timeTable.saveData = tic;
            this.fileSystem.saveResultsToDisk(saveRes)
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
                this.contLive = true;
 
                for i=1:this.measVars.AO.N
                    fprintf("AOI: Live AO: %d/%d.\n", i, this.measVars.AO.N);
                    this.fileSystem.setLiveAOInd(i);
                    if this.measVars.algo.general.useFrame
                        res(i) = this.runAcoustoOptics();
                    else
                        res = this.runAcoustoOptics();
                    end
                    pause(1)
                    if ~this.contLive
                        fprintf("AOI: Live AO was STOPPED.\n");
                        break;
                    end
                end
                
                this.fileSystem.closeLiveAoFileSystem();
                this.result = res;
            else %infinite number of measurements
                %saving is not allowed in unlimited live AO
                this.contLive = true;
                while(this.contLive)
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
                fprintf ("AOI: Generating synthetic data...");
                this.timeTable.acq = tic;
                this.rawData = this.algo.createVirtualData(this.measVars.AO.virtualDataNoiseSTD);
                this.timeTable.acq = toc(this.timeTable.acq);
           elseif this.periAvail.completeHardwareAvail
                % Get data from digitizer (Measure)
                fprintf ("AOI: Acquiring...");
                this.timeTable.acq = tic;
                if ~this.measVars.AO.splitMeas
                    this.IO.open();
                    pause(1);
                end
                this.rawData = this.digitizer.acquire();
                if ~this.measVars.AO.splitMeas
                    this.IO.close();
                end
                this.timeTable.acq = toc(this.timeTable.acq);
            else
                 % No source of data (Error)
                 fprintf("AOI: No source of data found.\n");  
                 res = [];
                 return;
            end
            
            this.timeTable.moveData = tic;
            this.algo.setRawData(this.rawData);
            this.rawData            = gather(this.rawData);
            this.timeTable.moveData = toc(this.timeTable.moveData);
            
            if this.measVars.AO.splitMeas
                this.timeTable.saveSplit = tic;
                this.fileSystem.saveMeasToDisk(this.rawData)
                this.timeTable.saveSplit = toc(this.timeTable.saveSplit);
                this.rawData = [];
            end
            
            fprintf ("Analyzing...");
            this.timeTable.reconstruct = tic;
            res = this.algo.reconstruct();
            this.timeTable.reconstruct = toc(this.timeTable.reconstruct);
            
            fprintf ("Done!\n");
        end
        
        % Misc
        function res = reCalcData(this, uVars)
            uVars.ao.exportData.meas = false;
            this.setVars(uVars);
            this.graphics.setGraphicsDynamicVars();
            keepSplits = this.measVars.AO.keepSplits;
            fprintf ("AOI: Analyzing...");
            if ~isempty(this.rawData)
                % non split, with raw data
                if uVars.ao.useGPU
                    this.algo.setRawData(gpuArray(this.rawData));
                else
                    this.algo.setRawData(this.rawData)
                end
                this.resultNew = this.algo.reconstruct();
            elseif ~this.measVars.AO.splitMeas && isempty(this.rawData) || this.measVars.AO.splitMeas &&  ~isempty(this.result) && ~keepSplits
                % non split without raw data OR split without splits data
                this.resultNew = this.algo.reconFromFFT(this.result);
            elseif  (this.measVars.AO.splitMeas && ~isempty(this.result)) && keepSplits %% old added condition: && ~isempty(this.result.splitRes(1).export.meas)
                % split with splits data 
                splitsNum = length(this.result.splitRes);
                for i=1:splitsNum
                    fprintf("%d...", i);
                    if this.measVars.AO.useGPU 
                        curData = gpuArray(this.result.splitRes(i).export.meas);
                    else
                        curData = this.result.splitRes(i).export.meas;
                    end
                    this.algo.setRawData(curData);
                    this.splits(:,i) = this.algo.reconstruct();
                end
                this.resultNew = this.algo.reconSplittedData(this.splits);
            else
                if isempty(this.result)
                    fprintf("\nAOI: Previous analysis was not complete. There is no results to reCalculate.\n");
                    beep();
                    return;
                else
                    this.resultNew = this.algo.reconExtFFTData(this.result);
                end
            end
            this.graphics.setData(gather(this.resultNew));
            this.plotAll();
            fprintf ("Done!\n");
            res = this.resultNew;
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
            
            if this.ownedByGUI
                this.ownerGUI.displayAOTimeTable();
            end
        end

        %Graphics
        function plotAll(this)
            this.graphics.plotSignal('meas');
            this.graphics.plotSignal('signal');
            this.graphics.plotSignal('deMul');
            this.graphics.plotReshapedSignal();
            this.graphics.plotFFT('rawFFT');
            this.graphics.plotFFT('fittingModel')
            this.graphics.plotFFT('fittedPowerFFT');
            this.graphics.plotFFT('finalFFT');
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
        function data = loadData(this, varargin)
            rawData = []; %#ok<PROPLC> 
            if isstring(varargin{1}) || ischar(varargin{1})
                %argument is file path
                dirPath = varargin{1};
                resPath = sprintf("%s/AO-Results.mat", dirPath);
                res = load(resPath);
                varsPath = sprintf("%s/AO-Vars.mat", dirPath);
                vars = load(varsPath);
                rawDataPath = sprintf("%s/AO-RawData.mat", dirPath);
                if exist(rawDataPath, "file")
                    if length(varargin) ==2 
                        loadRawData = varargin{2};
                    else
                        loadRawData = false;
                    end

                    if loadRawData
                        rawData = load(rawDataPath); %#ok<PROPLC> 
                    end
                end
            elseif isstruct(varargin{1})
                res = varargin{1};
                if isempty(varargin{2})
                    fprintf("Please supply a vars struct.\n");
                    return
                end
                vars = varargin{2};

                if ~isempty(varargin{3})
                    rawData = varargin{3}; %#ok<PROPLC> 
                end
            else
                fprintf("AOI: Unsupported loaded data format.\n");
                return
            end
            
            this.result = res;

            vars.uVarsNew            = this.copyUserVars(vars.uVars); % This is to make sure variables format is correct
            vars.uVarsNew.fileSystem = fileSystemAO.uVarsCreate();    % disable save
            vars.uVarsNew.figs       = AOGraphics.createUserVars();               % Take existing figures configurations
 
            this.rawData = rawData; %#ok<PROPLC> 

            this.setVars(vars.uVarsNew);
            this.graphics.setGraphicsDynamicVars();

            this.graphics.setData(gather(this.result));
            this.plotAll();

            data.res     = res;
            data.vars    = vars;
            data.rawData = rawData; %#ok<PROPLC> 
        end

        function uVars = copyUserVars(this, varsIn)
            uVars = this.createUserVars();
            groups = ["ao", "figs", "fileSystem"];
            for j=1:length(groups)
                fields = fieldnames(uVars.(groups(j)));
                for i=1:length(fields)
                    if isfield(varsIn.(groups(j)), fields{i})
                        uVars.(groups(j)).(fields{i}) = varsIn.(groups(j)).(fields{i});
                    end
                    % else: leave its default value.
                end
            end
        end
        
        function saveVarsToDisk(this, path)
            % path is relative to projPath
            this.fileSystem.saveVarsToDisk(this.getAOVars(), path);
        end
        
        function saveCurrentMeasManual(this, fsVars) 
            fsVars.saveResults = true;
            fsVars.saveVars    = true;
            
            this.setFileSystemVars(fsVars);
            this.fileSystem.configFileSystem();
            
            this.fileSystem.saveVarsToDisk();
            if isstruct(this.resultNew)
                res = this.resultNew;
            else
                res = this.result;
            end
            this.fileSystem.saveResultsToDisk(res);
            this.fileSystem.closeFileSystem();
        end
        
        function printProgress(this, T)
           i = this.vars.idxs(3);
           j = this.vars.idxs(2);
           totalNumOfLines = this.vars.scanSize(2)*this.vars.scanSize(3);
           curNumOfLines = (i-1)*this.vars.scanSize(2)+j;
           percent = curNumOfLines/totalNumOfLines * 100;
           this.vars.averageTime = (this.vars.averageTime *(curNumOfLines-1) +T )/curNumOfLines;
           totalTime = totalNumOfLines * this.vars.averageTime / 60;
           timeToNow = curNumOfLines * this.vars.averageTime / 60;
           str = sprintf("Done line %d/%d (%.2f [%%%%]). Progress: %.2f/%.2f mins.\n",...
                            curNumOfLines, totalNumOfLines, percent, timeToNow, totalTime);
           fprintf(str)
           if this.vars.tg.text && ~mod(this.vars.idxs(2)-1, this.vars.tg.rep)
               tgprintf(this.vars.tg.chatID, str);
           end 
        end
    end
end