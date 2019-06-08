classdef acoustoOptics < handle
    %CONSIST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Objects
        fGen;
        digitizer;
        IO;
        algo;
        graphics;

        % Data
        result   
        timeTable

        % Variables
        redVars  % Reduced from user to acoustoOptics
        extVars  % Extended from acoustoOptics to submodules
        measVars % Measurement vars (after calculation) given from submodules to acoustoOptics
        
        changeLog;
    end
    
    methods (Static)
       
        function redVars = uVarsCreate()
            % algo
            redVars.fSin              = [];              
            redVars.fTrain            = [];
            redVars.cycInPulse        = []; 
            redVars.phantomDepth      = [];
            redVars.distFromPhantom   = [];
            redVars.fExtClk           = []; %fs
            redVars.timeToSample      = [];
            redVars.extClkDcyc        = []; % [%]
            
            % algo & fGen
            redVars.fSclk             = [];     %update in fGen
            
            % algo & digitizer
            redVars.channels          = []; %update in digitizer
            
            % IO
            redVars.IOPort            = [];
            redVars.IOLine            = [];
            
            % General: algo & digitizer
            redVars.useGPU            = []; %algo, digitizer;
            redVars.exportRawData     = [];
        end
        
        function newExtVars = extVarsCreate()
            newExtVars.algo      = Algo.uVarsCreate();
            newExtVars.fGen      = fGen.uVarsCreate();
            newExtVars.digitizer = Digitizer.uVarsCreate();
            newExtVars.IO        = IO.uVarsCreate();
        end
            
        function plotReq = createGraphicsRequest()
            plotReq = Algo.createGraphicsRequest();
        end
        
        function log = checkParamsChanged(oldVars, newVars)
            fields = {'algo', 'fGen', 'digitizer', 'IO'};%, 'graphics'};
            changesFound = false;
            for i = 1:length(fields)
                names = fieldnames(newVars.(fields{i}));
                log.(fields{i}) = false;
                for j = 1:numel(names)
                    log.(fields{i}) = log.(fields{i}) ||... 
                        (isempty(oldVars.(fields{i}).(names{j}))) ||...
                        (~isequal(oldVars.(fields{i}).(names{j}), newVars.(fields{i}).(names{j})));
                        
                    if log.(fields{i}) && strcmp(fields{i}, 'algo') && ...
                             (strcmp(names{j}, 'fSin')       || strcmp(names{j}, 'fTrain') ||...
                              strcmp(names{j}, 'cycInPulse') || strcmp(names{j}, 'channels'))
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
            fprintf("------- Creating AcoustoOptic ----------\n")
            fprintf(" 1. Creating an Arbitrary Function Generator Object\n")
            this.fGen = fGen();
            fprintf(" 2. Creating a Digitizer Object\n")
            this.digitizer = Digitizer();
            fprintf(" 3. Creating an IO Object\n")
            this.IO = IO();
            fprintf(" 4. Creating an Algorithm Object\n")
            this.algo = Algo();
            
            this.redVars = this.uVarsCreate();
            this.extVars = this.extVarsCreate();
        end 

        function init(this)
            fprintf("------- Initiating AcoustoOptic ----------\n")
            fprintf(" 1. Resetting Instruments\n")
            instrreset;
            fprintf(" 2. Connecting to fGen\n")
            % Configure and Activate fGen
            this.fGen.connect()
            fprintf(" 3. Connecting to Digitizer\n")
            % Start the digitizer
            this.digitizer.connect();
            fprintf(" 4. Connecting to IO\n")
            % Start the IO
            this.IO.connect();       
        end
        
        function setMeasVars(this, reducedVars)
%             user gives reduced variables. this function convert these 
%             variables to extended variables suitable for each
%             submodule\peripheral. this function extend the variables by
%             adding default acousto optics values to the user reduced
%             variables
            fprintf("------- Sets Extended Vars From Reduced Vars ----------\n")
            fprintf(" 1. Sets Algorithm reduced Vars\n")
            % Algo 
            algoVars      = this.algo.uVarsCreate();

            algoVars.c               = 1550;
            algoVars.bufferSizeBytes = 8*2^20;
            algoVars.bytesPerSample  = 2;
            algoVars.fSin            = reducedVars.fSin;
            algoVars.fTrain          = reducedVars.fTrain;
            algoVars.cycInPulse      = reducedVars.cycInPulse;
            algoVars.phantomDepth    = reducedVars.phantomDepth;
            algoVars.distFromPhantom = reducedVars.distFromPhantom;
            algoVars.fExtClk         = reducedVars.fExtClk;
            algoVars.timeToSample    = reducedVars.timeToSample;
            algoVars.extClkDcyc      = reducedVars.extClkDcyc;
            algoVars.exportRawData   = reducedVars.exportRawData;

            fprintf(" 2. Sets Function Generator reduced Vars\n")
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
            algoVars.fSclk      = reducedVars.fSclk; %update in fGen
            fGenVars.ch{1}.Sclk = algoVars.fSclk;
            fGenVars.ch{2}.Sclk = algoVars.fSclk; 

            fprintf(" 3. Sets Digitizer reduced Vars\n")
            % digitizer 
            digitizerVars = this.digitizer.uVarsCreate();

            digitizerVars.mode = 'TS'; 
            digitizerVars.fs   = algoVars.fExtClk;
            
            % Algo & Digitizer
            algoVars.preTriggerSamples  = digitizerVars.preTriggerSamples;
            algoVars.useGPU             = reducedVars.useGPU;
            algoVars.channels           = reducedVars.channels;

            digitizerVars.useGPU            = algoVars.useGPU;
            digitizerVars.channels          = algoVars.channels;

            % IO
            IOVars = this.IO.uVarsCreate;

            IOVars.mode = 0; % Output Only
            IOVars.port = reducedVars.IOPort;
            IOVars.line = reducedVars.IOLine;

            % Calculate measurement variables and dimensions
            fprintf(" 4. Calculating Algorithm Variables\n")
            this.measVars.algo = this.algo.updateAlgoUserVars(algoVars);

            % Complete digitizer variables from calculated algorithm vars
            fprintf(" 5. Extracting Digitizer Variables from Algorith,\n")
           
            digitizerVars.numOfBuffers          = this.measVars.algo.digitizer.numOfBuffers;
            digitizerVars.bufferSizeBytes       = this.measVars.algo.digitizer.bufferSizeBytes;
            digitizerVars.samplesPerAcq         = this.measVars.algo.digitizer.samplesPerAcq;
            digitizerVars.samplesPerBuffer      = this.measVars.algo.digitizer.samplesPerBuffer;
            digitizerVars.samplesPerAcqAllCh    = this.measVars.algo.digitizer.samplesPerAcqAllCh;
            digitizerVars.samplesPerBufferAllCh = this.measVars.algo.digitizer.samplesPerBufferAllCh;
           
            % Writing Extended Vars to Acousto Optics Object
            newExtVars.algo      = algoVars;
            newExtVars.fGen      = fGenVars;
            newExtVars.digitizer = digitizerVars;
            newExtVars.IO        = IOVars;

            fprintf(" 6. Checking For Changes,\n")
            this.changeLog = this.checkParamsChanged(this.extVars, newExtVars);

            % set extended
            this.extVars = newExtVars;

        end
        
        function configPeripherals(this)
            fprintf("------- Configuring Peripheral ----------\n")

            % Download signal to fGen and Activate
            if this.changeLog.fGen
                fprintf(" ** Configuring fGen\n")
                this.configfGen();
            end
            
            % Config digitizer
            if this.changeLog.digitizer
                fprintf(" ** Configuring Digitizer\n")
                this.configDigitizer();
            end
            
             % Config the IO
            if this.changeLog.IO
                fprintf(" ** Configuring IO\n")
                this.configIO()
            end
            
            fprintf(" ** Done configuring Peripherals\n")
        end

        function res = measureAndAnlayse(this)
            this.timeTable.acq = tic;
            this.IO.open();
            bufferDataOut = this.digitizer.acquireDataTS();
            this.IO.close();
            this.timeTable.acq = toc(this.timeTable.acq);
            
            this.timeTable.analyse = tic;
            res = this.algo.analyse(bufferDataOut);
            this.result = res;
            this.timeTable.analyse = toc(this.timeTable.analyse);
            
            this.timeTable = this.getInnerTimeTables();
        end
        
        function res = fullMeasureAndAnalyse(this, uVars)
            this.setMeasVars(uVars);
            this.configPeripherals();
            res = this.measureAndAnlayse();
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
        end
    
        function configDigitizer(this)
            this.digitizer.setDigitizerVars(this.extVars.digitizer);
            this.digitizer.configure();
        end
        
        function configIO(this)
           this.IO.allocPorts(this.extVars.IO); 
        end
        
        function resetDigitizerMem(this)
         	this.digitizer.releaseBuffers();  
        end
        
        function updateSubmodulesVars(this)
            if log.fGen.changed
                this.configPeripherals()
            elseif log.algo.changed
                this.measVars.algo = this.algo.updateAlgoUserVars(this.extVars.algo);
                this.algo.resetAlgoArrays();
                this.configDigitizer();
            end
        end
        
        function setSamplingTime(this, timeToSample)
            this.algo.uVars.timeToSample = timeToSample;
            this.vars.algo = this.algo.updateAlgoUserVars(this.vars.algo.uVars);
        	this.configDigitizer();
        end

        function measVars = getMeasVars(this)
            fprintf("------- Downloading Measurement Vars ----------\n")
            measVars.algo      = this.measVars.algo;
            measVars.fGen      = this.measVars.fGen;
            measVars.digitizer = this.measVars.digitizer;
            measVars.IO        = this.measVars.IO;
        end
        
        function algoVars = getAlgoVars(this)
            algoVars = this.algo.getAlgoVars();
        end
        
        function timeTable = getTimeTable(this)
           timeTable = this.timeTable; 
        end
        
        function timeTable = getInnerTimeTables(this)
            timeTable.algo      = this.algo.getTimeTable();
            timeTable.digitizer = this.digitizer.getTimeTable();
        end
        
        function prepareGraphics(this)
            this.algo.directGraphics(); 
        end
        
    end
end

                % there is a request
                
                % what was previous:
                % if previous was internal - use the same figs, open and close
                % according to new request.
                % if it was external - create the figures
                %
                % if the current request is external
                % what was the previous?
                % if it was internal - close all figures, copy new
                % plotRequests
                % if it was external - copy new plot requests.

%         function prepareGraphics(this) 
%             fprintf("------- Preparing Algorithm Graphics Requests ----------\n")
%             if this.graphics.uVars.internal % if the current request is internal - 
%                 ch = this.graphics.uVars.ch;
%                 pos = this.graphics.uVars.pos;
%                 names = fieldnames(this.graphics.plotRequests);
%                 for i=1:length(names)
%                     if this.vars.uVarsOld.graphics.internal && this.vars.uVarsOld.graphics.mask(i) % if previous was internal
%                         close(this.graphics.plotRequests.(names{i}).getFigHandle()); %close old figs
%                     end
%                     if this.graphics.uVars.mask(i)
%                         fh = figure(); ax = axes(); %open new figs
%                         this.graphics.plotRequests.(names{i}).setParams(...
%                                              this.graphics.uVars.mask(i), fh, ax, ch, pos)
%                     end
%                 end
%             else % if the current request is external
%                 names = fieldnames(this.graphics.plotRequests);
%                 for i=1:length(names)
%                     if this.vars.uVarsOld.graphics.internal && this.vars.uVarsOld.graphics.mask(i) % if previous was internal
%                         close(this.graphics.plotRequests.getFigHandle()); %close old figs
%                     end
%                 end
%                 this.graphics.plotRequests = this.graphics.uVars.externalPlotRequests;
%             end
%             pause(0.1);
%             this.algo.setGraphics(this.graphics.plotRequests)
%         end