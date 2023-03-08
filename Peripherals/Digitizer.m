classdef Digitizer < handle
    %DIGITIZER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vars
        system
        buffers
        NPT
        TS
        CS
        
        bufferDataOut
        uData
        Data;
        croppedData;
        
        timeTable
        alazarDefs
        hardwareAvailable
    end
    
    methods (Static)
        
        function uVars = uVarsCreate() 
            uVars.mode                  = [];
            uVars.channels              = [];
            uVars.inputRange            = '1V';
            uVars.useGPU                = false;
            uVars.fs                    = [];
            uVars.bufferSizeBytes       = [];
            uVars.samplesPerMeas        = [];
            uVars.exportCropped         = [];
            uVars.croppedSamples        = [];
            uVars.coupling              = 1; %0=DC1=AC
            uVars.draw                  = false;
        end   
        
        function uVars = uVarsMonitorCreate()
            uVars.fs           = [];
            uVars.channels     = []; 
            uVars.timeToSample = [];
            uVars.triggerDelay = 0;
            uVars.coupling     = 1; %0=DC1=AC
            uVars.avgNum       = [];
            
            uVars.extClk  = true;
            uVars.extTrig = true;
        end
    end
    
    methods
        
        %Init Functions
        function this = Digitizer()
            this.system.boardHandle       = [];
            this.system.systemId          = int32(1);
            this.system.boardId           = int32(1);
            this.system.channelsMask      = 0;
            this.system.bufferCount       = 4;
            this.system.inputRange        = '1V';
            this.system.voltsRange        = 1;
            this.system.bytesPerSample    = 2;
            this.system.preTriggerSamples = 14;
            this.system.allocatedBuffers = false;
            this.buffers = cell(1, this.system.bufferCount);
            this.vars.figs.hFig    = 1;
            this.vars.figs.hAx     = 1;
            this.vars.figs.hPlot   = 1;
            this.vars.prevChannels = 0;
            this.vars.channels     = 0;
        end
        
        function initDefs(this)
            AlazarDefs;
            names = who();
            for i = 1:length(names)
                if strcmp(names{i}, 'this'); continue; end
                this.alazarDefs.(names{i}) = eval(names{i});
            end
        end
        
        function status = connect(this)
            status = true;
            this.timeTable.Connect = tic;
            
            try
                this.initDefs();
                this.system.hardwareAvailable = true;
            catch
                fprintf("DIGITIZER: No DEFS found, Assuming no digitizer on that system.\n");
                this.system.hardwareAvailable = false;
                status = false;
                return;
            end
            
            if ~alazarLoadLibrary() 
                fprintf('DIGITIZER: Error: ATSApi library not loaded\n'); 
                status = false;
                this.system.hardwareAvailable = false;
                return; 
            end
            
            this.system.boardHandle = AlazarGetBoardBySystemID(this.system.systemId, this.system.boardId);
            setdatatype(this.system.boardHandle, 'voidPtr', 1, 1);
            
            this.timeTable.Connect = toc(this.timeTable.Connect);
        end

        function resetDigitizer(this)
            unloadlibrary 'ATSApi';
            this.connect();
        end
        
        % Variables Functions
        function retVars = setVars(this, vars)
            this.vars.prevChannles = this.vars.channels;
            
            % General Vars
            this.system.bufferCount = vars.channels;
            this.vars.useGPU        = vars.useGPU;

            this.vars.fs       = vars.fs;
            this.vars.channels = vars.channels;
            
            % Utility vars
            channels          = this.vars.channels;
            bytesPerSample    = this.system.bytesPerSample;
            drawData = vars.draw;
            
            switch vars.mode
                case 'TS'
                    % Calculate How many samples to acquire according to
                    % user request and HW limitations.
                    
                    % Pre Trigger samples are bug.
%                     preTriggerSamples = this.calcPreTriggerSamplesTS();
                    preTriggerSamples = 0;
                    % User vars
                    bufferSizeBytes  = vars.bufferSizeBytes;          
                    samplesPerMeas   = vars.samplesPerMeas;
                    
                    acqSize = samplesPerMeas * bytesPerSample * channels;
                    
                    if acqSize < bufferSizeBytes
                        bufferSizeBytes = acqSize;
                    end

                    % Calculate TS vars
                    samplesPerMeas2       = samplesPerMeas        + preTriggerSamples;
                    samplesPerBufferAllCh = bufferSizeBytes       / bytesPerSample;
                    samplesPerBuffer      = samplesPerBufferAllCh / channels;
                    numOfBuffers          = ceil(samplesPerMeas2  / samplesPerBuffer);
                    samplesPerAcq         = samplesPerBuffer      * numOfBuffers;
                    samplesPerAcqAllCh    = samplesPerAcq         * channels;
                    
                    % Excessive samples
                    excSamplesPre   = preTriggerSamples;
                    excSamplesPost  = samplesPerAcq - samplesPerMeas - excSamplesPre;
                    excSamplesIdx   = samplesPerMeas2;
                    
                    exportData = true;
                    
                    % Collect TS Vars
                    this.TS.samplesPerMeas     = samplesPerMeas;
                    this.TS.samplesPerMeas2    = samplesPerMeas2;
                    this.TS.samplesPerAcq      = samplesPerAcq;
                    this.TS.samplesPerAcqAllCh = samplesPerAcqAllCh;
                    
                    this.TS.excSamplesPre  = excSamplesPre;
                    this.TS.excSamplesPost = excSamplesPost;
                    this.TS.excSamplesIdx  = excSamplesIdx;

                    this.TS.exportCropped  = vars.exportCropped;
                    this.TS.croppedSamples = vars.croppedSamples;
                    
                    drawData = false;
                case 'NPT'
                    preTriggerSamples = 0;
                    preRecordSamples = this.calcPreRecordSampleNPT();
                    
                    %User vars
                    postTriggerSamplesUser = round(vars.timeToSample*vars.fs);
                    postTriggerSamplesArt  = preTriggerSamples + preRecordSamples + postTriggerSamplesUser;
                    postTriggerSamples     = max(256, ceil(postTriggerSamplesArt/128)*128); % Align to 128;
                    
                    recordsPerBuffer = vars.avgNum;
    
                    % Calculate NPT vars
                    samplesPerRecord      = postTriggerSamples;
                    samplesPerBuffer      = samplesPerRecord * recordsPerBuffer;
                    samplesPerBufferAllCh = samplesPerBuffer * channels;
                    bufferSizeBytes       = bytesPerSample   *  samplesPerBufferAllCh;

                    if vars.numMeas == inf
                        numOfBuffers = hex2dec('7FFFFFFF');
                        recordsPerAcuisition = hex2dec('7FFFFFFF');
                        exportData = false;
                        samplesPerAcq = 0;
                        samplesPerAcqAllCh = 0;
                    else
                        numOfBuffers = vars.numMeas; %non stop acuisition = hex2dec('7FFFFFFF')
                        recordsPerAcuisition = recordsPerBuffer *  numOfBuffers;
                        exportData = true;
                        samplesPerAcq = samplesPerBuffer * numOfBuffers;
                        samplesPerAcqAllCh = samplesPerAcq * channels;
                    end
                    
                    postRecordSamples = samplesPerRecord - postTriggerSamplesArt;

                    % Collect NPT Vars
                    this.NPT.postTriggerSamplesUser = postTriggerSamplesUser;
                    this.NPT.preTriggerSamples      = preTriggerSamples;
                    this.NPT.postTriggerSamples     = postTriggerSamples;
                    this.NPT.samplesPerRecord       = samplesPerRecord;
                    this.NPT.recordsPerBuffer       = recordsPerBuffer;
                    this.NPT.recordsPerAcuisition   = recordsPerAcuisition;
                    
                    this.NPT.preRecordSamples       = preRecordSamples;
                    this.NPT.postRecordSamples      = postRecordSamples;

                    this.NPT.cutIdxs = [1:preRecordSamples, (postTriggerSamplesArt+1):samplesPerRecord];
                    
                    if drawData
                        this.vars.figs.tVec = (0:1:(samplesPerRecord-1))/this.vars.fs *1e6;
                        this.vars.figs.fVec = (this.vars.fs/samplesPerRecord) *  ( (-samplesPerRecord/2) : 1 : (samplesPerRecord/2)-1 ) /1e6;
                    end

                case 'ContNPT'
                    % Measure continuesly numerous triggers and perform no
                    % averaging.
                    
                    % General predefined variables:
%                     maxBufferSize = 2^23; % 8MB
                    maxBufferSize = 2^25; % 64MB
                    preRecordSamples = this.calcPreRecordSampleNPT();

                    % User vars:
                    postTriggerSamplesUser = round(vars.timeToSample*vars.fs);
                    userRecordSizeBytes = postTriggerSamplesUser*bytesPerSample*channels;

                    % Calculate how mant triggers acquisition are in a
                    % buffer:
                    recordsPerBufferRaw    = floor(maxBufferSize / userRecordSizeBytes);

                    if recordsPerBufferRaw == 0 
                        fprintf("DIGITIZER: NOTICE: measured time is longer than max suggested buffer size.\n")
                        recordsPerBufferRaw = 1;
                    end
                    
                    % How many buffer are required to capture total number 
                    % of triggers:
                    numOfBuffers = ceil(vars.numMeas/recordsPerBufferRaw);
                    
                    % How many redundant record will be captured:
                    recordsPerAcqRaw = numOfBuffers*recordsPerBufferRaw;
                    redundantRecord = recordsPerAcqRaw - vars.numMeas;

                    % Make many short records into several long continues
                    % with size equivalent to buffer size:
                    % (to spare the missed triggers between NPTs)
                    postTriggerSamplesCont = recordsPerBufferRaw*postTriggerSamplesUser;

                    postTriggerSamplesArt = postTriggerSamplesCont + preRecordSamples;
                    samplesPerRecord      = max(256, ceil(postTriggerSamplesArt/128)*128);
                    samplesPerRecordAllCh = samplesPerRecord * channels;
                    recordSizeBytes       = samplesPerRecordAllCh * bytesPerSample;
                    postRecordSamples     = samplesPerRecord - postTriggerSamplesArt;

                    % Common variables:
                    recordsPerBuffer      = 1;
                    bufferSizeBytes       = recordSizeBytes;
                    samplesPerBuffer      = samplesPerRecord; 
                    samplesPerBufferAllCh = samplesPerRecordAllCh;
                    recordsPerAcuisition  = numOfBuffers;
                    
                    samplesPerAcq = samplesPerBuffer * numOfBuffers;
                    samplesPerAcqAllCh = samplesPerAcq * channels;

                    exportData = true;

                    % Collect NPT Vars
                    this.NPT.postTriggerSamplesUser = postTriggerSamplesUser;
                    this.NPT.preTriggerSamples      = 0;
                    this.NPT.postTriggerSamples     = samplesPerRecord;
                    this.NPT.samplesPerRecord       = samplesPerRecord;
                    this.NPT.recordsPerBuffer       = recordsPerBuffer;
                    this.NPT.recordsPerAcuisition   = recordsPerAcuisition;
                    this.NPT.preRecordSamples       = preRecordSamples;
                    this.NPT.postRecordSamples      = postRecordSamples;
                    this.NPT.redundantRecord        = redundantRecord;

                    this.NPT.cutIdxs = [1:preRecordSamples, postTriggerSamplesArt+1:samplesPerRecord];
                    this.NPT.samplesPerAcqAllRecords = postTriggerSamplesUser * recordsPerAcqRaw;
                    this.NPT.extraRecordIdxs = (postTriggerSamplesUser*vars.numMeas+1) :  this.NPT.samplesPerAcqAllRecords;
            end
            
            this.calcTriggerDelayAlignment()

            this.vars.bufferSizeBytes       = bufferSizeBytes;
            this.vars.numOfBuffers          = numOfBuffers;
            this.vars.samplesPerBuffer      = samplesPerBuffer;
            this.vars.samplesPerBufferAllCh = samplesPerBufferAllCh;
            this.vars.samplesPerAcq         = samplesPerAcq;
            this.vars.samplesPerAcqAllCh    = samplesPerAcqAllCh;

            this.vars.exportData = exportData;
            this.vars.drawData   = drawData;
            
            % Set System Vars
            this.system.extClk                = vars.extClk;
            this.system.extTrig               = vars.extTrig;
            this.system.triggerLevel_volts    = 3; % external trigger level
            this.system.triggerRange_volts    = 5; % external trigger input range 
            this.system.triggerLevel_code     =(128 + 127 * this.system.triggerLevel_volts / this.system.triggerRange_volts);
            this.system.mode                  = vars.mode;
            
            this.system.triggerDelay_sec_user      = vars.triggerDelay;
            this.system.triggerDelay_samples_user  = floor(this.system.triggerDelay_sec_user * this.vars.fs + 0.5);
            this.system.triggerDelay_samples       = uint32(max(this.system.triggerDelayAlignment, 2^(floor(log2(this.system.triggerDelay_samples_user)))));
            this.system.triggerDelay_samples       = this.system.triggerDelay_samples_user;
            
            this.vars.triggerDelayManualSamples    = this.system.triggerDelay_samples_user - this.system.triggerDelay_samples +1;
            
            this.system.triggerTimeout_sec    = 0;
            this.system.triggerTimeout_clocks = uint32(floor(this.system.triggerTimeout_sec / 10.e-6 + 0.5));
            this.system.inputRange            = vars.inputRange;
            
            this.system.coupling = vars.coupling;
            
            this.setInputRange();

            switch vars.mode
                case 'TS'
                    this.TS.samplesPerBuffer      = samplesPerBuffer;
                    this.TS.samplesPerBufferAllCh = samplesPerBufferAllCh;
                    this.TS.numOfBuffers          = numOfBuffers;
                    this.CS                       = [];
                    this.NPT                      = [];
                case 'NPT'
                    this.TS = [];
                    this.CS = [];
                case 'ContNPT'
                    this.TS = [];
                    this.CS = [];
                case 'CS'
                    this.TS = [];
                    this.NPT = [];
            end
            
            retVars = this.vars;
        end
        
        function preRecordSamples = calcPreRecordSampleNPT(this)
            switch this.vars.channels
                case 1
                     preRecordSamples = 7;
                case 2
                     preRecordSamples = 6;
                case 4
                     preRecordSamples = 10;
                case 8
                     preRecordSamples = 10;
                case 16
                     preRecordSamples = 11;
            end
        end

        function preTriggerSamples = calcPreTriggerSamplesTS(this)
                if this.vars.channels == 1
                   preTriggerSamples = 7;
                elseif this.vars.channels == 2
                    preTriggerSamples = 6;
                elseif this.vars.channels == 4
                    preTriggerSamples = 10;
                elseif this.vars.channels == 8
                    preTriggerSamples = 10;
                elseif this.vars.channels == 16
                    preTriggerSamples = 11;
                end
            this.system.preTriggerSamples = preTriggerSamples;
        end
        
        function calcTriggerDelayAlignment(this)
            channels = this.vars.channels;
            if channels ==1
                this.system.triggerDelayAlignment = 16;
            elseif channels == 2
                this.system.triggerDelayAlignment = 8;
            elseif channels == 4
                this.system.triggerDelayAlignment = 4;
            elseif channels == 8
                this.system.triggerDelayAlignment = 2;
            elseif channels == 16
                this.system.triggerDelayAlignment = 1;
            end
        end
        
        function setInputRange(this)
            switch this.system.inputRange
                case '20mV'
                    this.system.rangeDef = this.alazarDefs.INPUT_RANGE_PM_20_MV;
                    this.system.voltsRange = 40e-3;
                case '50mV'
                    this.system.rangeDef  = this.alazarDefs.INPUT_RANGE_PM_50_MV;
                    this.system.voltsRange = 100e-3;
                case '100mV'
                    this.system.rangeDef  = this.alazarDefs.INPUT_RANGE_PM_100_MV;
                    this.system.voltsRange = 200e-3;
                case '200mV'
                    this.system.rangeDef  = this.alazarDefs.INPUT_RANGE_PM_200_MV;
                    this.system.voltsRange = 400e-3;
                case '500mV'
                    this.system.rangeDef  = this.alazarDefs.INPUT_RANGE_PM_500_MV;
                    this.system.voltsRange = 1000e-3;
                case '1V'
                    this.system.rangeDef  = this.alazarDefs.INPUT_RANGE_PM_1_V;
                    this.system.voltsRange = 2;
                case '2V'
                    this.system.rangeDef  = this.alazarDefs.INPUT_RANGE_PM_2_V;
                    this.system.voltsRange = 4;
                case '5V'
                    this.system.rangeDef  = this.alazarDefs.INPUT_RANGE_PM_5_V;
                    this.system.voltsRange = 10;
                otherwise
                    this.system.rangeDef  = this.alazarDefs.INPUT_RANGE_PM_1_V;
                    this.system.voltsRange = 2;
                    fprintf("'DIGITIZER: Does not support in requested input range. Default input range: 1V\n");
            end
        end
        
        function vars = getVars(this)
           vars.system  = this.system;
           vars.system.boardHandle = [];

           vars.vars = this.vars;
           vars.TS   = this.TS;
           vars.NPT  = this.NPT;
           vars.CS   = this.CS;
        end
        
        %Time Table Functions
        function initTimeTable(this)
            this.timeTable = struct();
        end
        
        function timeTable = getTimeTable(this)
            timeTable = this.timeTable;
        end
        
        %Config Functions
        function status = configure(this)
            %--- Configure Clk & Trig ---
            status = this.confExtClk();  if(~status); return; end
            status = this.confChnls();   if(~status); return; end
            status = this.confTrig();    if(~status); return; end

            %--- Configure Dimensions ---
            status = this.setRecordSize(); if(~status); return; end
            status = this.setADMA(); if(~status); return; end

            % --- Create pointers to CPU memory Buffers ---
            this.allocateDataMat();
            status = this.allocateBuffers();
        end
            
        function status = confExtClk(this)
            this.timeTable.confExtClk = tic;
            
            if this.system.extClk
                clkSource = this.alazarDefs.FAST_EXTERNAL_CLOCK;
                clkRate   = this.alazarDefs.SAMPLE_RATE_USER_DEF;
            else
                clkSource = this.alazarDefs.INTERNAL_CLOCK;
                if this.vars.fs == 1e6
                    clkRate = this.alazarDefs.SAMPLE_RATE_1MSPS;
                elseif this.vars.fs == 2e6
                    clkRate = this.alazarDefs.SAMPLE_RATE_2MSPS;
                elseif this.vars.fs == 5e6
                    clkRate = this.alazarDefs.SAMPLE_RATE_5MSPS;
                elseif this.vars.fs == 10e6
                    clkRate = this.alazarDefs.SAMPLE_RATE_10MSPS;
                elseif this.vars.fs == 20e6
                    clkRate = this.alazarDefs.SAMPLE_RATE_20MSPS;
                elseif this.vars.fs == 50e6
                    clkRate = this.alazarDefs.SAMPLE_RATE_50MSPS;
                elseif this.vars.fs == 100e6
                    clkRate = this.alazarDefs.SAMPLE_RATE_100MSPS;
                end
            end
            
            status = true;
            retCode = AlazarSetCaptureClock(this.system.boardHandle,  ... % HANDLE -- board handle
                                            clkSource,   ... % U32 -- clock source id
                                            clkRate,     ... % U32 -- sample rate id
                                            this.alazarDefs.CLOCK_EDGE_RISING,        ... % U32 -- clock edge id
                                            0);                          % U32 -- clock decimation                                
            if retCode ~= this.alazarDefs.ApiSuccess; fprintf('DIGITIZER: Error: AlazarSetCaptureClock failed -- %s\n', errorToText(retCode)); status = false; return; end
            
            this.timeTable.confExtClk = toc(this.timeTable.confExtClk);
            fprintf('DIGITIZER: Done configurations.\n')
        end
        
        function status = confChnls(this)
            this.timeTable.confChnls = tic;
            status = true;

            inputRange  = this.system.rangeDef;
            impedance   = this.alazarDefs.IMPEDANCE_50_OHM;
            successCode = this.alazarDefs.ApiSuccess;
            
            if this.system.coupling
                coupling    = this.alazarDefs.AC_COUPLING;
            else
                coupling    = this.alazarDefs.DC_COUPLING;
            end
            
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_A, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_B, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_C, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_D, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_E, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_F, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_G, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_H, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_I, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_J, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_K, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_L, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_M, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_N, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_O, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_P, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('DIGITIZER: Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end

            this.system.channelMask = 0;
            for i = 1:this.vars.channels
                this.system.channelMask = bitor(bitshift(this.system.channelMask,1),1);
            end
            
            this.timeTable.confChnls = toc(this.timeTable.confChnls);
        end
        
        function status = confTrig(this)
            this.timeTable.confExtTrig = tic;
            
            if this.system.extTrig
                trigSource = this.alazarDefs.TRIG_EXTERNAL;
            else
                trigSource = this.alazarDefs.TRIG_CHAN_A;
            end
            
            status = true;
            % Disable internal Triggers and set External Trigger
            retCode = AlazarSetTriggerOperation(this.system.boardHandle, this.alazarDefs.TRIG_ENGINE_OP_J,...
                                                this.alazarDefs.TRIG_ENGINE_J, trigSource, this.alazarDefs.TRIGGER_SLOPE_POSITIVE, 128,... % TODO: Test level code
                                                this.alazarDefs.TRIG_ENGINE_K, this.alazarDefs.TRIG_DISABLE, this.alazarDefs.TRIGGER_SLOPE_POSITIVE, 128 );                                                    
            if retCode ~= this.alazarDefs.ApiSuccess ; fprintf('DIGITIZER: Error: AlazarSetTriggerOperation failed -- %s\n', errorToText(retCode)); status = false; return; end
            
            if this.system.extTrig
                % Set External trigger Configurations
                retCode = AlazarSetExternalTrigger( this.system.boardHandle, this.alazarDefs.DC_COUPLING, this.alazarDefs.ETR_TTL );
                if retCode ~= this.alazarDefs.ApiSuccess; fprintf('DIGITIZER: Error: AlazarSetExternalTrigger failed -- %s\n', errorToText(retCode)); status = false; return; end
            end
            
            % Trigger Delay set to 0
            retCode = AlazarSetTriggerDelay(this.system.boardHandle, this.system.triggerDelay_samples);
            if retCode ~= this.alazarDefs.ApiSuccess; fprintf('DIGITIZER: Error: AlazarSetTriggerDelay failed -- %s\n', errorToText(retCode)); status = false; return; end

            % Set Trigger Timeout 
            retCode = AlazarSetTriggerTimeOut(this.system.boardHandle, this.system.triggerTimeout_clocks); % U32 -- timeout_sec / 10.e-6 (0 == wait forever)       
            if retCode ~= this.alazarDefs.ApiSuccess; fprintf('DIGITIZER: Error: AlazarSetTriggerTimeOut failed -- %s\n', errorToText(retCode)); status = false; return; end    
            
            this.timeTable.confExtTrig = toc(this.timeTable.confExtTrig);
        end
        
        function status = setADMA(this)
            this.timeTable.setADMA = tic;
            
            switch this.system.mode
                case 'TS'
                    admaFlags = this.alazarDefs.ADMA_EXTERNAL_STARTCAPTURE + this.alazarDefs.ADMA_TRIGGERED_STREAMING;% + this.alazarDefs.ADMA_FIFO_ONLY_STREAMING;
                    preTriggerSamples     = 0;
                    samplesPerRecord      = this.TS.samplesPerBuffer;
                    recordsPerBuffer      = 1;
                    recordsPerAcquisition =  hex2dec('7FFFFFFF');
                case {'NPT', 'ContNPT'}
                    admaFlags = this.alazarDefs.ADMA_EXTERNAL_STARTCAPTURE + this.alazarDefs.ADMA_NPT; % Select AutoDMA flags as required
                    preTriggerSamples     = -int32(this.NPT.preTriggerSamples);
                    samplesPerRecord      = this.NPT.samplesPerRecord;
                    recordsPerBuffer      = this.NPT.recordsPerBuffer;
                    recordsPerAcquisition = this.NPT.recordsPerAcuisition;
            end
            
            status = true;
            retCode = AlazarBeforeAsyncRead(this.system.boardHandle, this.system.channelMask, preTriggerSamples, samplesPerRecord, recordsPerBuffer, recordsPerAcquisition, admaFlags);
            if retCode ~= this.alazarDefs.ApiSuccess; fprintf('DIGITIZER: Error: AlazarBeforeAsyncRead failed -- %s\n', errorToText(retCode)); status = false; return; end
            this.timeTable.setADMA = toc(this.timeTable.setADMA);
        end

        function status = setRecordSize(this) 
            status = true;
            if strcmp(this.system.mode, 'NPT') || strcmp(this.system.mode, 'ContNPT')
                retCode = AlazarSetRecordSize(this.system.boardHandle, this.NPT.preTriggerSamples, this.NPT.postTriggerSamples);
                if retCode ~= this.alazarDefs.ApiSuccess; fprintf('DIGITIZER: Error: AlazarSetRecordSize failed -- %s\n', errorToText(retCode)); status = false; return; end
            end
        end
        
        % Memory Manage Functions       
        function status = allocateBuffers(this)
            this.timeTable.allocateBuffers = tic;

            if this.system.allocatedBuffers
                this.releaseBuffers();
            end

            status = true;
            for j = 1 : this.system.bufferCount
                pbuffer = AlazarAllocBuffer(this.system.boardHandle, this.vars.bufferSizeBytes);
%                 pbuffer1 = libpointer('uint8Ptr', this.bufferDataOut(j,:));
                if pbuffer == 0; fprintf('DIGITIZER: Error: AlazarAllocBuffer %u samples failed\n', SamplingCard.samplesPerBuffer); status = false; return; end
                this.buffers(1, j) = { pbuffer };
            end

            this.system.allocatedBuffers = true;
            
            this.timeTable.allocateBuffers = toc(this.timeTable.allocateBuffers);
        end
        
        function allocateDataMat(this)
            this.timeTable.allocMemoryBuffers = tic;
            this.bufferDataOut = [];
            switch this.system.mode
                case 'TS'
                    this.bufferDataOut = zeros(this.TS.numOfBuffers, this.TS.samplesPerBufferAllCh, 'uint16');
                case 'NPT'
                    if this.vars.exportData
                        this.bufferDataOut = zeros(this.NPT.postTriggerSamplesUser, this.vars.channels, this.vars.numOfBuffers, 'single');
                    end
                case 'ContNPT'
                    this.bufferDataOut = zeros(this.vars.numOfBuffers, this.vars.samplesPerBufferAllCh, 'uint16');
            end
            this.timeTable.allocMemoryBuffers = toc(this.timeTable.allocMemoryBuffers);  
        end
        
        function status = postBuffersToBoard(this)
            this.timeTable.postBuffersToBoard = tic;
            status = true;
            for bufferIndex = 1 : this.system.bufferCount
                pbuffer = this.buffers{1, bufferIndex};
                retCode = AlazarPostAsyncBuffer(this.system.boardHandle, pbuffer, this.vars.bufferSizeBytes);
                if retCode ~= this.alazarDefs.ApiSuccess; fprintf('DIGITIZER: Error: AlazarPostAsyncBuffer failed -- %s\n', errorToText(retCode)); status = false; return; end
            end
            
            this.timeTable.postBuffersToBoard = toc(this.timeTable.postBuffersToBoard);
        end
        
        function status = postSingleBufferToBoards(this, pbuffer)
            status = true;
            retCode = AlazarPostAsyncBuffer(this.system.boardHandle, pbuffer, this.vars.bufferSizeBytes);
            if retCode ~= this.alazarDefs.ApiSuccess; fprintf('DIGITIZER: Error: AlazarPostAsyncBuffer failed -- %s\n', errorToText(retCode)); status = false; end 
        end
        
        function status = armBoard(this)
            this.timeTable.armBoard = tic; 
            status = true;
            retCode = AlazarStartCapture(this.system.boardHandle);
            if retCode ~= this.alazarDefs.ApiSuccess; fprintf('DIGITIZER: Error: AlazarStartCapture failed -- %s\n', errorToText(retCode)); status = false; return; end
            
            this.timeTable.armBoard = toc(this.timeTable.armBoard);
        end
        
        function data   = getData(this)
           data = this.bufferDataOut;
        end
        
        function [dataOut, status] = acquire(this)
            this.timeTable.fullAcquisition = tic;
           
            %----- Configure AutoDMA -----
            status = this.setADMA(); if(~status); return; end
             
            % Post the buffers to the board
            status = this.postBuffersToBoard();  if(~status); return; end
            
            %------ Should Start Aquisition From This Point -----%
            % Arm the board system to wait for triggers
            status = this.armBoard(); if(~status); return; end
            
            % Init control variables
            this.vars.continueLooping = true;
            buffersCompleted = 0;
            captureDone = false;
            buf = 1;

            this.initPlots();
            
            this.timeTable.netAcquisition = tic;
            while ~captureDone && this.vars.continueLooping
                this.timeTable.loop(buf) = tic;
                bufferIndex = mod(buffersCompleted, this.system.bufferCount) + 1;
                pbuffer = this.buffers{1, bufferIndex};

                % Wait for the first available buffer to be filled by the board
                a = tic;
                [retCode, this.system.boardHandle, bufferOut] = AlazarWaitAsyncBufferComplete(this.system.boardHandle, pbuffer, 5000);
                this.timeTable.waitAsyncBufferComplete(buf) = toc(a);
                
                if retCode == this.alazarDefs.ApiSuccess
                    % This buffer is full
                    bufferFull = true;
                    captureDone = false;
                elseif retCode == this.alazarDefs.ApiWaitTimeout
                    % The wait timeout expired before this buffer was filled.
                    % The board may not be triggering, or the timeout period may be too short.
                    fprintf('DIGITIZER: Error: AlazarWaitAsyncBufferComplete timeout -- Verify trigger!\n');
                    bufferFull = false;
                    captureDone = true;
                else
                    % The acquisition failed
                    fprintf('DIGITIZER: Error: AlazarWaitAsyncBufferComplete failed -- %s\n', errorToText(retCode));
                    bufferFull = false;
                    captureDone = true;
                end
                
                if bufferFull
                    %cast Data and save it to RAM
                    a = tic;
                    setdatatype(bufferOut, 'uint16Ptr', 1, this.vars.samplesPerBufferAllCh);
                    this.timeTable.setDataType(buf) = toc(a);
                    
                    % Real-Time Post Processing (Only for NPT)
                    ppData = this.rtPostProcessing(bufferOut);
                    
                    if this.vars.drawData
                        this.drawData(ppData);
                    end
                    
                    a = tic;
                    if this.vars.exportData
                        this.exportData(buf, ppData);
                    end
                    this.timeTable.copyBufferData(buf) = toc(a);
                    
                    buf = buf +1;

                    % Make the buffer available to be filled again by the board
                    a = tic;
                    status = postSingleBufferToBoards(this, pbuffer); if(~status); captureDone = true; end
                    this.timeTable.postSingleBufferToBoards(buf) = toc(a);
                    
                    % Update progress
                    buffersCompleted = buffersCompleted + 1;
                    if buffersCompleted >= this.vars.numOfBuffers
                        captureDone = true;
                    end
                    
                end
                this.timeTable.loop(buf-1) = toc(this.timeTable.loop(buf-1));
            end
            this.timeTable.netAcquisition = toc(this.timeTable.netAcquisition);
            
            this.timeTable.abortRead = tic;
            AlazarAbortAsyncRead(this.system.boardHandle);
            this.timeTable.abortRead = toc(this.timeTable.abortRead);
            
            % Offline Post Processing (Only for TS)
            switch this.system.mode
                case 'TS'
                    if this.vars.useGPU
                        this.Data = gpuArray(this.bufferDataOut);
                    else
                        this.Data = this.bufferDataOut;
                    end
                    this.convertUnsignedSampleToVolts();
                    this.rawDataDeMultiplexing();
                    this.cropData();
                    dataOut = this.Data;
    
                    this.Data = [];
                case 'NPT'
                    if this.vars.exportData
                        dataOut = this.bufferDataOut(this.vars.triggerDelayManualSamples:this.NPT.postTriggerSamplesUser, :, :);
                    else
                        dataOut = [];
                    end
                case 'ContNPT'
                    if this.vars.useGPU
                        this.Data = gpuArray(this.bufferDataOut);
                    else
                        this.Data = this.bufferDataOut;
                    end
%                     this.bufferDataOut = [];
                    this.convertUnsignedSampleToVolts();
                    this.rawDataDeMultiplexing();
                    dataOut = this.Data;
                    
                    this.Data = [];
            end
                
            this.timeTable.fullAcquisition = toc(this.timeTable.fullAcquisition); 
        end

        function monitor(this, vars)
            this.vars.monitor = true;
            
            if nargin < 2
                vars.fs        = 100e6;
                vars.channels  = 1; 

                vars.triggerDelay = 0;
                vars.extClk  = false;
                vars.extTrig = true;

                vars.timeToSample = 2^14/vars.fs;
                vars.avgNum       = 100;
                
                vars.coupling = 1;
            end
            
            vars.mode      = 'NPT';
            vars.useGPU    = false;
            vars.numMeas   = inf;
            vars.draw      = true;
            vars.inputRange = '1V';
            
            this.setVars(vars);
            this.configure();

            this.acquire();
        end
        
        %Data Post Processing Functions
        function ppData = rtPostProcessing(this, buffer)
            switch this.system.mode
                case 'TS'
                    ppData = buffer.Value;
                case 'NPT'
                    ppData = buffer.Value;
                    
                    % Convert to Volts
                    shiftFact = 2;
                    bitsRange = 2^16/2;
                    ppData    = ((cast(ppData, 'single') - bitsRange) * this.system.voltsRange / (bitsRange * shiftFact));
                    
                    % Demultiplex
                    ppData = reshape(ppData, this.vars.channels, this.NPT.samplesPerRecord,  this.NPT.recordsPerBuffer);
                    
                    % Average
                    ppData = mean(ppData, 3);

                    % Organize dimensions
                    ppData = permute(ppData, [2,1]);

                    % Cut Redundant
                    ppData(this.NPT.cutIdxs, :) = [];

                case 'ContNPT'
                    ppData = buffer.Value;
            end
        end

        function exportData(this, idx, buffer)
            switch this.system.mode
                case {'TS', 'ContNPT'}
                    this.bufferDataOut(idx,:)   = buffer; 
                case 'NPT'
                    this.bufferDataOut(:,:,idx) = buffer;
            end
        end

        function convertUnsignedSampleToVolts(this)
            this.timeTable.convertUnsignedSampleToVolts = tic;
            
            shiftFact = 2;
            bitsRange = 2^16/2;
            this.Data = (cast(this.Data, 'single') - bitsRange) * this.system.voltsRange / (bitsRange * shiftFact);
            
            this.timeTable.convertUnsignedSampleToVolts = toc(this.timeTable.convertUnsignedSampleToVolts);
        end
        
        function rawDataDeMultiplexing(this)
            % bufferDataOut(input)  - [numOfBuffers x samplesPerBufferAllCh]
            % bufferDataOut(output) - [ch x samplesPerAcq]
            this.timeTable.rawDataDeMultiplexing = tic;
            switch this.system.mode
                case 'TS'
                    this.Data = reshape(reshape(this.Data', this.TS.samplesPerAcqAllCh, 1),...
                                        this.vars.channels,  this.TS.samplesPerAcq);
                case 'NPT'

                case 'ContNPT'
                    this.Data = reshape(this.Data, this.vars.numOfBuffers, this.vars.channels, this.vars.samplesPerBuffer);
                    this.Data(:,:,this.NPT.cutIdxs) = [];
                    this.Data = reshape(permute(this.Data, [2,3,1]), this.vars.channels, this.NPT.samplesPerAcqAllRecords);
                    this.Data(:, this.NPT.extraRecordIdxs) = [];
            end
            this.timeTable.rawDataDeMultiplexing = toc(this.timeTable.rawDataDeMultiplexing);             
        end
        
        function cropData(this)
            switch this.system.mode
                case 'TS'
                    if this.TS.exportCropped
                        this.croppedData = this.Data(:, 1:this.TS.croppedSamples);
                    end
                    this.Data    = this.Data(:,(this.TS.excSamplesPre+1):this.TS.excSamplesIdx); 
                case'NPT'
            end
        end
        
        function releaseBuffers(this)
            for bufferIndex = 1:size(this.buffers, 2)
                pbuffer =  this.buffers{1, bufferIndex};
                if ~isempty(pbuffer)
                    retCode = AlazarFreeBuffer(this.system.boardHandle, pbuffer);
                    if retCode ~= this.alazarDefs.ApiSuccess
                        fprintf('DIGITIZER: Error: AlazarFreeBuffer failed -- %s\n', errorToText(retCode));
                    end
                    clear pbuffer;
                end
            end
            this.system.allocatedBuffers = false;
        end

        %Figures Plots
        function initPlots(this)
            if this.vars.drawData
                channelChanged = this.vars.prevChannels ~= this.vars.channels;
                this.vars.prevChannles = this.vars.channels;
                if channelChanged || ~isgraphics(this.vars.figs.hFig) || sum(~isgraphics(this.vars.figs.hAx)) || sum(~isgraphics(this.vars.figs.hPlotSig)) || sum(~isgraphics(this.vars.figs.hPlotFFT))
                    if isgraphics(this.vars.figs.hFig)
                        close(this.vars.figs.hFig);
                    end
                    this.vars.figs.hFig  = figure();
                    this.vars.figs.hAx(1)   = subplot(1,2,1);
                    for i = 1:this.vars.channels
                        this.vars.figs.hPlotSig(i) = plot(this.vars.figs.hAx(1), this.vars.figs.tVec,zeros(1, this.NPT.samplesPerRecord)); hold on
                        legStr{i} = sprintf("ch-%d", i);
                    end
                    xlabel("Time [\mus]");
                    ylabel ("Voltage [mV]");
                    title("Signal")
                    
                    this.vars.figs.hAx(2)   = subplot(1,2,2);
                    for i = 1:this.vars.channels
                        this.vars.figs.hPlotFFT(i) = plot(this.vars.figs.hAx(2), this.vars.figs.fVec, zeros(1, this.NPT.samplesPerRecord)); hold on
                        legStr{i} = sprintf("ch-%d", i);
                    end
                    xlabel("Freq [MHz]");
                    ylabel ("Spectrum [AU]");
                    title("FFT")
                    
                end
            end
        end
        
        function drawData(this, ppData)
            try
                freqMat = abs(fftshift(fft(ppData, [], 1),1)).^2;
                for i=1:this.vars.channels
                    set(this.vars.figs.hPlotSig(i), 'YData', ppData(:,i));
                    set(this.vars.figs.hPlotFFT(i), 'YData', freqMat(:,i));
                end
                drawnow();
            catch
                this.vars.monitor = false;
                this.vars.continueLooping = false;
            end
        end
        
    end
end
