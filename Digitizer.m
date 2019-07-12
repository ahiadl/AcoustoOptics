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

        timeTable
        alazarDefs
    end
    
    methods (Static)
        
        function uVars = uVarsCreate() 
            uVars.mode                  = [];
            uVars.channels              = [];
            uVars.preTriggerSamples     = 14;
            uVars.voltsRange            = 2;
            uVars.useGPU                = false;
            uVars.fs                    = [];
            uVars.numOfBuffers          = [];
            uVars.bufferSizeBytes       = [];
            uVars.samplesPerAcq         = [];
            uVars.samplesPerBuffer      = [];
            uVars.samplesPerAcqAllCh    = [];
            uVars.samplesPerBufferAllCh = [];
        end     
        
    end
    
    methods
        
        function this = Digitizer()
            this.system.boardHandle  = [];
            this.system.systemId     = int32(1);
            this.system.boardId      = int32(1);
            this.system.channelsMask = 0;
            this.system.bufferCount  = 4;
            
            this.buffers = cell(1, this.system.bufferCount);
            
            this.initDefs();
        end
        
        function initDefs(this)
            AlazarDefs;
            names = who();
            for i = 1:length(names)
                if strcmp(names{i}, 'this'); continue; end
                this.alazarDefs.(names{i}) = eval(names{i});
            end
        end
        
        function setDigitizerVars(this, vars)
            this.vars.fs               = vars.fs;
            this.vars.channels         = vars.channels;
            this.vars.numOfBuffers     = vars.numOfBuffers;
            this.vars.bufferSizeBytes  = vars.bufferSizeBytes;
            this.vars.samplesPerAcq    = vars.samplesPerAcq;
            this.vars.samplesPerBuffer = vars.samplesPerBuffer;
            
            this.vars.samplesPerAcqAllCh    = vars.samplesPerAcqAllCh;
            this.vars.samplesPerBufferAllCh = vars.samplesPerBufferAllCh;
            
            this.vars.useGPU = vars.useGPU;
            
            this.system.triggerLevel_volts = 3; % external trigger level
            this.system.triggerRange_volts = 5; % external trigger input range 
            this.system.triggerLevel_code =(128 + 127 * this.system.triggerLevel_volts / this.system.triggerRange_volts);
            this.system.mode             = vars.mode;
            this.system.triggerDelay_sec = 0;
            this.system.triggerDelay_samples = uint32(floor(this.system.triggerDelay_sec * this.vars.fs + 0.5));
            this.system.triggerTimeout_sec = 0;
            this.system.triggerTimeout_clocks = uint32(floor(this.system.triggerDelay_sec / 10.e-6 + 0.5));
            this.system.voltsRange = vars.voltsRange;

            switch vars.mode
                case 'TS'
                    this.TS.samplesPerBuffer = vars.samplesPerBuffer;
                    this.TS.samplesPerBufferAllCh = vars.samplesPerBufferAllCh;
                    this.TS.numOfBuffers = vars.numOfBuffers;
                    this.CS = [];
                    this.NPT = [];
                case 'NPT'
                    this.TS = [];
                    this.CS = [];
                case 'CS'
                    this.TS = [];
                    this.NPT = [];
            end 
        end
        
        function initTimeTable(this)
            this.timeTable = struct();
        end
        
        function timeTable = getTimeTable(this)
            timeTable = this.timeTable;
        end
        
        function status = configure(this)
%             this.setDigitizerVars(vars);
            
            %--- Configure Clk & Trig ---
            status = this.confExtClk();  if(~status); return; end
            status = this.confChnls();   if(~status); return; end
            status = this.confExtTrig(); if(~status); return; end

            %--- Configure Dimensions ---
            status = setRecordSize(this); if(~status); return; end
            status = setADMA(this); if(~status); return; end

            % --- Create pointers to CPU memory Buffers ---
            status = allocateBuffers(this);
            this.allocateBufferOutMemory();
            
        end
        
        function allocateBufferOutMemory(this)
            this.timeTable.allocMemoryBuffers = tic;
            this.bufferDataOut = [];
            this.Data = [];
            this.uData = [];
            if this.vars.useGPU
                this.bufferDataOut = gpuArray(single(zeros(1, this.TS.samplesPerBufferAllCh)));
                this.Data          = gpuArray(single(zeros(this.vars.channels,  this.vars.samplesPerAcq)));
            else
                this.bufferDataOut = uint16(zeros(this.TS.numOfBuffers, this.TS.samplesPerBufferAllCh));
                this.uData         = single(zeros(this.vars.channels,  this.vars.samplesPerAcq));
                this.Data          = single(zeros(this.vars.channels,  this.vars.samplesPerAcq));
            end
            this.timeTable.allocMemoryBuffers = toc(this.timeTable.allocMemoryBuffers);  
        end
        
        function status = connect(this)
            this.timeTable.Connect = tic;
            
            status = true; 
%             AlazarDefs;
            if ~alazarLoadLibrary(); fprintf('Error: ATSApi library not loaded\n'); status = false; return; end
            this.system.boardHandle = AlazarGetBoardBySystemID(this.system.systemId, this.system.boardId);
            setdatatype(this.system.boardHandle, 'voidPtr', 1, 1);
            
            this.timeTable.Connect = toc(this.timeTable.Connect);
        end
            
        function status = confExtClk(this)
            this.timeTable.confExtClk = tic;
            
%             AlazarDefs;
            status = true;
            retCode = AlazarSetCaptureClock(this.system.boardHandle,  ... % HANDLE -- board handle
                                            this.alazarDefs.EXTERNAL_CLOCK,           ... % U32 -- clock source id
                                            this.alazarDefs.SAMPLE_RATE_USER_DEF,     ... % U32 -- sample rate id
                                            this.alazarDefs.CLOCK_EDGE_RISING,        ... % U32 -- clock edge id
                                            0);                          % U32 -- clock decimation                                
            if retCode ~= this.alazarDefs.ApiSuccess; fprintf('Error: AlazarSetCaptureClock failed -- %s\n', errorToText(retCode)); status = false; return; end
            
            this.timeTable.confExtClk = toc(this.timeTable.confExtClk);
        end
        
        function status = confChnls(this)
            this.timeTable.confChnls = tic;
%             AlazarDefs; 
            status = true;
            inputRange = this.alazarDefs.INPUT_RANGE_PM_1_V;
            coupling = this.alazarDefs.DC_COUPLING;
            impedance = this.alazarDefs.IMPEDANCE_50_OHM;
            successCode = this.alazarDefs.ApiSuccess;
            
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_A, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_B, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_C, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_D, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_E, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_F, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_G, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_H, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_I, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_J, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_K, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_L, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_M, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_N, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_O, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
            retCode = AlazarInputControlEx( this.system.boardHandle, this.alazarDefs.CHANNEL_P, coupling, inputRange, impedance);
            if (retCode ~= successCode); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end

            this.system.channelMask = 0;
            for i = 1:this.vars.channels
                this.system.channelMask = bitor(bitshift(this.system.channelMask,1),1);
            end
            
            this.timeTable.confChnls = toc(this.timeTable.confChnls);
        end
        
        function status = confExtTrig(this)
            this.timeTable.confExtTrig = tic;
            
%             AlazarDefs; 
            status = true;
            % Disable internal Triggers and set External Trigger
            retCode = AlazarSetTriggerOperation(this.system.boardHandle, this.alazarDefs.TRIG_ENGINE_OP_J,...
                                                this.alazarDefs.TRIG_ENGINE_J, this.alazarDefs.TRIG_EXTERNAL,...
                                                this.alazarDefs.TRIGGER_SLOPE_POSITIVE, this.system.triggerLevel_code,... %
                                                this.alazarDefs.TRIG_ENGINE_K, this.alazarDefs.TRIG_DISABLE, this.alazarDefs.TRIGGER_SLOPE_POSITIVE, 128 );                                                    
            if retCode ~= this.alazarDefs.ApiSuccess ; fprintf('Error: AlazarSetTriggerOperation failed -- %s\n', errorToText(retCode)); status = false; return; end

            % Set External trigger Configurations
            retCode = AlazarSetExternalTrigger( this.system.boardHandle, this.alazarDefs.DC_COUPLING, this.alazarDefs.ETR_TTL );
            if retCode ~= this.alazarDefs.ApiSuccess; fprintf('Error: AlazarSetExternalTrigger failed -- %s\n', errorToText(retCode)); status = false; return; end

            % Trigger Delay set to 0
            retCode = AlazarSetTriggerDelay(this.system.boardHandle, this.system.triggerDelay_samples);
            if retCode ~= this.alazarDefs.ApiSuccess; fprintf('Error: AlazarSetTriggerDelay failed -- %s\n', errorToText(retCode)); status = false; return; end

            % Set Trigger Timeout 
            retCode = AlazarSetTriggerTimeOut(this.system.boardHandle, this.system.triggerTimeout_clocks); % U32 -- timeout_sec / 10.e-6 (0 == wait forever)       
            if retCode ~= this.alazarDefs.ApiSuccess; fprintf('Error: AlazarSetTriggerTimeOut failed -- %s\n', errorToText(retCode)); status = false; return; end    
            
            this.timeTable.confExtTrig = toc(this.timeTable.confExtTrig);
        end
        
        function status = setADMA(this)
            this.timeTable.setADMA = tic;
            
%             AlazarDefs; 
            status = true;
            if (strcmp(this.system.mode, 'TS')) 
                admaFlags = this.alazarDefs.ADMA_EXTERNAL_STARTCAPTURE + this.alazarDefs.ADMA_TRIGGERED_STREAMING;
                retCode = AlazarBeforeAsyncRead(this.system.boardHandle, this.system.channelMask, 0, this.TS.samplesPerBuffer, 1, hex2dec('7FFFFFFF'), admaFlags);

            elseif (strcmp(this.system.mode, 'NPT'))
                 admaFlags = this.alazarDefs.ADMA_EXTERNAL_STARTCAPTURE + this.alazarDefs.ADMA_NPT; % Select AutoDMA flags as required
                 retCode = AlazarBeforeAsyncRead(this.system.boardHandle, this.system.channelMask, -int32(this.NPT.preTriggerSamples), this.NPT.samplesPerRecord, this.NPT.recordsPerBuffer, hex2dec('7FFFFFFF'), admaFlags);
            end
            if retCode ~= this.alazarDefs.ApiSuccess; fprintf('Error: AlazarBeforeAsyncRead failed -- %s\n', errorToText(retCode)); status = false; return; end
            
            this.timeTable.setADMA = toc(this.timeTable.setADMA);
        end
        
        function status = allocateBuffers(this)
            this.timeTable.allocateBuffers = tic;
            
%             AlazarDefs; 
            status = true;
            for j = 1 : this.system.bufferCount
                pbuffer = AlazarAllocBuffer(this.system.boardHandle, this.vars.bufferSizeBytes);
                if pbuffer == 0; fprintf('Error: AlazarAllocBuffer %u samples failed\n', SamplingCard.samplesPerBuffer); status = false; return; end
                this.buffers(1, j) = { pbuffer };
            end
            
            this.timeTable.allocateBuffers = toc(this.timeTable.allocateBuffers);
        end
        
        function status = setRecordSize(this)
%             AlazarDefs; 
            status = true;
            if strcmp(this.system.mode, 'NPT')
                retCode = AlazarSetRecordSize(this.system.boardHandle, this.NPT.preTriggerSamples, this.NPT.postTriggerSamples);
                if retCode ~= this.alazarDefs.ApiSuccess; fprintf('Error: AlazarSetRecordSize failed -- %s\n', errorToText(retCode)); status = false; return; end
            end
        end
        
        function [bufferDataOut, status] = acquire(this)
            switch this.system.mode
                case 'TS'
                    [bufferDataOut, status] = acquireDataTSdig();
                case 'NPT' 
                    [bufferDataOut, status] = acquireDataNPTdig();
                case 'CS'
                    bufferDataOut = [];
                    status = false;
            end
        end
        
        function status = postBuffersToBoard(this)
            this.timeTable.postBuffersToBoard = tic;
            
%             AlazarDefs; 
            
            status = true;
            for bufferIndex = 1 : this.system.bufferCount
                pbuffer = this.buffers{1, bufferIndex};
                retCode = AlazarPostAsyncBuffer(this.system.boardHandle, pbuffer, this.vars.bufferSizeBytes);
                if retCode ~= this.alazarDefs.ApiSuccess; fprintf('Error: AlazarPostAsyncBuffer failed -- %s\n', errorToText(retCode)); status = false; return; end
            end
            
            this.timeTable.postBuffersToBoard = toc(this.timeTable.postBuffersToBoard);
        end
        
        function status = postSingleBufferToBoards(this, pbuffer)
%             AlazarDefs; 
            
            status = true;
            retCode = AlazarPostAsyncBuffer(this.system.boardHandle, pbuffer, this.vars.bufferSizeBytes);
            if retCode ~= this.alazarDefs.ApiSuccess; fprintf('Error: AlazarPostAsyncBuffer failed -- %s\n', errorToText(retCode)); status = false; end 
        end
        
        function status = armBoard(this)
            this.timeTable.armBoard = tic;
%             AlazarDefs; 
            status = true;
            % startTime = tic;
            retCode = AlazarStartCapture(this.system.boardHandle);
            if retCode ~= this.alazarDefs.ApiSuccess; fprintf('Error: AlazarStartCapture failed -- %s\n', errorToText(retCode)); status = false; return; end
            % startTime = toc(startTime);
            
            this.timeTable.armBoard = toc(this.timeTable.armBoard);
        end
        
        function data = getData(this)
           data = this.bufferDataOut;
%            this.bufferDataOut = [];
        end
        
        function [dataOut, status] = acquireDataTS(this)
            this.timeTable.fullAcquisition = tic;
           
            %----- Configure AutoDMA -----
            status = this.setADMA(); if(~status); return; end
             
            % Post the buffers to the board
            status = this.postBuffersToBoard();  if(~status); return; end
            
            %------ Should Start Aquisition From This Point -----%
            % Arm the board system to wait for triggers
            status = this.armBoard(); if(~status); return; end
            
            % Init control variables
            buffersCompleted = 0;
            captureDone = false;
            buf = 1;
            
            this.timeTable.netAcquisition = tic;
            while ~captureDone
                bufferIndex = mod(buffersCompleted, this.system.bufferCount) + 1;
                pbuffer = this.buffers{1, bufferIndex};

                % Wait for the first available buffer to be filled by the board
                [retCode, this.system.boardHandle, bufferOut] = AlazarWaitAsyncBufferComplete(this.system.boardHandle, pbuffer, 5000);
                if retCode == this.alazarDefs.ApiSuccess
                    % This buffer is full
                    bufferFull = true;
                    captureDone = false;
                elseif retCode == this.alazarDefs.ApiWaitTimeout
                    % The wait timeout expired before this buffer was filled.
                    % The board may not be triggering, or the timeout period may be too short.
                    fprintf('Error: AlazarWaitAsyncBufferComplete timeout -- Verify trigger!\n');
                    bufferFull = false;
                    captureDone = true;
                else
                    % The acquisition failed
                    fprintf('Error: AlazarWaitAsyncBufferComplete failed -- %s\n', errorToText(retCode));
                    bufferFull = false;
                    captureDone = true;
                end
                
                if bufferFull
                    %cast Data and save it to RAM
                    setdatatype(bufferOut, 'uint16Ptr', 1, this.TS.samplesPerBufferAllCh);
                    if this.vars.useGPU
                        this.bufferDataOut(:) = cast(bufferOut.Value, 'single'); %TODO: check that curData is still gpuArray 
                        this.convertUnsignedSampleToVoltsGPU();
                        this.rawDataDeMultiplexingGPU(buf);
                    else
                        this.bufferDataOut(buf,:) = bufferOut.Value;
                    end
                    buf = buf +1;

                    a = tic;
                    % Make the buffer available to be filled again by the board
                    status = postSingleBufferToBoards(this, pbuffer); if(~status); captureDone = true; end
                    
                    this.timeTable.postSingleBufferToBoards(buf) = toc(a);
                    
                    % Update progress
                    buffersCompleted = buffersCompleted + 1;
                    if buffersCompleted >= this.vars.numOfBuffers
                        captureDone = true;
                    end
                end
            end
            this.timeTable.netAcquisition = toc(this.timeTable.netAcquisition); 
            if ~this.vars.useGPU
                this.convertUnsignedSampleToVolts();
                this.rawDataDeMultiplexing();
            end
            dataOut = this.Data;
            this.timeTable.fullAcquisition = toc(this.timeTable.fullAcquisition); 
        end

        function convertUnsignedSampleToVolts(this)
            this.timeTable.convertUnsignedSampleToVolts = tic;
            
            shiftFact = 2;
            bitsRange = 2^16/2;

            this.uData = (cast(this.bufferDataOut, 'single') - bitsRange) * this.system.voltsRange / (bitsRange * shiftFact);
            
            this.timeTable.convertUnsignedSampleToVolts = toc(this.timeTable.convertUnsignedSampleToVolts);
        end
        
        function rawDataDeMultiplexing(this)
            % bufferDataOut(input)  - [numOfBuffers x samplesPerBufferAllCh]
            % bufferDataOut(output) - [ch x samplesPerAcq]
            this.timeTable.rawDataDeMultiplexing = tic;

            this.Data = reshape(reshape(this.uData', this.vars.samplesPerAcqAllCh, 1),...
                                      this.vars.channels,  this.vars.samplesPerAcq);

            this.timeTable.rawDataDeMultiplexing = toc(this.timeTable.rawDataDeMultiplexing);             
        end

        function convertUnsignedSampleToVoltsGPU(this)
            this.timeTable.convertUnsignedSampleToVolts = tic;
            
            shiftFact = 2;
            bitsRange = 2^16/2;

            this.bufferDataOut = (this.bufferDataOut - bitsRange) * this.system.voltsRange / (bitsRange * shiftFact);
            
            this.timeTable.convertUnsignedSampleToVolts = toc(this.timeTable.convertUnsignedSampleToVolts);
        end
        
        function rawDataDeMultiplexingGPU(this, buf)
            this.timeTable.rawDataDeMultiplexing = tic;
            
            this.Data(:, ((buf-1)*this.vars.samplesPerBuffer+1) : buf*this.vars.samplesPerBuffer) =...
                            reshape(this.bufferDataOut', this.vars.channels,  this.vars.samplesPerBuffer);
            
            this.timeTable.rawDataDeMultiplexing = toc(this.timeTable.rawDataDeMultiplexing);             
        end
        
        function releaseBuffers(this)
            for bufferIndex = 1:this.system.bufferCount
                pbuffer =  this.buffers{1, bufferIndex};
                retCode = AlazarFreeBuffer(this.system.boardHandle, pbuffer);
                if retCode ~= this.alazarDefs.ApiSuccess
                    fprintf('Error: AlazarFreeBuffer failed -- %s\n', errorToText(retCode));
                end
                clear pbuffer;
            end
        end
    end
end

