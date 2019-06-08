classdef fGen < handle
    
    properties
        inst
        me
        ChVars
        ChData
        valVars
    end
    
    methods (Static)       
        function chVars = getChVarsStruct()
            chVars.id              = [];
            chVars.triggerOwner    = [];
            chVars.daqFilter       = [];
            chVars.triggerOwner    = [];
            chVars.triggerWidth    = [];
            chVars.Sclk            = [];
            chVars.amp             = [];
            chVars.bias            = [];
            chVars.useExtSclkSrc   = [];
        end
        
        function chData = getChDataStruct
            chData.data            = [];
            chData.readyData       = [];
            chData.transmittedData = [];
            chData.dataLen         = [];
        end
        
        function uVars = uVarsCreate()
           uVars.ch{1} =   fGen.getChVarsStruct();
           uVars.ch{1}.id = 1;
           uVars.ch{2} =   fGen.getChVarsStruct();
           uVars.ch{2}.id = 2;
        end

    end
    
    methods
        function this = fGen()
            this.ChData{1} = this.getChDataStruct();
            this.ChData{2} = this.getChDataStruct();
            this.ChVars{1} = this.getChVarsStruct();
            this.ChVars{2} = this.getChVarsStruct();

            this.valVars = false;
        end
        
        function init(this)
            this.connect();
            this.reset()
        end
        
        function connect(this)
            this.inst = instrfind('Type', 'gpib', 'BoardIndex', 8, 'PrimaryAddress', 4, 'Tag', '');
            if isempty(this.inst)
                this.inst = gpib('AGILENT', 8, 4);
            else
                fclose(this.inst);
                this.inst = this.inst(1);
            end
            this.me.blocksize = 10000;
            set(this.inst, 'OutputBufferSize', this.me.blocksize);
            
            fopen(this.inst);
            if strcmp(this.inst.Status,'close')
                fprintf("Could not open fGen. please check device is ON and CONNECTED.")
            end
        end
        
        function reset(this)
            % Clearing all the old errors from error's buffer.
            fprintf(this.inst, '*CLS');
            % Reseting the device.
            fprintf(this.inst, '*RST');
            % Whoami
            this.me.idn = query(this.inst, '*IDN?');
            fprintf('Instrument IDN: %s', this.me.idn);
        end
        
        function setProperties(this, ch1, ch2)
           this.ChVars{1} = ch1;
           this.ChVars{2} = ch2;
        end
        
        function setData(this, ch1, ch2)
           this.ChData{1}.data    = ch1.data;
           this.ChData{1}.dataLen = ch1.dataLen;
           this.ChData{2}.data    = ch2.data;
           this.ChData{2}.dataLen = ch2.dataLen;
        end
        
        function configChannel(this, i)
            prepareSignal(this, i);
            err = zeros(1,4);
            while ~strcmp(err(1:4), '-100')
                % ------ Configure the Output ----------
                fprintf(this.inst, [':INSTrument ', num2str(this.ChVars{i}.id)]);
                if this.valVars; disp(query(this.inst, ':INST?')); end
                fprintf(this.inst, ':OUTP OFF');
                if this.valVars; disp(query(this.inst, ':OUTP?'));  end
                % Set a filter of 25Mhz at the output
                fprintf(this.inst, [':OUTP:FILT ', num2str(this.ChVars{i}.daqFilter), 'M']);
                if this.valVars; disp(query(this.inst, ':OUTP:FILT?')); end
                % ------ Configure the Sync Signal (output Trigger) ----------
                % Set a sync signal of a 16 samples bit
                fprintf(this.inst, ':OUTP:SYNC:SOUR BIT');
                if this.valVars; disp(query(this.inst, ':OUTP:SYNC:SOUR?')); end
                fprintf(this.inst, [':OUTP:SYNC:WIDT ',  this.ChVars{i}.triggerWidth]);
                if this.valVars; disp(query(this.inst, ':OUTP:SYNC:WIDT?')); end
                if (this.ChVars{i}.triggerOwner)
                    fprintf(this.inst, ':OUTP:SYNC ON');
                    if this.valVars; disp(query(this.inst, ':OUTP:SYNC?')); end
                end
%                 pause(0.1)
                
                % ----- Define the Samples CLK ---------
                % Set arbitrary waveform SCLK to 100 S/s
                sclkStr = num2str(this.ChVars{i}.Sclk, '%10.0e');
                sclkStr = strrep(sclkStr, '+', '');
                sclkStr = strrep(sclkStr, '-', '');
                fprintf(this.inst, [':FREQ:RAST ', sclkStr]);
                if this.valVars; disp(query(this.inst, ':FREQ:RAST?')); end
                % Select internal source for the SCLK
                if ~this.ChVars{i}.useExtSclkSrc
                    fprintf(this.inst, ':FREQ:RAST:SOUR INT');
                    if this.valVars; disp(query(this.inst, ':FREQ:RAST:SOUR?')); end
                else
                    % TODO: complete external sclk.
                end
                pause(0.1);
                
                % ----- Set Signal Properties ---------
                % Set signal source to as given by user
                fprintf(this.inst, ':FUNC:MODE USER');
                % Setting the amplitude and bias
                fprintf(this.inst, [':VOLT:LEV:AMPL ',  num2str(this.ChVars{i}.amp)] );
                if this.valVars; disp(query(this.inst, ':VOLT:LEV:AMPL? ')); end
                fprintf(this.inst, [':VOLT:LEV:OFFS ',  num2str(this.ChVars{i}.bias)]);
                if this.valVars; disp(query(this.inst, ':VOLT:LEV:OFFS?')); end
                pause(0.1);
                
                % ----- Set External Trigger Mode ---------
                % Set continous mode on (waveform is generated continously
                % and independently on external trigger)
                fprintf(this.inst, ':INIT:CONT 1');
                if this.valVars; disp(query(this.inst, ':INIT:CONT?')); end
                % Select channel 1 as the active channel
                fprintf(this.inst, ':INST:SEL CH1'); %TODO: check if necessary;
                if this.valVars; disp(query(this.inst, ':INST:SEL?')); end
                pause(0.1);
                
                % ------ Data Source ---------
                % Delete all segments from memory
                fprintf(this.inst, ':TRAC:DEL:ALL');
                % Define trace
                fprintf(this.inst, [':TRAC:DEF 1,' num2str(this.ChData{i}.dataLen)]);
                cLength = strcat(num2str(length(num2str(2*this.ChData{i}.dataLen))),num2str(2*this.ChData{i}.dataLen));
                samples = sprintf('#%s',cLength);
                % Select trace
                fprintf(this.inst, ':TRAC:SEL 1');
                % Write data blockwise
                fprintf(this.inst, '%s', [':TRAC:DATA' samples]);
                n = numel(this.ChData{i}.readyData);
                num_blocks = ceil(n/this.me.blocksize);
                pause(0.1);
                for k=1:num_blocks
                    pause(0.2);
                    fprintf('downloading chunk %d\n', k);
                    fwrite(this.inst, this.ChData{i}.readyData((k-1)*this.me.blocksize+1:min(n,k*this.me.blocksize)), 'uint8');
                end
                pause(0.1);
                % Connect the Output
                fprintf(this.inst, ':OUTP ON');
                if this.valVars; disp(query(this.inst, ':OUTP?')); end
                % Checking for Errors
                err = query(this.inst, ':SYST:ERR?');
                % Wait for completion
                query(this.inst, '*OPC?');
                pause(0.1);            
            end
            fprintf(this.inst, ':ENAB');
            pause(0.1);
            this.ChData{i}.transmittedData = this.ChData{i}.readyData;
        end
        
        function prepareSignal(this, i)
           % Convert data to 14 bit, as described in manual section 4-66
            data = uint16(((2^13-1)* this.ChData{i}.data)+2^13); % 14 bit data
            % Convert data to byte array [low byte 0, high byte 0, low byte 1, ...]
            bytes = zeros(2*length(data),1, 'uint8');
            % Low byte
            bytes(1:2:2*length(data)) = uint8(bitand(data, uint16(hex2dec('00FF'))));
            % High byte
            bytes(2:2:2*length(data)) = uint8(bitshift(bitand(data, uint16(hex2dec('FF00'))), -8));
            
            this.ChData{i}.readyData = bytes;
        end
        
        function disableOutput(this, ch)     
            fprintf(this.inst, [':INSTrument ', num2str(this.CH{ch}.id)]);
            fprintf(this.inst, ':OUTP OFF');
        end
        
        function disableAllOutputs(this)
            for ch = 1:2
                fprintf(this.inst, [':INSTrument ', num2str(this.CH{ch}.id)]);
                fprintf(this.inst, ':OUTP OFF');
            end
            fprintf(this.inst, ':OUTP:SYNC OFF');
        end
        
        function enableOutput(this, ch)     
            fprintf(this.inst, [':INSTrument ', num2str(this.CH{ch}.id)]);
            fprintf(this.inst, ':OUTP ON');
        end
        
        function enableAllOutputs(this)
            for ch = 1:2
                fprintf(this.inst, [':INSTrument ', num2str(this.CH{ch}.id)]);
                fprintf(this.inst, ':OUTP ON');
                if (this.CH{ch}.triggerOwner)
                    fprintf(this.inst, ':OUTP:SYNC ON');
%                     disp(query(this.inst, ':OUTP:SYNC?'))
                end
            end
        end
        
    end
end

