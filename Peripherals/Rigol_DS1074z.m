classdef Rigol_DS1074z < handle
    %Class object fo the Rigol_DS1054z Scope 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scope automation: Rigol_DS1054z Subclass
% =======================================================================
% - Writen by Ofry Livney 2019.08
% =========================================================================
% GENERAL NOTES: helpful insights for scope automation
% - 
% - IMPORTANT NOTE: for further elaboration and use for specific scopes in
% - the lab make sure you are using the correct programmer's guide
% - i.e. for the KEYSIGHT DSOX1102G the Programmer's Guide can be found at
% - https://www.keysight.com/upload/cmc_upload/All/1000_X-Series_prog_guide.pdf
%
% - To find the compatible VISA address port number use the company's connection 
% - expert. i.e. for the KEYSIGHT DSOX1102G use the Keysight connection
% - expert program which can be found at 
% - https://www.keysight.com/en/pd-1985909/io-libraries-suite?nid=-33330.977662.00&cc=IL&lc=eng&cmpid=zzfindiolib
% - Keysight Connection Expert is enough to get the VISA address port
% - number for non-Keysight instruments as well.
% - 
% - For more info about the visaobj options type instrhelp visa in the
% - Command Window. 
% - 
% HOW TO USE THIS AUTOMATION:
% - Step 1: Setting the scope's default properties:
% - all scope parameters (channels,trigger,timebase,waveform)
% - should be set in the Constructor method. This determines the default
% - properties of our scope.
% 
% - Step 2: initiallize the scope instance: 
% - i.e. myScope = Rigol_DS1054z ('ni','USB0::0x1AB1::0x04CE::DS1ZA181706779::0::INSTR') 
% 
% - Step 3: Use the myScope.SetupScope() method to communicate all scope object
% - properties to the real scope. As can be seen in the implementation, 
% - this method calls all configuring methods.   
% 
% - Step 4: Measuring the scope's input:
% - simply use myScope.TakeMeasurement(NumofCh) method by specifing the input 
% - channel in the argument. Data will be stored in the myScope.Data 
% - property (which is defined as a property for the superclass). 
% - important note: using this method again will run over the last  
% - measurement, so save the data in another variable before repeating 
% - this step.
% 
% - step 5: Changing scope settings after initiallization:
% - i.e. if you would like to change the triger level - just change the
% - value of the object's property - this.Trigger.trigLevel= 0.5; 
% - and afterwards you should use the myScope.ConfTrigger() method to 
% - communicate the new settings to the scope. It is possible to use the 
% - myScope.SetupScope() again, but that will communicate all current
% - configurations and that isn't needed for a specific change. In the same
% - manner one could update the settings for Waveform/Timebase/specific
% - Channel, with the corresponding method.
% 
% - step 6: Adding functionalities to the automation:
% - this can be done by updating the methods for specific updates of the 
% - scope's active configurations: Channels(ConfChannel(Ch)),
% - Trigger(ConfTrigger()),Timebase(ConfTimebase()),Waveform(ConfWaveform()).
% - i.e. I'm using the KEYSIGHT DSOX1102G for my experiment and I would  
% - like to add a setting for the trigger of the scope. I should find the
% - relevant method - ConfTrigger() and update it according to the
% - programmer's guide.for example, I wanted to add a control of the 
% - TRIGger:SWEap mode(can recieve two optional values: {AUTO|NORMal}). 
% - I would add the following line to ConfTrigger():
% - fprintf(this.VisaObj, [':TRIGger:SWEap, trigger.sweapmode)]);
% - Then, all I have to do is add my setting (sweapmode) to the Constructor
% - method and the next time I initiallize an instance of a
% - Keysight DSOX 1102G scope -> it'll have the sweapmode property and the  
% - trigger will be configured accordingly for my experiment.
% =========================================================================
% SPECIFIC NOTES: Rigol DS1054z scope automation
% =========================================================================  
% SCOPE PROPERTIES & SETTINGS: RIGOL DS1054z scope automation
% - NOTE: Full explanation about the commands&settings of this scope is 
% - available in the relevant programmer's guide.
%
% - In case you are initiallizing a RIGOL DS1054z their fields
% - should be as such:
% - Channels prop: its fields are Ch1,Ch2,Ch3,Ch4.
% - and for each channel struct the fields are: 
% - 1) status,id,coup,invert,offset,scale - which should be given
% - for proper initiallization of the scope setup.
% - 2) v0,v_ref,dv - are added during the SetupScope method or when directly
% - configuring a channel with the ConfChannel method.
%
% - Trigger prop: the fields of this struct: mode, trigCh, edge, trigLevel
% 
% - Timebase prop: the fields of this struct: mode, scale, offset
%
% - Waveform prop: the fields of this struct: chToSample, WavPointsMode, Format  
% =========================================================================
    properties (SetAccess = public)
        %initiallization is implemented in each  class constructor
        VisaObj; % The c'tor creates a visa object for matlab use -
         % a controller obj which recieves the commands using fprintf,fread,query functions 
         % i.e fprintf(visaobj,'WAV:POINTS 1000000'), query(visaobj,'WAV:POINTS?')  
        Vendor; %Company that created the scope. i.e. 'Keysight'
        VisaAddr; % prop for communication intiallization.
%         rawData;
%         data;     % holds the captured waveform data from the last measurement
        Channels; % prop holds settings for all needed channels
        Trigger;  % prop holds all the settings for scope's trigger.         
        Timebase; % prop holds the settings of the scope's screen Timebase                  
        Waveform; % prop holds the waveform settings for the planned data aquisition                                            
        Acq;
        res;
        
        connected
    end
    
    methods
        function this = Rigol_DS1074z (VendorName, VisaAddress)
            %Constructor
            
            %Interface configuration and instrument connection
            this.VisaObj = visadev(VendorName,VisaAddress);%initiallizes the visa object
            this.VisaObj.InputBufferSize = 5000;  % Set the buffer size at 5e7 - note: it is not verified if this setting this input buffer size is actually needed
            this.VisaObj.Timeout = 5;   % Set the timeout value
            this.VisaObj.ByteOrder = 'littleEndian'; % Set the Byte order
            this.Vendor = VendorName;
            this.VisaAddr = VisaAddress; 
            this.OpenConnection();  
            
            %initiallizng default settings for class properties
            %Channels property fields
            this.Channels.NumOfChannels=4; %the keysight DSOX1102G has 2 channels
            %CH1
            this.Channels.Ch1.id= 1;
            this.Channels.Ch1.status= 'OFF';
            this.Channels.Ch1.coup= 'AC';
            this.Channels.Ch1.invert= 'OFF';
            this.Channels.Ch1.offset= 0;
            this.Channels.Ch1.scale= 2;

            %CH2
            this.Channels.Ch2.id= 2;
            this.Channels.Ch2.status= 'OFF';
            this.Channels.Ch2.coup= 'AC';
            this.Channels.Ch2.invert= 'OFF';
            this.Channels.Ch2.offset= -10;
            this.Channels.Ch2.scale= 5;

            %CH3
            this.Channels.Ch3.status = 'OFF';
            this.Channels.Ch3.id = 3;
            this.Channels.Ch3.probe = 1;
            this.Channels.Ch3.coup = 'AC';
            this.Channels.Ch3.offset = 0;
            this.Channels.Ch3.scale = 0.05;
            this.Channels.Ch3.invert = 'OFF';
            
            %CH4
            this.Channels.Ch4.status = 'OFF';
            this.Channels.Ch4.id = 4;
            this.Channels.Ch4.probe = 1;
            this.Channels.Ch4.coup = 'AC';
            this.Channels.Ch4.offset = 0;
            this.Channels.Ch4.scale = 1;
            this.Channels.Ch4.invert = 'OFF';
             
            %Trigger property fields
            this.Trigger.mode = 'EDGE';
            this.Trigger.trigCh= 1;
            this.Trigger.edge = 'POSitive';
            this.Trigger.trigLevel= 0; %[v]
            
            %Timebase property fields
            this.Timebase.mode= 'MAIN';
            this.Timebase.scale= 100e-6;
            this.Timebase.offset= 0;
 
            %Waveform property fields
            this.Waveform.chToSample= 1;
            this.Waveform.WavPointsMode= 'NORMal';
            this.Waveform.Format='BYTE';
             
        end  
        
        function this = SetupScope(this)
            %Sets all scope configurations at the same instance.
            %fprintf(this.VisaObj,':*RST; :AUTOSCALE'); % Reset the instrument and autoscale and stop
            %channel configurations 
            %DS1054z has 4 channels
            this = this.SetupChannel(1);
            this = this.SetupChannel(2);
            this = this.SetupChannel(3);
            this = this.SetupChannel(4);
            %Trigger Configurations
            this.ConfTrigger();
            %TimeBase Configurations
            this.ConfTimebase();
            %WaveForm Configurations
            this.ConfWaveform();
        end
        
        function this = TakeMeasurement(this, Ch)
                visaObj = this.VisaObj;

                %     S1. :STOP Set the instrument to STOP state (you can only read
                %     the waveform data in the internal memory when the
                %     oscilloscope is in STOP state)
                    fprintf(visaObj, ':STOP');
                %     S2. :WAV:SOUR CHAN1 Set the channel source to CH1
                    fprintf(visaObj, [':WAV:SOUR CHAN', num2str(Ch)]);
                %     S3. :WAV:MODE RAW Set the waveform reading mode to RAW
                    fprintf(visaObj, ':WAV:MODE RAW');
                %    S4. :WAV:FORM BYTE Set the return format of the waveform data to BYTE
                    fprintf(visaObj, ':WAV:FORM BYTE');
                %    S5. Acquire data and convert it volts    
                    fprintf(visaObj, 'ACQuire:MDEPth?' );      
                    pause(0.2);
                    memDepth=char(fread(visaObj)');
                    if (strcmp(memDepth(1:4),'AUTO'))
                        numGrids=12;
                        fprintf(visaObj, ':TIMebase:MAIN:SCALe?');
                        pause(0.2);
                        dt=str2double(char(fread(visaObj)'));
                        fprintf(visaObj, ':ACQuire:SRATe?');
                        pause(0.2);
                        SampleRate=str2double(char(fread(visaObj)'));
                        res.totalNumOfSamples = SampleRate*dt*numGrids;   
                    else
                        res.totalNumOfSamples = str2double(memDepth);
                    end
                    samplesInBuffer = 250000;
                    
                    numOfBuffers = ceil(res.totalNumOfSamples / 250000);
                    data = zeros(samplesInBuffer+12, numOfBuffers-2);%???????????12?????
                    lastBuffer = [];

                    for i = 1:numOfBuffers
                        startIdx = (i-1)*samplesInBuffer + 1;
                        endIdx   = i*samplesInBuffer;
                        samplesPerBuffer = endIdx-(startIdx-1);
                        sizeofByte = 8;
                        bytesPerBuffer=samplesPerBuffer*sizeofByte;
                        fprintf(visaObj, [':WAV:STAR ', num2str(startIdx)]);
                        fprintf(visaObj, [':WAV:STOP ', num2str(endIdx)]);
                        fprintf(visaObj, ':WAV:DATA?');
                        if i == 1
                            [firstBuffer, ~] = fread(visaObj, bytesPerBuffer); 
                            firstBuffer = firstBuffer(12:length(firstBuffer)-1);
                        elseif i == numOfBuffers
                            [lastBuffer, ~] = fread(visaObj, bytesPerBuffer);
                            lastBuffer = lastBuffer(12:length(lastBuffer)-1);
                        else
                            [data(:,i-1), ~] = fread(visaObj, bytesPerBuffer);   
                        end
                    end
                    data=data(12:length(data)-1,:);
                    rawData = [firstBuffer(:); data(:); lastBuffer(:)];

                    fprintf(visaObj, ':WAVeform:YORigin?' );      
                    res.v0    = str2num(char(fread(visaObj)'));
                    fprintf(visaObj, ':WAVeform:YREFerence?' );   
                    res.v_ref = str2num(char(fread(visaObj)'));
                    fprintf(visaObj, ':WAVeform:YINCrement?' );   
                    res.dv    = str2num(char(fread(visaObj))'); 

                    fprintf(visaObj, ':WAVeform:XORigin?' );      
                    res.to    = str2num(char(fread(visaObj)'));
                    fprintf(visaObj, ':WAVeform:XREFerence?' );      
                    res.tRef  = str2num(char(fread(visaObj)'));
                    fprintf(visaObj, ':WAVeform:XINCrement?' );      
                    res.dt   = str2num(char(fread(visaObj)'));

                    res.WaveVolts = (rawData-res.v0-res.v_ref).*res.dv;   % Convertion of Samples to Volts
                    res.timeVec = (0:1:(length(res.WaveVolts)-1))*res.dt + res.to;

                    fprintf(visaObj, ':RUN');
                % Storing waveform data for post-processing and anallizing tasks
                    this.Data.RawData = rawData;
                    this.Data.CapturedWaveform = res;
        end


        function res = acquire(this)
            chStr = ['Ch', num2str(this.Waveform.chToSample)];
            ch = this.Channels.(chStr);
            
%             this.channelOn(ch.id);
%             fprintf(visaObj, ':ACQuire:SRATe?');
%             res.fs = str2double(char(fread(visaObj)'));

            fprintf(this.VisaObj, ':DISPlay:CLEar' );              % Clear scope display
            pause(this.Acq.waitAfterClear);                                 % Data aquisition time
%             pause(2);
            
            fprintf(this.VisaObj, ':WAV:DATA?' );               % Read the waveform data
            pause(0.2);
            [rawData,len] = fread(this.VisaObj, 2048);          % Request the data
            WaveSamples   = rawData(12:len-1);                  % Sampled data
            data          = ((WaveSamples-ch.v0-ch.v_ref).*ch.dv).';   % Convertion of Samples to Volts
            
            res.rawData = rawData;
            res.data = data;
            
            res.tVec = this.Timebase.dispTimeVec;
            res.t0   = this.Timebase.t0;
            res.tRef = this.Timebase.tRef;
            res.dt   = this.Timebase.dt;
            
            res.dv   = ch.dv;
            res.v0   = ch.v0;
            res.vRef = ch.v_ref;
            this.res = res;
        end

        function channelOn(this, id)
            for i=1:4
                chStr = ['Ch', num2str(id)];
                if i == id
                    status = 'ON';
                else
                    status = 'OFF';
                end
                this.Channels.(chStr).status = status;
                fprintf(this.VisaObj, [':CHANnel',num2str(i),':DISPlay ', status]);
            end
        end
        
        function this = SetupChannel(this, NumofCh)
             %this method calls upon ConfChannel method to configure the
             %given NumofCh settings.
             switch NumofCh
                case 1
                    this = ConfChannel(this, this.Channels.Ch1);
                case 2
                    this = ConfChannel(this, this.Channels.Ch2);
                case 3
                    this = ConfChannel(this, this.Channels.Ch3);
                case 4
                    this = ConfChannel(this, this.Channels.Ch4);
                otherwise
                    disp('Error in SetupChannel method: Channel number does not exist!');
             end
        end 
        
        function setChannelVars(this, ch)
%             this.Channels.(['Ch',num2str(ch.id)]).status = ch.status;
            this.Channels.(['Ch',num2str(ch.id)]).coup    = ch.coup;
            this.Channels.(['Ch',num2str(ch.id)]).invert  = ch.invert;
            this.Channels.(['Ch',num2str(ch.id)]).offset  = ch.offset;
            this.Channels.(['Ch',num2str(ch.id)]).fineRes = ch.fineRes;
            this.Channels.(['Ch',num2str(ch.id)]).scale   = ch.scale;
        end
        
        function vars = getVars(this)
           vars.acq      = this.Acq;
           vars.timebase = this.Timebase;
           vars.waveform = this.Waveform;
           vars.trigger  = this.Trigger;
           vars.channels = this.Channels;
        end
        
        function setGlobalConf(this, conf)
            this.Trigger.mode      = conf.trigger.mode;
            this.Trigger.trigCh    = conf.trigger.trigCh;
            this.Trigger.edge      = conf.trigger.edge;
            this.Trigger.trigLevel = conf.trigger.trigLevel;
            
            this.Timebase.mode   = conf.timebase.mode;
            this.Timebase.scale  = conf.timebase.scale;
            this.Timebase.offset = conf.timebase.offset;
            
            this.Waveform.chToSample    = conf.waveform.chToSample;
            this.Waveform.WavPointsMode = conf.waveform.WavPointsMode;
            this.Waveform.Format        = conf.waveform.Format;
            
            this.Acq.mode           = conf.acq.mode;
            this.Acq.numOfAvg       = conf.acq.numOfAvg;
            this.Acq.waitAfterClear = conf.acq.waitAfterClear;
        end
        
        function this = ConfChannel(this, Ch)
            %Method for Configuring all the settings of the given channel
            %Ch at the same instance. It is possible to change these
            %settings separately by accessing them directly.
            fprintf(this.VisaObj, [':CHANnel',num2str(Ch.id),':DISPlay ', Ch.status]);
            if (strcmp(Ch.status, 'ON'))
                fprintf(this.VisaObj, [':CHANnel',num2str(Ch.id),':COUPling ', Ch.coup]);
                fprintf(this.VisaObj, [':CHANnel',num2str(Ch.id),':INVert ', Ch.invert]);
                fprintf(this.VisaObj, [':CHANnel',num2str(Ch.id),':OFFSet ', num2str(Ch.offset)]);
                if Ch.fineRes
                    fprintf(this.VisaObj, [':CHANnel',num2str(Ch.id),':VERNier ON']);
                else
                    fprintf(this.VisaObj, [':CHANnel',num2str(Ch.id),':VERNier OFF']);
                end
                fprintf(this.VisaObj, [':CHANnel',num2str(Ch.id),':SCALe ', num2str(Ch.scale)]);
                fprintf(this.VisaObj, [':CHANnel',num2str(Ch.id),':UNITs VOLT']); %{VOLT|AMP}
            end

            fprintf(this.VisaObj, ':WAV:SOUR?');
            fs = char(fread(this.VisaObj)');
            OrigCH    = str2double(fs);
            fprintf(this.VisaObj, [':WAV:SOUR CHAN', num2str(Ch.id)]);
            fprintf(this.VisaObj, ':WAVeform:YORigin?' );      
            Ch.v0    = str2double(char(fread(this.VisaObj)')); 
            fprintf(this.VisaObj, ':WAVeform:YREFerence?' );   
            Ch.v_ref = str2double(char(fread(this.VisaObj)'));  
            fprintf(this.VisaObj, ':WAVeform:YINCrement?' );   
            Ch.dv    = str2double(char(fread(this.VisaObj)')); 
            fprintf(this.VisaObj ,[':WAV:SOUR CHAN', num2str(OrigCH)]);
        end

        function this = ConfChannelFromVars(this, id)
            %Method for Configuring all the settings of the given channel
            %Ch at the same instance. It is possible to change these
            %settings separately by accessing them directly.
            chStr = ['Ch', num2str(id)];
            Ch = this.Channels.(chStr);
            
            fprintf(this.VisaObj, [':CHANnel',num2str(this.Channels.(chStr).id),':DISPlay ', 'ON']);
            fprintf(this.VisaObj, [':CHANnel',num2str(this.Channels.(chStr).id),':COUPling ', this.Channels.(chStr).coup]);
            fprintf(this.VisaObj, [':CHANnel',num2str(this.Channels.(chStr).id),':INVert ', this.Channels.(chStr).invert]);
            fprintf(this.VisaObj, [':CHANnel',num2str(this.Channels.(chStr).id),':OFFSet ', num2str(this.Channels.(chStr).offset)]);
            if this.Channels.(chStr).fineRes
                fprintf(this.VisaObj, [':CHANnel',num2str(this.Channels.(chStr).id),':VERNier ON']);
            else
                fprintf(this.VisaObj, [':CHANnel',num2str(this.Channels.(chStr).id),':VERNier OFF']);
            end
            fprintf(this.VisaObj, [':CHANnel',num2str(this.Channels.(chStr).id),':SCALe ', num2str(this.Channels.(chStr).scale)]);
            fprintf(this.VisaObj, [':CHANnel',num2str(this.Channels.(chStr).id),':UNITs VOLT']); %{VOLT|AMP}

            fprintf(this.VisaObj, ':WAV:SOUR?');
            OrigCH = str2double(char(fread(this.VisaObj)'));

            fprintf(this.VisaObj, [':WAV:SOUR CHAN', num2str(Ch.id)]);
%             fprintf(this.VisaObj, ':WAV:SOUR?');
%             pause(0.2);
            
            fprintf(this.VisaObj, ':WAVeform:YORigin?' );   
%             pause(0.5);
            tmp = char(fread(this.VisaObj)');
            this.Channels.(chStr).v0    = str2double(tmp); 
            fprintf(this.VisaObj, ':WAVeform:YREFerence?' );   
%             pause(0.2);
            this.Channels.(chStr).v_ref = str2double(char(fread(this.VisaObj)'));  
            fprintf(this.VisaObj, ':WAVeform:YINCrement?' );   
%             pause(0.2);
            this.Channels.(chStr).dv    = str2double(char(fread(this.VisaObj)'));
            
%             fprintf(this.VisaObj, [':CHANnel',num2str(this.Channels.(chStr).id),':DISPlay ', 'OFF']);
            
%             fprintf(this.VisaObj ,[':WAV:SOUR CHAN', num2str(OrigCH)]);
        end
        
        function ConfAcq(this)
            if strcmp(this.Acq.mode, 'AVER')
                fprintf(this.VisaObj, ':ACQuire:TYPE AVER' );                     % Averaging ON
                fprintf(this.VisaObj, [':ACQuire:AVER ', num2str(this.Acq.numOfAvg)]); % Averaging max # of captures (limited by the Delay parameter)
            else
                fprintf(this.VisaObj, ':ACQuire:TYPE NORM' ); 
            end
        end
        
        function ConfTrigger(this)
           %Method for Configuring the settings of the trigger
              fprintf(this.VisaObj, [':TRIGger:MODE ', this.Trigger.mode]);
                if (strcmp(this.Trigger.mode, 'EDGE'))
                    %External Trigger control - irrelevant for DS1054z scope
                        %if (strcmpi(this.Trigger.trigCh, 'EXTernal')||strcmpi(this.Trigger.trigCh, 'EXT'))
                        %    %for using external trigger as trigger source. Both EXTernal and EXT
                        %    %are viable options to communicate to the scope by visa comm 
                        %    fprintf(this.VisaObj, [':TRIGger:EDGe:SOURce ', this.Trigger.trigCh]);
                        %end
                    fprintf(this.VisaObj, [':TRIGger:EDGe:SOURce CHANnel', num2str(this.Trigger.trigCh)]);
                    fprintf(this.VisaObj, [':TRIGger:EDGe:SLOPe ', this.Trigger.edge]);
                    fprintf(this.VisaObj, [':TRIGger:EDGe:LEVel ', num2str(this.Trigger.trigLevel/2)]);
                end
        end
            
        function ConfTimebase(this)
            %Method for Configuring the settings of the display Timebase (dt)
            fprintf(this.VisaObj, [':TIMebase:MODE ',this.Timebase.mode]);
            fprintf(this.VisaObj, [':TIMebase:MAIN:SCALe ', num2str(this.Timebase.scale, '%.9f')]);
            fprintf(this.VisaObj, [':TIMebase:MAIN:OFFSet ', num2str(this.Timebase.offset, '%.9f')]);
             
            fprintf(this.VisaObj, ':WAVeform:XORigin?' );      
            this.Timebase.t0    = str2double(char(fread(this.VisaObj)'));
            fprintf(this.VisaObj, ':WAVeform:XREFerence?' );      
            this.Timebase.tRef  = str2double(char(fread(this.VisaObj)'));
            fprintf(this.VisaObj, ':WAVeform:XINCrement?' );      
            this.Timebase.dt   = str2double(char(fread(this.VisaObj)'));
            
            this.Timebase.dispTimeVec = (0:1:1199)*this.Timebase.dt + this.Timebase.t0;
        end
            
        function ConfWaveform(this)
%            %Method for Configuring the settings of the Waveform capture
             fprintf(this.VisaObj, [':WAV:SOUR CHANnel', num2str(this.Waveform.chToSample)]); 
             fprintf(this.VisaObj, ':WAV:SOUR?' );
             fs = fread(this.VisaObj); disp(char(fs'));
             fprintf(this.VisaObj, [':WAV:MODE ', this.Waveform.WavPointsMode] ); % NORM / MAX / RAW. NORM means read the waveform from the screen. For detais see page 234 in the manual
             fprintf(this.VisaObj, ':WAV:MODE?' ); 
             fs = fread(this.VisaObj); disp(char(fs')); % Check scan mode
             fprintf(this.VisaObj, [':WAVeform:FORMat ', this.Waveform.Format] ); % ASCii / WORD / BYTE
             fprintf(this.VisaObj, ':WAV:FORMat?' );
             fs = fread(this.VisaObj, 2048); disp(char(fs')); % Check scan mode
        end

        % methods for opening/closing the connection with the Visa object
        function OpenConnection (this)
            try
                fopen(this.VisaObj); % Open the connection
                this.connected = true;
            catch
                fprintf("SCOPE: Couldnot connect to scope\n");
                this.connected = false;
            end
        end 
        
        function this = CloseConnection(this)
            % Close the VISA connection.
            if this.connected
                fclose(this.VisaObj);  
            end
            delete(this.VisaObj);
            clear this.VisaObj;
        end 
    end
end

