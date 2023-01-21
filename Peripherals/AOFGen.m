classdef AOFGen < handle
    %AOFGEN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fg
        vars
    end

    methods (Static)
        function uVars = createUserVars()
            uVars.fgClk   =[];

            uVars.fSig    =[];
            uVars.fSqnc   =[];% Repetition rate [Hz]
            uVars.delay   =[];
            uVars.sigType =[];

            uVars.fsClk   =[];
            uVars.cycPerPulse = [];
            
            uVars.shift = 0;

            uVars.clkPower = [];
            uVars.sigPower = [];
        end
    end
    
    methods
        function this = AOFGen()
            this.fg = fGen();
            this.fg.connect();
        end
        
        function init(this)
            fgClk = 100e6;

            this.vars.fGenVars = this.fg.uVarsCreate();
            
            this.vars.fGenVars.ch{1}.daqFilter     = 50;
 
            this.vars.fGenVars.ch{1}.bias          = 0;
            this.vars.fGenVars.ch{1}.triggerOwner  = true;
            this.vars.fGenVars.ch{1}.triggerWidth  = 512;
            this.vars.fGenVars.ch{1}.useExtSclkSrc = false;

            this.vars.fGenVars.ch{2}.daqFilter     = 50;

            this.vars.fGenVars.ch{2}.bias          = 0;
            this.vars.fGenVars.ch{2}.triggerOwner  = false;
            this.vars.fGenVars.ch{2}.triggerWidth  = 16;
            this.vars.fGenVars.ch{2}.useExtSclkSrc = false;

            this.vars.fGenVars.ch{1}.Sclk = fgClk;
            this.vars.fGenVars.ch{2}.Sclk = fgClk;

            this.vars.dataCh1 = this.fg.getChDataStruct();
            this.vars.dataCh2 = this.fg.getChDataStruct();
        end
        
        function setVars(this, uVars)
            this.vars.fgClk = uVars.fgClk;
            
            this.vars.fSig    = uVars.fSig;
            this.vars.fSqnc   = uVars.fSqnc; % Repetition rate [Hz]
            this.vars.delay   = uVars.delay;
            this.vars.sigType = uVars.sigType;
            this.vars.cycPerPulse = uVars.cycPerPulse;
            this.vars.fsClk = uVars.fsClk; 
            
            if isempty(uVars.clkPower)
                this.vars.fGenVars.ch{2}.amp = 2;
            else
                this.vars.fGenVars.ch{2}.amp = 2 * uVars.clkPower/100;
            end

            if isempty(uVars.sigPower)
                this.vars.fGenVars.ch{1}.amp = 2;
            else
                this.vars.fGenVars.ch{1}.amp = 2 * uVars.sigPower/100 ;
            end

            this.vars.shift = uVars.shift;
            
            this.createSignals();
        end
        
        function curVars = getVars(this)
            curVars = this.vars;
        end
        
        function createClkSig(this)
            fs = this.vars.fsClk;
            fgClk = this.vars.fgClk;

            % The data loaded to the AFG must be a multiplication of 16 in manners of data length.
            % creating sClk with 16 cycles meets this requirement no matter
            % the fFgClk used.
            sClkCycles  = 16;
            sClkDutyCyc = 50/100;
            this.vars.fsNaive = fs;

            NgClkCyc = fgClk/fs;
            
            % If there is no round number of fgClk samples in sClk cycle, 
            % slow the sClk to the period where fgClk fits in.
            if mod(NgClkCyc,1) ~= 0
                NgClkCyc =  ceil(NgClkCyc);
                fs = fgClk / NgClkCyc;
                fprintf("Notice: The sampling frequence you have chosen cannot be genrated by the AFG. \n the closest sampling frequency is: %d\n", fs);
            end
            
            sClkT = 1/fs;
            NgClk = NgClkCyc * sClkCycles;
            
            % Create sClk Data
            cycleData                   = ones(NgClkCyc,1);
            dutyCycleIdx                = floor(NgClkCyc*(1-sClkDutyCyc))+1;
            cycleData(dutyCycleIdx:end) = 0;
            clkData                     = repmat(cycleData, sClkCycles, 1);
            
            %--------------------
            % Collect parameters
            %--------------------
            this.vars.sClk.fs    = fs;
            this.vars.sClk.sClkT = sClkT;
            
            this.vars.sClk.NgClkCyc = NgClkCyc;
            this.vars.sClk.NgClk    = NgClk;
            this.vars.clkData = clkData;
        end

        function createSig(this)
            % This function generates the US signal (single sequence) in manners of fgClk.
            % The US signal number of samples should be a multiply of:
            % Pulses, sClk, and 16 (fg demand).
            % This function is called after the function: calcSamplingClkDim
            % that makes sure that fgClk/sClk is a natural number.

            cycPerPulse = this.vars.cycPerPulse;
            fSig        = this.vars.fSig;
            fSqnc       = this.vars.fSqnc;
            fgClk       = this.vars.fgClk;

            % US Signal Vars
            fSqncNaive = fSqnc;
            NgSqnc     = floor(fgClk / fSqnc);   %this may vary later
            NgClkCyc   = this.vars.sClk.NgClkCyc;

            fPulse = fSig/cycPerPulse; 
            
            NgSin   = floor(fgClk / fSig);    % Hidden assumtion: fgClk/fSin is a Natural number.        
            NgPulse = floor(fgClk / fPulse); 

            factor = lcm( lcm(NgPulse, NgClkCyc), 16); %lcm = least common multiple
            
            NgSqnc = ceil(NgSqnc/factor)*factor;
            fSqnc  = fgClk / NgSqnc;
            
            pulsePerSqnc = NgSqnc / NgPulse;

            if strcmp(this.vars.sigType, 'Hadamard')
                sMatrix = createSMatrix(pulsePerSqnc, 'QR'); 
                NHad       = size(sMatrix, 1);
                
                %createSMatrix returns an S-matrix with minimal order of N 
                %that satisfies: N>=pulsePerSqnc 
                %In case N~=pulsePerSqnc, fSqnc needs to be changed to fit N pulses.
                if NHad ~= pulsePerSqnc
                    NgSqnc = NgSqnc + (NHad-pulsePerSqnc)*NgPulse;
                    fSqnc               = fgClk / NgSqnc;
                    pulsePerSqnc        = NHad;
                    fprintf("NOTICE: Multiplexing matrix size is different than num of pos\n")
                end
            end
            
            samplesVec = 0: 1 : NgPulse-1;
            
            sigVec = zeros(1, pulsePerSqnc);
            sigVec(1) = 1;

            switch this.vars.sigType
                case 'Hadamard'
                    sigVec = sMatrix(1, :);
                    sigSamples = sin(2*pi * (fSig/fgClk) * samplesVec);
                case 'SinePulse'
                    sigSamples = sin(2*pi * (fSig/fgClk) * samplesVec);
                case 'CW'
                    sigVec = ones(1, pulsePerSqnc);
                    sigSamples = sin(2*pi * (fSig/fgClk) * samplesVec);
                case 'DeltaPulse'
                    sigSamples    = zeros(1,NgPulse);
                    sigSamples(1) = 1;
                case 'RectPulse'
                    sigSamples = ones(1,NgPulse);
                case 'BiPolarRectPulse'
                    if mod(NgPulse, 2)
                        lenPos = floor(NgPulse/2);
                        lenNeg = floor(NgPulse/2)+1;
                    else
                        lenPos = NgPulse/2;
                        lenNeg = NgPulse/2;
                    end
                    sigSamples = [ones(1,lenPos), -1*ones(1,lenNeg)];
            end
            
            sigData = repmat(sigVec, NgPulse, 1);
            idxs = find(sigVec);
            sigData(:, idxs) = repmat(sigSamples', 1, length(idxs));
            sigData = sigData(:);

            sigData = circshift(sigData, this.vars.shift);
            %--------------------
            % Collect parameters
            %--------------------
            this.vars.sig.NgSin   = NgSin;
            this.vars.sig.NgPulse = NgPulse;
            this.vars.sig.NgSqnc  = NgSqnc;

            this.vars.sig.fSqncNaive = fSqncNaive;
            this.vars.sig.fSqnc      = fSqnc;
            
            this.vars.sig.pulsePerSqnc = pulsePerSqnc;

            this.vars.sigData = sigData;
        end

        function createSignals(this)
            this.createClkSig();
            this.vars.dataCh2.data    = this.vars.clkData;
            this.vars.dataCh2.dataLen = length(this.vars.dataCh2.data);

            this.createSig();
            this.vars.dataCh1.data    = this.vars.sigData;
            this.vars.dataCh1.dataLen = length(this.vars.dataCh1.data);
        end

        function config(this)
            this.fg.reset();
            this.fg.setProperties(this.vars.fGenVars.ch{1}, this.vars.fGenVars.ch{2});
            this.fg.setData(this.vars.dataCh1, this.vars.dataCh2);
            this.fg.configChannel(1); % this will enable the channel
            this.fg.configChannel(2); % this will enable the channel
        end

        function closeOutputs(this)
            this.fg.disableAllOutputs();
        end
    end
end

