classdef AORigolGen < handle
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
        function this = AORigolGen(input)
            address = 'USB0::0x1AB1::0x0642::DG1ZA231702060::0::INSTR';
            if nargin ==1; address = input; end

            this.fg = Rigol_DG1022A('NI', address);
        end
        
        function init(this)
            this.vars.fGenVars.ch{1}.syncOwner     = true;
            this.vars.fGenVars.ch{1}.coup          = 50;
            this.vars.fGenVars.ch{2}.syncOwner     = false;
            this.vars.fGenVars.ch{2}.coup          = 50;
        end
        
        function setVars(this, uVars)
            this.vars.fgClk       = uVars.fgClk;
            
            this.vars.fSig        = uVars.fSig;
            this.vars.fSqnc       = uVars.fSqnc; % Repetition rate [Hz]
            this.vars.delay       = uVars.delay;
            this.vars.sigType     = uVars.sigType;
            this.vars.cycPerPulse = uVars.cycPerPulse;
            this.vars.fsClk       = uVars.fsClk; 
            this.vars.amp         = uVars.amp;
            
            this.vars.fGenVars.ch{2}.srate = uVars.fgClk;
            this.vars.fGenVars.ch{1}.srate = uVars.fgClk;
            
            if isempty(uVars.clkPower)
                this.vars.fGenVars.ch{2}.amp = 1;
            else
                this.vars.fGenVars.ch{2}.amp = 1 * uVars.clkPower/100;
            end

            if isempty(uVars.sigPower)
                this.vars.fGenVars.ch{1}.amp = 1;
            else
                this.vars.fGenVars.ch{1}.amp = 1 * uVars.sigPower/100 ;
            end

            this.vars.shift = uVars.shift;
            
            this.createSignals();
            
            this.vars.fGenVars.ch{1}.coup =  50;
            
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

%             factor = lcm( lcm(NgPulse, NgClkCyc), 16); %lcm = least common multiple
            factor = 1;
            
            NgSqnc = ceil(NgSqnc/factor)*factor;
            fSqnc  = fgClk / NgSqnc;
            
            pulsePerSqnc = floor(NgSqnc / NgPulse);

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

            sigData = this.vars.amp*circshift(sigData, this.vars.shift);
            
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
            this.vars.fGenVars.ch{1}.data = this.vars.dataCh1.data;
            this.fg.GenArb(1, this.vars.fGenVars.ch{1});
        end

        function closeOutputs(this)
            this.fg.disableAllOutputs();
        end
        
        function setDebug(this, debug)
            this.fg.debug = debug;
        end
    end
end



git 