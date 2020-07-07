classdef Algo < handle
    %ALGO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        uVars;
        extClk;
        usSignal;
        geometry;
        digitizer;
        samples;
        timing;
        len;
        freq;
        hadamard;
        general;
        
        data
        reshapedData
        Wn
        
        measLimit
        timeTable;
        res;
    end
 
    methods (Static)   
        
        function uVars = uVarsCreate()
            uVars.fSin       = [];              
            uVars.fTrain     = [];
            uVars.cycInPulse = [];
            uVars.channels   = [];
            uVars.c                 = [];  
            uVars.phantomDepth      = [];
            uVars.distFromPhantom   = [];
            uVars.bufferSizeBytes   = [];
            uVars.bytesPerSample    = [];
            uVars.preTriggerSamples = [];
            uVars.fExtClk           = []; %fs
            uVars.fSclk             = [];
            uVars.timeToSample      = [];
            uVars.extClkDcyc        = []; % [%]
            uVars.envDC             = [];
            uVars.envUS             = [];
            uVars.useQuant          = false;
            uVars.quantTime         = [];
            uVars.fastAnalysis      = false;
            uVars.useGPU            = false;
            uVars.useHadamard       = false;
            
            uVars.export.netSignal      = true;
            uVars.export.deMultiplexed  = true;
            uVars.export.reshapedSignal = true;
            uVars.export.fftRes         = true;
            uVars.export.usCompCmplx    = true;
        end
       
        function gReq = createGraphicRequest()
            gReq = algoGraphics.createGraphicsRunVars();
        end
        
    end
    
    methods
        function this = Algo()
            this.initTimeTable();  
            this.measLimit = 20;
        end
        
        function setMeasLimit(this, limit)
           this.measLimit = limit; 
        end
        
        function initTimeTable(this)
            this.timeTable.copyFullRawData          = 0;
            
            this.timeTable.CopyMeasSignal           = 0;
            this.timeTable.CutDigiPreSamples        = 0;
            this.timeTable.CutPropPreSamples        = 0;
            this.timeTable.CutInPhantopmPropSamples = 0;
            this.timeTable.CutPostSamples           = 0;
            this.timeTable.CopyNetSignal            = 0;
            this.timeTable.ExtractNetSignal         = 0;
            
            this.timeTable.deMulReshape1            = 0;
            this.timeTable.deMultiplex              = 0; 
            this.timeTable.deMulReshape2            = 0;
            this.timeTable.CopyDeMultiplexed        = 0;
            
            this.timeTable.ReshapeQuant             = 0;
            this.timeTable.permuteQuant             = 0;
            this.timeTable.Reshape1                 = 0;
            this.timeTable.Reshape2                 = 0;
            this.timeTable.permuteReshape           = 0;
            this.timeTable.Reshape3                 = 0;
            this.timeTable.copyReshaped             = 0;
            this.timeTable.OverallReshape           = 0;
            
            this.timeTable.totalFastAnalysis        = 0;
            this.timeTable.FFT                      = 0;
            this.timeTable.copyFFT                  = 0;
            this.timeTable.savePhases               = 0;
            this.timeTable.meanQuantFFT             = 0;
            this.timeTable.meanChFFT                = 0;
            this.timeTable.fftShift                 = 0;
            this.timeTable.createFit                = 0;
            this.timeTable.subFit                   = 0;
            this.timeTable.phiExtract               = 0;
            this.timeTable.SignalProcessing         = 0;
            
            this.timeTable.FullAnalysis             = 0;
        end

        function [vars] = updateAlgoUserVars(this, user)
           this.uVars.fSin              = user.fSin; %V
           this.uVars.fTrain            = user.fTrain; %V
           this.uVars.cycInPulse        = user.cycInPulse; %V
           this.uVars.channels          = user.channels; %V
           this.uVars.timeToSample      = user.timeToSample;
           this.uVars.c                 = user.c;
           this.uVars.phantomDepth      = user.phantomDepth;
           this.uVars.distFromPhantom   = user.distFromPhantom;
           this.uVars.bufferSizeBytes   = user.bufferSizeBytes;
           this.uVars.bytesPerSample    = user.bytesPerSample;
           this.uVars.preTriggerSamples = user.preTriggerSamples;
           this.uVars.fSclk             = user.fSclk;
           this.uVars.fExtClk           = user.fExtClk; %raw fs
           this.uVars.extClkDcyc        = user.extClkDcyc;
           this.uVars.envDC             = user.envDC;
           this.uVars.envUS             = user.envUS;
           this.uVars.envHar            = user.envHar;
           this.uVars.useQuant          = user.useQuant;
           this.uVars.quantTime         = user.quantTime;
           this.uVars.fastAnalysis      = user.fastAnalysis;
           this.uVars.useGPU            = user.useGPU;
           this.uVars.useHadamard       = user.useHadamard;
           this.uVars.longMeas          = user.longMeas;
           this.uVars.export            = user.export;

           vars = this.calcDimensions();
        end
        
        function calcGeneral(this)
           this.general.longMeas = this.uVars.longMeas;  
           this.general.fitToShift = true;
           this.general.useGPU = this.uVars.useGPU;
           this.general.useHadamard = this.uVars.useHadamard;
           
           this.general.export.netSignal      = this.uVars.export.netSignal      && ~this.general.longMeas;
           this.general.export.deMultiplexed  = this.uVars.export.deMultiplexed  && ~this.general.longMeas ;
           this.general.export.reshapedSignal = this.uVars.export.reshapedSignal && ~this.general.longMeas;
           this.general.export.fftRes         = this.uVars.export.fftRes         && ~this.general.longMeas;
           this.general.export.usCompCmplx    = this.uVars.export.usCompCmplx    && ~this.general.longMeas;
        end
        
        function vars = calcDimensions(this)
            this.calcGeneral();
            this.calcGeometry();
            this.calcExtClkDim();
            this.calcUSSignalDim();
            this.calcSamples();
            this.calcDigitizerParams();
            this.calcTiming();
            this.calcLength();
            this.calcFreq();
            
            this.allocateReshapeMemory()
            
            vars = this.getVars();
        end
        
        function allocateReshapeMemory(this)
            this.reshapedData = [];
        end
        
        function calcGeometry(this)
            this.geometry.c = this.uVars.c;
            this.geometry.distFromPhantom = this.uVars.distFromPhantom;
            this.geometry.phantomDepth = this.uVars.phantomDepth;
        end
        
        function vars = getVars(this)
            vars.extClk    = this.extClk;
            vars.usSignal  = this.usSignal;
            vars.samples   = this.samples;
            vars.timing    = this.timing;
            vars.len       = this.len;
            vars.freq      = this.freq;
            vars.geometry  = this.geometry;
            vars.digitizer = this.digitizer;
            vars.uVars     = this.uVars;
        end
          
        function calcExtClkDim(this)
            % External Clock Params
            this.extClk.fExtClk             = this.uVars.fExtClk;
            this.extClk.fSclk               = this.uVars.fSclk;
            this.extClk.extClkSamplesPerCyc = this.extClk.fSclk/this.extClk.fExtClk;
            this.extClk.extClkCycles        = 16;
            this.extClk.fExtClkNaive        = this.extClk.fExtClk;
            this.extClk.extClkDcyc          = this.uVars.extClkDcyc;

            % if there is no round number of sclk in on extClk cycle, slow the extClk
            % to the period where sclk fits in.
            if mod(this.extClk.extClkSamplesPerCyc,1)~=0
                this.extClk.extClkSamplesPerCyc  =  ceil(this.extClk.extClkSamplesPerCyc);
                this.extClk.fExtClk              = this.extClk.fSclk / this.extClk.extClkSamplesPerCyc;
                fprintf("Notice: The clock rate you have chosen cannot be genrated by the AFG. the closest rate is: %d\n", this.extClk.extClkTrueFreq);
            end

            this.extClk.tExtClk          = 1/this.extClk.fExtClk; 
            this.extClk.extClkSigSamples = this.extClk.extClkSamplesPerCyc*this.extClk.extClkCycles;
            this.digitizer.fs = this.extClk.fExtClk;
        end
        
        function calcUSSignalDim(this)
            % US Excitation params
            % Make sure that the train is a multiple of (in manner of sclk): Pulses, extClk cycles and 16.
            this.usSignal.cycInPulse = this.uVars.cycInPulse;
            this.usSignal.fSin       = this.uVars.fSin;
            this.usSignal.fTrain     = this.uVars.fTrain;
 
            this.usSignal.SclkSamplesInPulse = this.usSignal.cycInPulse*(this.extClk.fSclk/this.usSignal.fSin);
            this.usSignal.fTrainNaive        = this.usSignal.fTrain;
            this.usSignal.SclkSamplesInTrain = this.extClk.fSclk/this.usSignal.fTrain;
            factor                           = (this.usSignal.cycInPulse * (this.extClk.fExtClk/this.usSignal.fSin)) * ...
                                                this.extClk.extClkSamplesPerCyc;
            this.usSignal.SclkSamplesInTrain = ceil(this.usSignal.SclkSamplesInTrain/factor)*factor;
            this.usSignal.fTrain             = this.extClk.fSclk / this.usSignal.SclkSamplesInTrain;
            
            if this.uVars.useHadamard
                % in future, n should be calculated out of 
                % phantom height/pulse len
                n =  this.usSignal.SclkSamplesInTrain / this.usSignal.SclkSamplesInPulse;
                this.hadamard.sMatrix = createSMatrix(n, 'QR');
%                 disp(det(this.hadamard.sMatrix))
                this.hadamard.sMatInv = inv(this.hadamard.sMatrix);
                this.hadamard.n = size( this.hadamard.sMatrix, 1);
                if this.hadamard.n ~= n
                    this.usSignal.SclkSamplesInTrain = this.usSignal.SclkSamplesInTrain + (size(this.hadamard.sMatrix, 1)-n)*this.usSignal.SclkSamplesInPulse;
                    this.usSignal.fTrain = this.extClk.fSclk / this.usSignal.SclkSamplesInTrain;
                    fprintf("NOTICE: Multiplexing matrix size is different than num of pos\n")
                end
            end

        end
        
        function calcSamples(this)
            if ~this.uVars.useQuant
                this.uVars.quantTime = this.uVars.timeToSample;
            end

            this.samples.samplesPerSin         = this.digitizer.fs / this.usSignal.fSin;
            this.samples.samplesPerPulse       = this.samples.samplesPerSin * this.usSignal.cycInPulse;
%                 this.samples.samplesPerTrain       = floor(this.digitizer.fs / this.usSignal.fTrain);
            this.samples.samplesPerTrain       = this.usSignal.SclkSamplesInTrain / (this.extClk.fSclk / this.digitizer.fs);

            this.samples.numOfTrains      = ceil(this.uVars.timeToSample * this.usSignal.fTrain);
            this.samples.trainsPerQuant   = ceil(this.uVars.quantTime *this.usSignal.fTrain);
            this.samples.numOfTrains      = ceil(this.samples.numOfTrains/this.samples.trainsPerQuant)*this.samples.trainsPerQuant;
            this.samples.numOfTrainsAllCh = this.samples.numOfTrains * this.uVars.channels;
            this.samples.numOfQuant       = this.samples.numOfTrains / this.samples.trainsPerQuant;
            
            this.samples.samplesPerQuant       = this.samples.trainsPerQuant*this.samples.samplesPerTrain;
            this.samples.samplesPerPos         = this.samples.samplesPerPulse * this.samples.trainsPerQuant; %this is practically perPosPerQuant
            this.samples.prePhantomSamples     = ceil( (this.geometry.distFromPhantom / this.geometry.c) * this.digitizer.fs );
            this.samples.inPhantomPropSamples  = this.samples.samplesPerTrain - this.samples.samplesPerPulse; 
            this.samples.samplesPerSignal      = this.samples.samplesPerTrain * this.samples.numOfTrains;
            
            if this.uVars.useHadamard
                this.samples.samplesPerMeas = this.uVars.preTriggerSamples + this.samples.prePhantomSamples + this.samples.inPhantomPropSamples + this.samples.samplesPerSignal;
                this.samples.preCut         = this.uVars.preTriggerSamples + this.samples.prePhantomSamples + this.samples.inPhantomPropSamples;
            else
                this.samples.samplesPerMeas = this.uVars.preTriggerSamples + this.samples.prePhantomSamples + this.samples.samplesPerSignal;
                this.samples.preCut         = this.uVars.preTriggerSamples + this.samples.prePhantomSamples;
            end
            this.samples.signalInd = this.samples.preCut + this.samples.samplesPerSignal;
            
            this.samples.samplesPerSignalAllCh = this.samples.samplesPerSignal * this.uVars.channels;
            this.samples.samplesPerMeasAllCh   = this.samples.samplesPerMeas   * this.uVars.channels;

            this.samples.numOfPos        = this.samples.samplesPerTrain / this.samples.samplesPerPulse;
            this.samples.samplesPerZAxis = this.samples.numOfPos * this.samples.samplesPerPulse;
            
            this.timing.trueQuantTime = this.samples.samplesPerQuant /this.uVars.fExtClk;
            this.timing.trainTime = this.samples.samplesPerTrain / this.uVars.fExtClk;
            
            if this.general.longMeas 
               this.samples.trainsPerLimit = floor(this.measLimit/this.timing.trainTime);
               this.samples.trainsPerChunk = this.uVars.channels * this.samples.trainsPerLimit;
               this.samples.numOfChunks    = ceil(this.samples.numOfTrainsAllCh/this.samples.trainsPerChunk);
               this.samples.trainsInLastChunk = this.samples.numOfTrainsAllCh - (this.samples.numOfChunks-1)*this.samples.trainsPerChunk;
            end
            
            if this.uVars.useHadamard  
                % create sparse matrix for single vec demultiplexing (single sVec)
                this.hadamard.nSVecPad = this.samples.numOfPos * this.samples.samplesPerPulse; %samplesPerTrain
                this.hadamard.singleSVecDeMulMat = zeros(this.samples.numOfPos, this.hadamard.nSVecPad);
                this.hadamard.singleSVecDeMulMat(:,1:(this.samples.samplesPerPulse):end) = this.hadamard.sMatInv;

                % Expand the matrix to keep order of demultiplexing
                this.hadamard.singleTrainMat = kron( flip(this.hadamard.singleSVecDeMulMat, 1), ones(  this.samples.samplesPerPulse,1) );

                % create sparse matrix for a single train
                this.hadamard.singleTrainDeMulMat = zeros(this.hadamard.nSVecPad, this.hadamard.nSVecPad);
                for i = 1:this.samples.samplesPerPulse
                    this.hadamard.singleTrainDeMulMat( i:this.samples.samplesPerPulse:end, :) = ...
                                                        circshift(this.hadamard.singleTrainMat( i:this.samples.samplesPerPulse:end, :), (i-1), 2);
                end
                
                this.hadamard.singleTrainDeMulMat = single(this.hadamard.singleTrainDeMulMat);
                
                if this.uVars.useGPU
                    this.hadamard.singleTrainDeMulMat = gpuArray(this.hadamard.singleTrainDeMulMat);
                end
            end
        end
        
        function calcTiming(this)
            this.timing.tSin   = 1/this.usSignal.fSin;                         %[s], duration of one period of the sin
            this.timing.tPulse = this.usSignal.cycInPulse * this.timing.tSin; %[s], duration of the pulse
            this.timing.tTrain = 1/this.usSignal.fTrain;                       %[s], duration of one period of the train
            
            this.timing.dts     = 1/this.digitizer.fs;                          %[s], duration of a single sample
            this.timing.dtSclk  = 1/this.extClk.fSclk;
            
            this.timing.tQuant            = this.samples.samplesPerQuant*this.timing.dts;
            this.timing.timeOfSample      = this.timing.dts*this.samples.samplesPerSignal; %  should it be samplesPerSignal -1?
            this.timing.actualSampledTime = this.digitizer.samplesPerAcq * this.timing.dts;
            
            this.timing.tExtClk    = (0:1:(this.extClk.extClkSigSamples-1)    ) * this.timing.dtSclk;
            this.timing.tUS        = (0:1:(this.usSignal.SclkSamplesInTrain-1)) * this.timing.dtSclk;
            if ~this.general.longMeas
                this.timing.tPulseVec  = (0:1:this.samples.samplesPerPulse  - 1)*this.timing.dts;
                this.timing.tTrainMeas = (0:1:this.samples.samplesPerTrain  - 1)*this.timing.dts;
                this.timing.tSigVec    = (0:1:this.samples.samplesPerSignal - 1)*this.timing.dts;
                this.timing.tMeasVec   = (0:1:this.samples.samplesPerMeas   - 1)*this.timing.dts;
                this.timing.tPosVec    = (0:1:this.samples.samplesPerPos    - 1)*this.timing.dts;
                this.timing.tAcqVec    = (0:1:this.digitizer.samplesPerAcq  - 1)*this.timing.dts;
                this.timing.tQuantVec  = (0:1:this.samples.samplesPerQuant  - 1)*this.timing.dts;
            else
                this.timing.tPulseVec  = [];
                this.timing.tTrainMeas = [];
                this.timing.tSigVec    = [];
                this.timing.tMeasVec   = [];
                this.timing.tPosVec    = [];
                this.timing.tAcqVec    = [];
                this.timing.tQuantVec  = [];
            end
        end
        
        function calcDigitizerParams(this)
            this.digitizer.channels          = this.uVars.channels;
            this.digitizer.bufferSizeBytes   = this.uVars.bufferSizeBytes;
            this.digitizer.bytesPerSample    = this.uVars.bytesPerSample;
            this.digitizer.preTriggerSamples = this.uVars.preTriggerSamples;
            
            this.digitizer.samplesPerBufferAllCh = this.digitizer.bufferSizeBytes / this.digitizer.bytesPerSample;
            this.digitizer.samplesPerBuffer      = this.digitizer.samplesPerBufferAllCh / this.digitizer.channels;
            this.digitizer.numOfBuffers          = ceil(this.samples.samplesPerMeas / this.digitizer.samplesPerBuffer);
            this.digitizer.samplesPerAcq         = this.digitizer.samplesPerBuffer * this.digitizer.numOfBuffers;
            this.digitizer.samplesPerAcqAllCh    = this.digitizer.samplesPerAcq * this.digitizer.channels;
        end
        
        function calcLength(this)
            this.len.dzs      = this.geometry.c*this.timing.dts;
            this.len.sinLen   = this.geometry.c * this.timing.tSin;
            this.len.pulseLen = this.geometry.c*this.timing.tPulse;
            this.len.zRes     = this.len.pulseLen; 
            this.len.trainLen = this.geometry.c*this.timing.tTrain;
            
            this.len.zLen =  this.samples.samplesPerZAxis * this.len.dzs; 
            this.len.zVec = (0:1:(this.samples.samplesPerZAxis-1))*this.len.dzs;
            this.len.zVecUSRes = (0:1:(this.samples.numOfPos-1)) * this.len.zRes;
            this.len.zIdx =  1:1:this.samples.numOfPos;
            this.len.zIdxLen = this.samples.numOfPos;
        end
        
        function calcFreq(this)
            k                             = 0:this.samples.samplesPerPos-1;
            this.freq.frequencyBar        = this.digitizer.fs * k / this.samples.samplesPerPos;
            this.freq.frequencyBarShifted = [this.freq.frequencyBar( (this.samples.samplesPerPos/2+1):end ) - this.digitizer.fs,...
                                             this.freq.frequencyBar( 1:(this.samples.samplesPerPos/2) )] ;
            this.freq.df                  = abs(this.freq.frequencyBar(2) - this.freq.frequencyBar(1));
            this.freq.fs                  = this.extClk.fExtClk;
%             this.Wn = exp((-1j*(2*pi)/this.samples.samplesPerPos)*(this.freq.fUSIdx-1)*(0:(this.samples.samplesPerPos-1)))';
            this.freq.envDC  =  this.uVars.envDC;
            this.freq.envUS  =  this.uVars.envUS;
            this.freq.envHar =  this.uVars.envHar;
            this.freq.samplesFFT = this.samples.samplesPerPos;
            this.freq.fUS = this.usSignal.fSin;
            
            this.freq.fIdxPos = ((this.freq.fUS * this.freq.samplesFFT)/this.freq.fs)+1;
            this.freq.fIdxNeg = ( ( (this.freq.fUS + this.freq.fs/2) * this.freq.samplesFFT)/this.freq.fs)+1;
            this.freq.DcIdx   = 1;
            
            this.freq.shiftFactor = floor(this.freq.samplesFFT/2);
            
            this.freq.fIdxPosShift = this.freq.fIdxPos + this.freq.shiftFactor;
            this.freq.fIdxNegShift = this.freq.fIdxNeg - this.freq.shiftFactor;
            this.freq.DcIdxShift   = this.freq.DcIdx   + this.freq.shiftFactor;
            this.freq.envUsIdx     = sum(abs(this.freq.frequencyBar-this.freq.fUS) < this.freq.envUS);
            this.freq.envDcIdx     = sum(abs(this.freq.frequencyBar) < this.freq.envDC);
            this.freq.envHarIdx    = sum(abs(this.freq.frequencyBar - this.freq.fs) < this.freq.envHar); 
            if this.general.fitToShift
                this.freq.fUSIdx = this.freq.fIdxPosShift;
            else
                this.freq.fUSIdx = this.freq.fIdxPos;
            end
            
            this.freq.fitIdxShift = [ this.freq.envHarIdx                          :(this.freq.fIdxNegShift - this.freq.envUsIdx),...
                                     (this.freq.fIdxNegShift + this.freq.envUsIdx) :(this.freq.DcIdxShift   - this.freq.envDcIdx),...
                                     (this.freq.DcIdxShift   + this.freq.envDcIdx) :(this.freq.fIdxPosShift - this.freq.envUsIdx), ...
                                     (this.freq.fIdxPosShift + this.freq.envUsIdx) : this.freq.samplesFFT   - this.freq.envHarIdx];
            
            this.freq.fitIdx = [ this.freq.envDcIdx                      : (this.freq.fIdxPos    - this.freq.envUsIdx),...
                                (this.freq.fIdxPos + this.freq.envUsIdx) : (this.freq.fIdxNeg    - this.freq.envUsIdx),...
                                (this.freq.fIdxNeg + this.freq.envUsIdx) : (this.freq.samplesFFT - this.freq.envDcIdx)];
            
                            
            if this.uVars.useGPU
                this.freq.frequencyBar        = gpuArray(this.freq.frequencyBar);
                this.freq.frequencyBarShifted = gpuArray(this.freq.frequencyBarShifted);
                this.freq.fitIdxShift         = gpuArray(this.freq.fitIdxShift);
                this.freq.fitIdx              = gpuArray(this.freq.fitIdx);
            end
            
            this.freq.linearft = fittype(  {'1',...
                                            'x',...
                                            '(1/2)*(3*x^2 - 1)',...
                                            '(1/2)*(5*x^3 - 3*x)',...
                                            '(1/8)*(35*x^4 - 30*x^2 + 3)',...
                                            '(1/8)*(63*x^5 - 70*x^3 + 15*x)',...
                                            '(1/16)*(231*x^6 - 315*x^4 + 105*x^2 - 5)',...
                                            '(1/16)*(429*x^7 - 693*x^5 + 315*x^3 - 35*x)',...
                                            '(1/128)*(6435*x^8 - 12012*x^6 + 6930*x^4 - 1260*x^2 + 35)',...
                                            '(1/128)*(12155*x^9 - 25740*x^7 + 18018*x^5 - 4620*x^3 + 315*x)',...
                                            '(1/256)*(46189*x^10 - 109395*x^8 + 90090*x^6 - 30030*x^4 + 3465*x^2 - 63)'});
                                        
            warning('off','curvefit:fit:equationBadlyConditioned');
        end
        
        function [sigData, clkData] = createSignalsForfGen(this)
            % Generate the Signal
            dt = 1/this.extClk.fSclk;
            tSig = (0:1:(this.usSignal.SclkSamplesInTrain-1)) * dt;
            
            if ~this.uVars.useHadamard
                sigData = zeros(this.usSignal.SclkSamplesInTrain, 1);
                sigData(1:this.usSignal.SclkSamplesInPulse) = sin(2*pi*this.usSignal.fSin*tSig(1:this.usSignal.SclkSamplesInPulse));
            else
                sigData = reshape(repmat(this.hadamard.sMatrix(1,:), this.usSignal.SclkSamplesInPulse, 1), 1, this.usSignal.SclkSamplesInTrain);
                sigData = sigData.*sin(2*pi*this.usSignal.fSin*tSig);
            end
            %Generate The Clock
            clkData      = ones(this.extClk.extClkSamplesPerCyc, 1);
            dutyCycleIdx = floor(this.extClk.extClkSamplesPerCyc*(this.extClk.extClkDcyc/100))+1;
            clkData(dutyCycleIdx:end) = 0;
            clkData = repmat(clkData, this.extClk.extClkCycles, 1);

            this.extClk.data = clkData;
            this.usSignal.data = sigData;  
        end      
        
        function res = analyse(this)
            % rawData(input)  - [ch x samplesPerAcq]
            % res(output) - struct:
            %               rawData - [ch x samplesPerSignal]
            %               signal  - [ch x samplesPerPos x numOfPos]
            %               fftRes  - [ch x samplesPerPos x numOfPos]
            %               phiCh   - [ch x 1 x numOfPos] 
            %               phi     - [1  x 1 x numOfPos]          
            % TODO: check if the instances of the signal contribute to
            % algorith, slowdown.
            
            this.timeTable.FullAnalysis = tic;
            this.res = [];
            
            this.timeTable.ExtractNetSignal = tic;
            this.extractNetSignal();
            this.timeTable.ExtractNetSignal = toc(this.timeTable.ExtractNetSignal);
            
            if this.uVars.useHadamard
                this.demultiplexSignal();
            end
            
            this.timeTable.OverallReshape = tic;
            this.reshapeSignal();
            this.timeTable.OverallReshape = toc(this.timeTable.OverallReshape);

            this.timeTable.SignalProcessing = tic;
            this.signalProcessing();
            this.timeTable.SignalProcessing = toc(this.timeTable.SignalProcessing);

            res = this.res;
            
             %all data in res is CPU
            this.data = [];
            if this.general.longMeas
                this.res = [];
            end
            this.timeTable.FullAnalysis = toc(this.timeTable.FullAnalysis);
        end
        
        function extractNetSignal(this)
            % rawData(input)  - [ch x samplesPerAcq]
            % rawData(output) - [ch x samplesPerSignal]
           
            % Chop unneccesarry cycles 
%             digiPreSamplesIdx = this.digitizer.preTriggerSamples;
%             propPreSamplesIdx = this.samples.prePhantomSamples;
%             propInPhantomIdx  = this.samples.inPhantomPropSamples;
%             
%             if this.uVars.useHadamard
%                 preCut    = digiPreSamplesIdx + propPreSamplesIdx + propInPhantomIdx;
%                 signalInd = preCut + this.samples.samplesPerSignal;
%             else
%                 preCut    = digiPreSamplesIdx + propPreSamplesIdx;
%                 signalInd = preCut + this.samples.samplesPerSignal;
%             end
            this.data = this.data(:,this.samples.preCut+1:this.samples.signalInd);
            
            if this.general.export.netSignal
                this.timeTable.CopyNetSignal = tic;
                this.res.netSignal = gather(this.data);
                this.timeTable.CopyNetSignal = toc(this.timeTable.CopyNetSignal);
            end

        end
        
        function demultiplexSignal(this)
            this.timeTable.deMulReshape1 = tic;
            this.data = reshape(this.data', this.samples.samplesPerTrain, this.samples.numOfTrainsAllCh);
            this.timeTable.deMulReshape1 = toc(this.timeTable.deMulReshape1);
            
            this.timeTable.deMultiplex = tic;
            if this.general.longMeas
                
                for i=1:this.samples.numOfChunks
                    if i == this.samples.numOfChunks
                        this.data(:, (i-1)*this.samples.trainsPerChunk+1: end) =...
                            this.hadamard.singleTrainDeMulMat * this.data(:, (i-1)*this.samples.trainsPerChunk+1: end);
                    else
                        this.data(:, (i-1)*this.samples.trainsPerChunk+1: i*this.samples.trainsPerChunk) =...
                            this.hadamard.singleTrainDeMulMat * this.data(:, (i-1)*this.samples.trainsPerChunk+1: i*this.samples.trainsPerChunk);
                    end
                end
            else
                this.data = (this.hadamard.singleTrainDeMulMat * this.data);
            end
            this.timeTable.deMultiplex = toc(this.timeTable.deMultiplex);
            
            this.timeTable.deMulReshape2 = tic;
            this.data = reshape(this.data, this.samples.samplesPerSignal, this.digitizer.channels)';
            this.timeTable.deMulReshape2 = toc(this.timeTable.deMulReshape2);
            
            this.timeTable.CopyDeMultiplexed = tic;
            if this.general.export.deMultiplexed
                this.res.deMultiplexed = gather(this.data);
            end
            this.timeTable.CopyDeMultiplexed = toc(this.timeTable.CopyDeMultiplexed);
                       
        end
        
        function reshapeSignal(this)
           % signal(input)  - [ch x samplesPerSignal]
           % signal(output) - [ch x samplesPerPos x numOfPos]
                % Separate into Quants
                this.timeTable.ReshapeQuant = tic;
                this.data = reshape(this.data, this.digitizer.channels,...
                                               this.samples.samplesPerQuant, ...
                                               this.samples.numOfQuant);
                this.timeTable.ReshapeQuant = toc(this.timeTable.ReshapeQuant);
                
                %Bring quant-dim to be first (averaging convinient)
                this.timeTable.permuteQuant = tic;
                this.data = permute(this.data, [3,1,2]);
                this.timeTable.permuteQuant = toc(this.timeTable.permuteQuant);
                
                this.timeTable.Reshape1 = tic;
                this.data = reshape(this.data, this.samples.numOfQuant,...
                                               this.digitizer.channels,...
                                               this.samples.samplesPerPulse,...
                                               this.samples.trainsPerQuant * this.samples.numOfPos);
                this.timeTable.Reshape1 = toc(this.timeTable.Reshape1);
                
                this.timeTable.Reshape2 = tic;
                this.data = reshape(this.data, this.samples.numOfQuant,...
                                               this.digitizer.channels,... 
                                               this.samples.samplesPerPulse,...
                                               this.samples.numOfPos,this.samples.trainsPerQuant);
                this.timeTable.Reshape2 = toc(this.timeTable.Reshape2); 
                
                this.timeTable.permuteReshape = tic;
                this.data = permute(this.data, [1,2,3,5,4]);
                this.timeTable.permuteReshape = toc(this.timeTable.permuteReshape); 
                
                this.timeTable.Reshape3 = tic;
                this.data = reshape(this.data, this.samples.numOfQuant,...
                                               this.digitizer.channels,...
                                               this.samples.samplesPerPos,...
                                               this.samples.numOfPos);
                this.timeTable.Reshape3 = toc(this.timeTable.Reshape3);                           
                 
                this.timeTable.copyReshaped = tic;
                if this.general.export.reshapedSignal
                    this.res.reshapedSignal = gather(this.data);
                end
                this.timeTable.copyReshaped = toc(this.timeTable.copyReshaped);
        end
        
        function signalProcessing(this)
            % data(input)              - [quants x ch x samplesPerPos x numOfPos]
            % fftRes(output)           - [quants x ch x samplesPerPos x numOfPos]
            % usCompCmplx(output)      - [quants x ch x numOfPos]
            % qAvgChFFT(output)        - [ch x samplesPerPos x numOfPos]
            % unFittedFFT(output)      - [samplesPerPos x numOfPos]
            % unFittedFFTShift(output) - [samplesPerPos x numOfPos]
            % fitModel(output)         - [samplesPerPos x numOfPos]
            % fittedFFT(output)        - [samplesPerPos x numOfPos]
            % phi(output)              - [1 x numOfPos]
            
            this.timeTable.FFT = tic;
            this.data = (2./this.samples.samplesPerPos) * fft(this.data,[],3);

            this.timeTable.FFT = toc(this.timeTable.FFT);

            this.timeTable.copyFFT = tic;
            if this.general.export.fftRes
                this.res.fftRes = gather(this.data);
            end
            this.timeTable.copyFFT = toc(this.timeTable.copyFFT);
            
            this.timeTable.savePhases = tic;
            if this.general.export.usCompCmplx
                this.res.usCompCmplx = squeeze(this.data(:,:, this.freq.fIdxPos, :)); %quite large
            end
            this.timeTable.savePhases = toc(this.timeTable.savePhases);
            
            this.timeTable.meanQuantFFT = tic;
            this.res.qAvgChFFT = squeeze(mean(abs(this.data).^2, 1));
            this.timeTable.meanQuantFFT = toc(this.timeTable.meanQuantFFT);
            
%             if ~this.general.longMeas
                this.timeTable.meanChFFT = tic;
                this.res.unFittedFFT = squeeze(mean(this.res.qAvgChFFT, 1));
                this.timeTable.meanChFFT = toc(this.timeTable.meanChFFT);

                this.timeTable.fftShift = tic;
                this.res.unFittedFFTShift = fftshift(this.res.unFittedFFT, 1);
                this.timeTable.fftShift = toc(this.timeTable.fftShift);

                this.timeTable.createFit = tic;                  
                this.res.fitModel = this.calcFittedFFT();
                this.timeTable.createFit = toc(this.timeTable.createFit);

                this.timeTable.subFit = tic;
                if this.general.fitToShift
                    this.res.fittedFFT = this.res.unFittedFFTShift - this.res.fitModel;
                else
                    this.res.fittedFFT = this.res.unFittedFFT - this.res.fitModel;
                end
                this.timeTable.subFit = toc(this.timeTable.subFit);

                this.timeTable.phiExtract = tic;
                this.res.phi = sqrt(abs(this.res.fittedFFT(this.freq.fUSIdx, :)));
                this.timeTable.phiExtract = toc(this.timeTable.phiExtract);
%             end
        end        
                
        function res = analyseSplittedData(this, resArr, splitNum)
            % usCompCmplx(output)      - [quants x ch x numOfPos]
            % qAvgFFT(output)          - [ch x samplesPerPos x numOfPos]
            % unFittedFFT(output)      - [samplesPerPos x numOfPos]
            % unFittedFFTShift(output) - [samplesPerPos x numOfPos]
            % fitModel(output)         - [samplesPerPos x numOfPos]
            % fittedFFT(output)        - [samplesPerPos x numOfPos]
            % phi(output)              - [1 x numOfPos]
            
            fprintf("ALGO: Analyzing splitted data.\n");
            this.res = [];
            splitQAVGFFT = zeros(this.uVars.channels, this.samples.samplesPerPos, this.samples.numOfPos, splitNum);
            
            this.timeTable.meanQuantFFT = tic;
            for i = 1:splitNum
                splitQAVGFFT(:,:,:,i) = resArr(i).qAvgChFFT;
            end
            this.res.qAvgChFFT          = mean(splitQAVGFFT, 4);
            this.timeTable.meanQuantFFT = toc(this.timeTable.meanQuantFFT);
            
            this.timeTable.meanChFFT = tic;
            this.res.unFittedFFT     = squeeze(mean(this.res.qAvgChFFT, 1));
            this.timeTable.meanChFFT = toc(this.timeTable.meanChFFT);
            
            this.timeTable.fftShift   = tic;
            this.res.unFittedFFTShift = fftshift(this.res.unFittedFFT, 1);
            this.timeTable.fftShift   = toc(this.timeTable.fftShift);
            
            this.timeTable.createFit = tic;  
            this.res.fitModel        = this.calcFittedFFT();
            this.timeTable.createFit = toc(this.timeTable.createFit);
            
            this.timeTable.subFit = tic;
            if this.general.fitToShift
                this.res.fittedFFT = this.res.unFittedFFTShift - this.res.fitModel;
            else
                this.res.fittedFFT = this.res.unFittedFFT - this.res.fitModel;
            end
            this.timeTable.subFit = toc(this.timeTable.subFit);
            
            this.timeTable.phiExtract = tic;
            this.res.phi = sqrt(abs(this.res.fittedFFT(this.freq.fUSIdx, :)));
            this.timeTable.phiExtract = toc(this.timeTable.phiExtract);
            
            res = this.res;
        end
        
        function fittedData = calcFittedFFT(this)
            fittedData = zeros(size(this.res.unFittedFFT));
            
            if this.general.useGPU
                fittedData = gpuArray(fittedData);
            end 
            
            if this.general.fitToShift
                fBar = this.freq.frequencyBarShifted;
                idxs = this.freq.fitIdxShift;
                dataToFit = double(this.res.unFittedFFTShift);
            else
                fBar = this.freq.frequencyBar;
                idxs = this.freq.fitIdx;
                dataToFit = double(this.res.unFittedFFT);
            end
            
            for i = 1:size(fittedData, 2)
                fitModel = fit(fBar(idxs)', dataToFit(idxs,i), this.freq.linearft);
                    
                fittedData(:, i) = fitModel.a + ...
                    fitModel.b*fBar + ...
                    fitModel.c*(1/2)*(3*fBar.^2-1) + ...
                    fitModel.d*(1/2)*(5*fBar.^3-3*fBar) + ...
                    fitModel.e*(1/8)*(35*fBar.^4-30*fBar.^2+3) + ...
                    fitModel.f*(1/8)*(63*fBar.^5-70*fBar.^3+15*fBar)+ ...
                    fitModel.g*(1/16)*(231*fBar.^6-315*fBar.^4+105*fBar.^2-5) + ...
                    fitModel.h*(1/16)*(429*fBar.^7-693*fBar.^5+315*fBar.^3-35*fBar) + ...
                    fitModel.k*(1/128)*(6435*fBar.^8 - 12012*fBar.^6 + 6930*fBar.^4 - 1260*fBar.^2 + 35) + ...
                    fitModel.l*(1/128)*(12155*fBar.^9 - 25740*fBar.^7 + 18018*fBar.^5 - 4620*fBar.^3 + 315*fBar) + ...
                    fitModel.m*(1/256)*(46189*fBar.^10 - 109395*fBar.^8 + 90090*fBar.^6 - 30030*fBar.^4 + 3465*fBar.^2 - 63);
            end
        end
        
        function timeTable = getTimeTable(this)
            timeTable = this.timeTable;
        end
        
        function vars = getAlgoVars(this)
            vars.extClk    = this.extClk;
            vars.usSignal  = this.usSignal;
            vars.geometry  = this.geometry;
            vars.digitizer = this.digitizer;
            vars.samples   = this.samples;
            vars.timing    = this.timing;
            vars.len       = this.len;
            vars.freq      = this.freq;
        end
        
        function resetAlgoArrays(this)
            this.res=struct();
        end
        
        function setRawData (this, rawData)
            this.data = rawData;
        end
        
        function res = analyseLoadedData(this, resArr, splitNum)
            % usCompCmplx(output)      - [quants x ch x numOfPos]
            % qAvgFFT(output)          - [ch x samplesPerPos x numOfPos]
            % unFittedFFT(output)      - [samplesPerPos x numOfPos]
            % unFittedFFTShift(output) - [samplesPerPos x numOfPos]
            % fitModel(output)         - [samplesPerPos x numOfPos]
            % fittedFFT(output)        - [samplesPerPos x numOfPos]
            % phi(output)              - [1 x numOfPos]
            
            fprintf("ALGO: Analyzing splitted data.\n");
            this.res = [];
            splitQAVGFFT = zeros(this.uVars.channels, this.samples.samplesPerPos, this.samples.numOfPos, splitNum);
            
            this.timeTable.meanQuantFFT = tic;
            for i = 1:splitNum
                splitQAVGFFT(:,:,:,i) = resArr(i).qAvgChFFT;
            end
            this.res.qAvgChFFT          = mean(splitQAVGFFT, 4);
            this.timeTable.meanQuantFFT = toc(this.timeTable.meanQuantFFT);
            
            this.res.fitModel        = this.calcFittedFFTNew();
            
            this.timeTable.subFit = tic;
            if this.general.fitToShift
                this.res.fittedFFT = fftshift(this.res.qAvgChFFT, 2) - this.res.fitModel;
            else
                this.res.fittedFFT = this.res.qAvgChFFT - this.res.fitModel;
            end
            this.timeTable.subFit = toc(this.timeTable.subFit);
            
            this.res.fittedMeanFFT = squeeze(mean( this.res.fittedFFT,1));
            
            this.timeTable.phiExtract = tic;
            this.res.phi = sqrt(abs(this.res.fittedMeanFFT(this.freq.fUSIdx, :)));
            this.timeTable.phiExtract = toc(this.timeTable.phiExtract);
            
            res = this.res;
        end
        
        function fittedData = calcFittedFFTNew(this)
            fittedData = zeros(size(this.res.qAvgChFFT));
            
            if this.general.useGPU
                fittedData = gpuArray(fittedData);
            end 
            
            if this.general.fitToShift
                fBar = this.freq.frequencyBarShifted;
                idxs = this.freq.fitIdxShift;
                dataToFit = double(fftshift(this.res.qAvgChFFT, 2));
            else
                fBar = this.freq.frequencyBar;
                idxs = this.freq.fitIdx;
                dataToFit = double(this.res.unFittedFFT);
            end
            
            for i = 1:size(fittedData, 1)
                for j = 1:size(fittedData, 3)
                    fitModel = fit(fBar(idxs)', dataToFit(i, idxs, j)', this.freq.linearft);

                    fittedData(i, :, j) = fitModel.a + ...
                        fitModel.b*fBar + ...
                        fitModel.c*(1/2)*(3*fBar.^2-1) + ...
                        fitModel.d*(1/2)*(5*fBar.^3-3*fBar) + ...
                        fitModel.e*(1/8)*(35*fBar.^4-30*fBar.^2+3) + ...
                        fitModel.f*(1/8)*(63*fBar.^5-70*fBar.^3+15*fBar)+ ...
                        fitModel.g*(1/16)*(231*fBar.^6-315*fBar.^4+105*fBar.^2-5) + ...
                        fitModel.h*(1/16)*(429*fBar.^7-693*fBar.^5+315*fBar.^3-35*fBar) + ...
                        fitModel.k*(1/128)*(6435*fBar.^8 - 12012*fBar.^6 + 6930*fBar.^4 - 1260*fBar.^2 + 35) + ...
                        fitModel.l*(1/128)*(12155*fBar.^9 - 25740*fBar.^7 + 18018*fBar.^5 - 4620*fBar.^3 + 315*fBar) + ...
                        fitModel.m*(1/256)*(46189*fBar.^10 - 109395*fBar.^8 + 90090*fBar.^6 - 30030*fBar.^4 + 3465*fBar.^2 - 63);
                end
            end
        end
        
    end
end
