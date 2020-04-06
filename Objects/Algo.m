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
        
        data
        reshapedData
        Wn
        
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
            
            uVars.useQuant          = false;
            uVars.quantTime         = [];
            uVars.fastAnalysis      = false;
            uVars.useGPU            = false;
            uVars.useHadamard       = false;
            
            % Do i need this? should be tested on a PC with GPU
            % see if gather takes long time.
%             uVars.exportRawData.rawData       = false;
%             uVars.exportRawData.netSignal     = false;
%             uVars.exportRawData.deMultiplexed = false;
%             uVars.exportRawData.reshape       = false;
%             uVars.exportRawData.FFT           = false;

            %uVars.gReq              = Algo.createGraphicRequest();
        end
       
        function gReq = createGraphicRequest()
            gReq = algoGraphics.createGraphicsRunVars();
        end
        
    end
    
    methods
        function this = Algo()
            this.initTimeTable();          
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
            
            this.timeTable.deMultiplex              = 0;
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
            this.timeTable.meanFFT                  = 0;
            this.timeTable.extractFreq              = 0;
            this.timeTable.squeezeCmplxPhi          = 0;
            this.timeTable.absPhi                   = 0;
            this.timeTable.PowerSpectrunm           = 0;
            this.timeTable.ChAvg                    = 0;
            this.timeTable.QuantAvg                 = 0;
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
           this.uVars.useQuant          = user.useQuant;
           this.uVars.quantTime         = user.quantTime;
           this.uVars.fastAnalysis      = user.fastAnalysis;
           this.uVars.useGPU            = user.useGPU;
%            this.uVars.exportRawData     = user.exportRawData;
           this.uVars.useHadamard       = user.useHadamard;
           
           vars = this.calcDimensions();
           
%            this.graphics.gReq = user.gReq;
%            this.graphics.obj.setGlobalReq(user.gReq);
%            this.setGraphicsDynamicVars();
           
        end
        
        function vars = calcDimensions(this)
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
                
                this.samples.numOfTrains     = ceil(this.uVars.timeToSample * this.usSignal.fTrain);
                this.samples.trainsPerQuant  = ceil(this.uVars.quantTime *this.usSignal.fTrain);
                this.samples.numOfTrains     = ceil(this.samples.numOfTrains/this.samples.trainsPerQuant)*this.samples.trainsPerQuant;
                this.samples.numOfQuant      = this.samples.numOfTrains / this.samples.trainsPerQuant;
                
                this.samples.samplesPerQuant       = this.samples.trainsPerQuant*this.samples.samplesPerTrain;
                this.samples.samplesPerPos         = this.samples.samplesPerPulse * this.samples.trainsPerQuant;
                this.samples.prePhantomSamples     = ceil( (this.geometry.distFromPhantom / this.geometry.c) * this.digitizer.fs );
                this.samples.inPhantomPropSamples  = this.samples.samplesPerTrain - this.samples.samplesPerPulse; 
                this.samples.samplesPerSignal      = this.samples.samplesPerTrain * this.samples.numOfTrains;
                if this.uVars.useHadamard
                    this.samples.samplesPerMeas        = this.uVars.preTriggerSamples + this.samples.prePhantomSamples + this.samples.inPhantomPropSamples + this.samples.samplesPerSignal;
                else
                    this.samples.samplesPerMeas        = this.uVars.preTriggerSamples + this.samples.prePhantomSamples + this.samples.samplesPerSignal;
                end
                this.samples.samplesPerSignalAllCh = this.samples.samplesPerSignal * this.uVars.channels;
                this.samples.samplesPerMeasAllCh   = this.samples.samplesPerMeas   * this.uVars.channels;
            
                this.samples.numOfPos        = this.samples.samplesPerTrain / this.samples.samplesPerPulse;
                this.samples.samplesPerZAxis = this.samples.numOfPos * this.samples.samplesPerPulse;
                
                if this.uVars.useHadamard
                    A = speye(this.samples.numOfTrains);
%                     eye(this.samples.numOfTrains);
                    % create sparse matrix for single vec demultiplexing
                    % (single sVec)
                    this.hadamard.nSVecPad = this.samples.numOfPos * this.samples.samplesPerPulse;
                    this.hadamard.singleSVecDeMulMat = zeros(this.samples.numOfPos, this.hadamard.nSVecPad);
                    this.hadamard.singleSVecDeMulMat(:,1:(this.samples.samplesPerPulse):end) = this.hadamard.sMatInv;
                                        
                    % Expand the matrix to keep order of demultiplexing
                    this.hadamard.singleTrainMat = kron( flip(this.hadamard.singleSVecDeMulMat, 1), ones(  this.samples.samplesPerPulse,1) ); 
%                     this.hadamard.singleTrainMat = kron( this.hadamard.singleSVecDeMulMat, ones(  this.samples.samplesPerPulse,1) ); 
                    
                    % create sparse matrix for a single train
                    this.hadamard.singleTrainDeMulMat = zeros(this.hadamard.nSVecPad, this.hadamard.nSVecPad);
                    for i = 1:this.samples.samplesPerPulse
                        this.hadamard.singleTrainDeMulMat( i:this.samples.samplesPerPulse:end, :) = ...
                                                            circshift(this.hadamard.singleTrainMat( i:this.samples.samplesPerPulse:end, :), (i-1), 2);
                    end

                    % expand the matrix to demultiplex all of the trains
                    this.hadamard.invSMatFullSig = kron( sparse(A),...
                                                         sparse(this.hadamard.singleTrainDeMulMat));
                                                     
%                     this.hadamard.invSMatFullSig = kron( sparse( A ),...
%                                                          sparse( this.hadamard.singleTrainDeMulMat ) );
                    if this.uVars.useGPU
                        this.hadamard.invSMatFullSig = gpuArray(this.hadamard.invSMatFullSig);
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
            this.timing.tPulseVec  = (0:1:this.samples.samplesPerPulse  - 1)*this.timing.dts;
            this.timing.tTrainMeas = (0:1:this.samples.samplesPerTrain  - 1)*this.timing.dts;
            this.timing.tSigVec    = (0:1:this.samples.samplesPerSignal - 1)*this.timing.dts;
            this.timing.tMeasVec   = (0:1:this.samples.samplesPerMeas   - 1)*this.timing.dts;
            this.timing.tPosVec    = (0:1:this.samples.samplesPerPos    - 1)*this.timing.dts;
            this.timing.tAcqVec    = (0:1:this.digitizer.samplesPerAcq  - 1)*this.timing.dts;
            this.timing.tQuantVec  = (0:1:this.samples.samplesPerQuant  - 1)*this.timing.dts;
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
            this.freq.fSinIdx             = (this.usSignal.fSin / this.digitizer.fs) * this.samples.samplesPerPos + 1;
            this.freq.fs                  = this.extClk.fSclk;
            this.Wn = exp((-1j*(2*pi)/this.samples.samplesPerPos)*(this.freq.fSinIdx-1)*(0:(this.samples.samplesPerPos-1)))';
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
            this.timeTable.FullAnalysis = toc(this.timeTable.FullAnalysis);
            
        end
        
        function extractNetSignal(this)
            % rawData(input)  - [ch x samplesPerAcq]
            % rawData(output) - [ch x samplesPerSignal]
           
            % Chop unneccesarry cycles
            digiPreSamplesIdx = this.digitizer.preTriggerSamples+1;
            propPreSamplesIdx = this.samples.prePhantomSamples +1;
            propInPhantomIdx  = this.samples.inPhantomPropSamples+1;
            
            this.timeTable.CopyMeasSignal1 = tic;
            this.res.measSignal = this.data(:, 1:this.samples.samplesPerMeas);
            this.timeTable.CopyMeasSignal1 = toc(this.timeTable.CopyMeasSignal1);
            
            this.timeTable.CutDigiPreSamples = tic;
            this.data = this.data(:,digiPreSamplesIdx:end); % Card Bug 
            this.timeTable.CutDigiPreSamples = toc(this.timeTable.CutDigiPreSamples);
            
            this.timeTable.CutPropPreSamples = tic;
            this.data      = this.data(:,propPreSamplesIdx:end); % Spread untill first pulse hitting the phantom
            this.timeTable.CutPropPreSamples = toc(this.timeTable.CutPropPreSamples);
            
            if this.uVars.useHadamard
                this.timeTable.CutInPhantopmPropSamples = tic;
                this.data = this.data(:,propInPhantomIdx:end);
                this.timeTable.CutInPhantopmPropSamples = toc(this.timeTable.CutInPhantopmPropSamples);
            end
            
            this.timeTable.CutPostSamples = tic;
            this.data = this.data(:,1:this.samples.samplesPerSignal); % extra trains to fill a buffer
            this.timeTable.CutPostSamples = toc(this.timeTable.CutPostSamples);
            
            this.timeTable.CopyNetSignal = tic;
            this.res.netSignal = gather(this.data);
            this.timeTable.CopyNetSignal = toc(this.timeTable.CopyNetSignal);

        end
        
        function demultiplexSignal(this)
            
            this.timeTable.deMultiplex = tic;
            this.data = (this.hadamard.invSMatFullSig * double(this.data'))';
            this.timeTable.deMultiplexed = toc(this.timeTable.deMultiplex);
            
            this.timeTable.CopyDeMultiplexed = tic;
            this.res.deMultiplexed = gather(this.data);
            this.timeTable.CopyDeMultiplexed = toc(this.timeTable.CopyDeMultiplexed);
            
            
%             figure();
%             plot(this.data())
%             if this.graphics.gReq.validStruct.deMul
%                 this.graphics.obj.plotSignal('deMul', this.timing.tSigVec*1e6, this.data) 
%             end
            
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
                this.res.reshapedSignal = gather(this.data);
                this.timeTable.copyReshaped = toc(this.timeTable.copyReshaped);
        end
        
        function signalProcessing(this)
            % signal(input)      - [quants x ch x samplesPerPos x numOfPos]
            % fftRes(output)     - [quants x ch x samplesPerPos x numOfPos]
            % qAvgFFT(output)    - [1 x ch x samplesPerPos x numOfPos]
            % phiChCmplx(output) - [quants x ch x numOfPos]
            % phiCh(output)      - [quants x ch x numOfPos]
            % phiQuant(output)   - [quant x numOfPos]
            % phi(output)        - [numOfPos]
            
            if this.uVars.fastAnalysis
                this.timeTable.totalFastAnalysis = tic;
                for k = 1:this.samples.numOfQuant
                    for j = 1:this.digitizer.channels
                        for i=1:this.samples.numOfPos
                        a = tic;
                        this.res.phiCh(k,j,i) = (2./this.samples.samplesPerPos) * abs(reshape(this.data(k,j,:,i), 1, this.samples.samplesPerPos, []) * this.Wn);
                        this.timeTable.fastAnalysis(i) = toc(a);
                        end
                    end
                end
                this.timeTable.totalFastAnalysis = toc(this.timeTable.totalFastAnalysis); 
            else
                this.timeTable.FFT = tic;
                this.data = (2./this.samples.samplesPerPos) * fft(cast(this.data, 'single'),[],3); %casting is not needed, already single
                this.timeTable.FFT = toc(this.timeTable.FFT); 

                this.timeTable.copyFFT = tic;
                this.res.fftRes = gather(this.data);
                this.timeTable.copyFFT = toc(this.timeTable.copyFFT);
                
                this.timeTable.meanFFT = tic;
                this.res.qAvgFFT = mean(abs(this.res.fftRes), 1);
                this.timeTable.meanFFT = toc(this.timeTable.meanFFT);
                
                this.timeTable.PowerSpectrunm = tic;

                this.timeTable.extractFreq = tic;
                this.res.phiChCmplx = gather(this.data(:,:,this.freq.fSinIdx,:));
                this.timeTable.extractFreq = toc(this.timeTable.extractFreq);

                this.timeTable.squeezeCmplxPhi = tic;
                this.res.phiChCmplx = permute(this.res.phiChCmplx, [1,2,4,3]);
                this.timeTable.squeezeCmplxPhi = toc(this.timeTable.squeezeCmplxPhi);
                
                this.timeTable.absPhi = tic;
                this.res.phiCh = abs(this.res.phiChCmplx); %TODO: square?
                this.timeTable.absPhi = toc(this.timeTable.absPhi);

                this.timeTable.PowerSpectrunm = toc(this.timeTable.PowerSpectrunm);
            end

            this.timeTable.ChAvg = tic;
            this.res.phiQuant = gather(permute(mean(this.res.phiCh,2), [1,3,2]));    
            this.timeTable.ChAvg = toc(this.timeTable.ChAvg);   
            
            this.timeTable.QuantAvg = tic;
            this.res.phi    = gather(squeeze(mean(this.res.phiQuant,1)));
            this.res.phiStd = gather(squeeze(std(this.res.phiQuant, 0 ,1))); 
            this.timeTable.QuantAvg = toc(this.timeTable.QuantAvg);        
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
        
    end
end
