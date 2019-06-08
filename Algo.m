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
        
        graphics;
        
        timeTable;
        res;
    end
 
    methods (Static)   
        function gNames = getGraphicsNames()
            gNames = {'extClk';'usSignal';'fullSignal';'measSamples';'netSignal';'reshapedSignal';'FFT';'phiCh';'phi'};
        end
        
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
            uVars.useGPU            = false;
            uVars.exportRawData     = true;
            uVars.gReq              = Algo.createGraphicRequest();
        end
       
        function gReq = createGraphicRequest()
            gReq = Graphics.createGraphicsRunVars(Algo.getGraphicsNames());
        end
        
    end
    
    methods
        function this = Algo()
            this.initTimeTable();
            
            this.graphics.graphicsNames = this.getGraphicsNames();
            this.graphics.obj  = Graphics(this.graphics.graphicsNames);
            this.graphics.gReq = Graphics.createGraphicsRunVars(this.graphics.graphicsNames);
            this.setGraphicsStaticVars();
        end
        
        function initTimeTable(this)
            this.timeTable.copyFullRawData = 0;

            this.timeTable.Reshape1 = 0;
            this.timeTable.Reshape2 = 0;
            this.timeTable.OverallReshape = 0;
            
            this.timeTable.CopyNetSignal1 = 0;
            this.timeTable.CutDigiPreSamples = 0;
            this.timeTable.CutPropPreSamples = 0;
            this.timeTable.CutPostSamples = 0;
            this.timeTable.ExtractNetSignal = 0;
            
            this.timeTable.FFT = 0;
            this.timeTable.PowerSpectrunm = 0;
            this.timeTable.ChAvg = 0;
            this.timeTable.SignalProcessing = 0;
            
            this.timeTable.FullAnalysis = 0;
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
           this.uVars.useGPU            = user.useGPU;
           this.uVars.exportRawData     = user.exportRawData;
           vars = this.calcDimensions();
           
           this.graphics.gReq = user.gReq;
           this.graphics.obj.setGlobalReq(user.gReq);
           this.setGraphicsDynamicVars()
           
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
            
            vars = this.getVars();
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
            vars.uVars      = this.uVars;
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
        end
        
        function calcSamples(this)
            this.samples.samplesPerSin         = this.digitizer.fs / this.usSignal.fSin;
            this.samples.samplesPerPulse       = this.samples.samplesPerSin * this.usSignal.cycInPulse;
            this.samples.samplesPerTrain       = this.digitizer.fs / this.usSignal.fTrain;
            this.samples.numOfTrains           = ceil(this.uVars.timeToSample * this.usSignal.fTrain);
            this.samples.prePhantomSamples     = ceil( (this.geometry.distFromPhantom / this.geometry.c) * this.digitizer.fs );
            this.samples.samplesPerSignal      = this.samples.samplesPerTrain * this.samples.numOfTrains;
            this.samples.samplesPerMeas        = this.uVars.preTriggerSamples + this.samples.prePhantomSamples + this.samples.samplesPerSignal;
            this.samples.samplesPerSignalAllCh = this.samples.samplesPerSignal * this.uVars.channels;
            this.samples.samplesPerMeasAllCh   = this.samples.samplesPerMeas   * this.uVars.channels;

            this.samples.samplesPerPos   = this.samples.samplesPerPulse * this.samples.numOfTrains;
            this.samples.numOfPos        = this.samples.samplesPerTrain / this.samples.samplesPerPulse;
            this.samples.samplesPerZAxis = this.samples.numOfPos * this.samples.samplesPerPulse;
        end
        
        function calcTiming(this)
            this.timing.tSin  = 1/this.usSignal.fSin;                         %[s], duration of one period of the sin
            this.timing.tPulse = this.usSignal.cycInPulse * this.timing.tSin; %[s], duration of the pulse
            this.timing.tTrain = 1/this.usSignal.fTrain;                       %[s], duration of one period of the train
            this.timing.dts    = 1/this.digitizer.fs;                          %[s], duration of a single sample

            this.timing.timeOfSample      = this.timing.dts*this.samples.samplesPerSignal; %  should it be samplesPerSignal -1?
            this.timing.actualSampledTime = this.digitizer.samplesPerAcq * this.timing.dts;
            
            this.timing.tPulseVec  = (0:1:this.samples.samplesPerPulse  - 1)*this.timing.dts;
            this.timing.tTrainMeas = (0:1:this.samples.samplesPerTrain  - 1)*this.timing.dts;
            this.timing.tSigVec    = (0:1:this.samples.samplesPerSignal - 1)*this.timing.dts;
            this.timing.tMeasVec   = (0:1:this.samples.samplesPerMeas   - 1)*this.timing.dts;
            this.timing.tPosVec    = (0:1:this.samples.samplesPerPos    - 1)*this.timing.dts;
            this.timing.tAcqVec    = (0:1:this.digitizer.samplesPerAcq  - 1)*this.timing.dts;
            
            
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
            this.freq.fSinIdx             = (this.usSignal.fSin / this.digitizer.fs) * this.samples.samplesPerPos + 1;        
        end
        
        function [sigData, clkData] = createSignalsForfGen(this)
            % Generate the Signal
            dt = 1/this.extClk.fSclk;
            tSig = (0:1:(this.usSignal.SclkSamplesInTrain-1)) * dt;

            sigData = zeros(this.usSignal.SclkSamplesInTrain, 1);
            sigData(1:this.usSignal.SclkSamplesInPulse) = sin(2*pi*this.usSignal.fSin*tSig(1:this.usSignal.SclkSamplesInPulse));
            
            %Generate The Clock
            clkData      = ones(this.extClk.extClkSamplesPerCyc, 1);
            dutyCycleIdx = floor(this.extClk.extClkSamplesPerCyc*(this.extClk.extClkDcyc/100))+1;
            clkData(dutyCycleIdx:end) = 0;
            clkData = repmat(clkData, this.extClk.extClkCycles, 1);
            
            tExtClk = (0:1:(this.extClk.extClkSigSamples-1)) * dt;
            
            this.extClk.data = clkData;
            this.usSignal.data = sigData;
          
            if this.graphics.gReq.validStruct.extClk
                this.graphics.obj.plotSingleSignal('extClk', tExtClk*1e6, clkData);
            end
            
            if this.graphics.gReq.validStruct.usSignal
                this.graphics.obj.plotSingleSignal('usSignal', tSig*1e6, sigData);
            end
            
        end      
        
        function res = analyse(this, data)
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
            if this.uVars.exportRawData
                this.timeTable.copyFullRawData = tic;
                this.res.rawData = gather(data);
                this.timeTable.copyFullRawData = toc(this.timeTable.copyFullRawData);
            end
            
            this.timeTable.ExtractNetSignal = tic;
            data = this.extractNetSignal(data);
            this.timeTable.ExtractNetSignal = toc(this.timeTable.ExtractNetSignal);
            
            this.timeTable.OverallReshape = tic;
            data = this.reshapeSignal(data);
            this.timeTable.OverallReshape = toc(this.timeTable.OverallReshape);
            
            this.timeTable.SignalProcessing = tic;
            this.signalProcessing(data);
            this.timeTable.SignalProcessing = toc(this.timeTable.SignalProcessing);
            
            res = this.res;
            
            this.timeTable.FullAnalysis = toc(this.timeTable.FullAnalysis);
            
        end
        
        function data = extractNetSignal(this, data)
            % rawData(input)  - [ch x samplesPerAcq]
            % rawData(output) - [ch x samplesPerSignal]
            
            if this.graphics.gReq.validStruct.fullSignal
                this.graphics.obj.plotSelectedSignals('fullSignal', this.timing.tAcqVec*1e6, data);
            end
            
            if this.graphics.gReq.validStruct.measSamples
                this.graphics.obj.plotSelectedSignals('measSamples', this.timing.tMeasVec*1e6, data(:, 1:this.samples.samplesPerMeas));
            end
           
            % Chop unneccesarry cycles
            this.timeTable.CutDigiPreSamples = tic;
            data(:,1:this.digitizer.preTriggerSamples)    = []; % Card Bug 
            this.timeTable.CutDigiPreSamples = toc(this.timeTable.CutDigiPreSamples);
            
            this.timeTable.CutPropPreSamples = tic;
            data(:,1:this.samples.prePhantomSamples)      = []; % Spread untill first pulse hitting the phantom
            this.timeTable.CutPropPreSamples = toc(this.timeTable.CutPropPreSamples);
            
            this.timeTable.CutPostSamples = tic;
            data(:,(this.samples.samplesPerSignal+1):end) = []; % extra trains to fill a buffer
            this.timeTable.CutPostSamples = toc(this.timeTable.CutPostSamples);
            
            if this.graphics.gReq.validStruct.netSignal
                this.graphics.obj.plotSelectedSignals('netSignal', this.timing.tSigVec*1e6, data);
            end
            
            if this.uVars.exportRawData
                this.timeTable.CopyNetSignal1 = tic;
                this.res.netSignal = gather(data);
                this.timeTable.CopyNetSignal1 = toc(this.timeTable.CopyNetSignal1);
            end
        end
        
        function reshapedData = reshapeSignal(this,data)
           % signal(input)  - [ch x samplesPerSignal]
           % signal(output) - [ch x samplesPerPos x numOfPos]
               
            this.timeTable.Reshape1 = tic;
            data = reshape(data, this.digitizer.channels,...
                                              this.samples.samplesPerTrain,...
                                              this.samples.numOfTrains);
            this.timeTable.Reshape1 = toc(this.timeTable.Reshape1);
            
            this.timeTable.allocReshapeMem = tic;
            if this.uVars.useGPU
                reshapedData = gpuArray(zeros(this.digitizer.channels, this.samples.samplesPerPos, this.samples.numOfPos));
            else
                reshapedData = zeros(this.digitizer.channels, this.samples.samplesPerPos, this.samples.numOfPos);
            end
            this.timeTable.allocReshapeMem = toc(this.timeTable.allocReshapeMem);
            
            this.timeTable.Reshape2 = tic;   
            for i=1:this.samples.numOfPos
                reshapedData(:,:,i) = reshape( data(:, ((i-1)*this.samples.samplesPerPulse+1) : (i*this.samples.samplesPerPulse), :),...
                             this.digitizer.channels,...
                             this.samples.samplesPerPos,...
                             1);

            end
            this.timeTable.Reshape2 = toc(this.timeTable.Reshape2);
            
            if this.uVars.exportRawData 
                this.timeTable.copyReshaped = tic;
                this.res.reshapedSignal = gather(reshapedData);
                this.timeTable.copyReshaped = toc(this.timeTable.copyReshaped);
            end
           
            if this.graphics.gReq.validStruct.reshapedSignal
                this.graphics.obj.plotSelectedSignals('reshapedSignal', this.timing.tPosVec*1e6, reshapedData);
            end

        end
        
        function signalProcessing(this, data)
            % signal(input)  - [ch x samplesPerPos x numOfPos]
            % fftRes(output) - [ch x samplesPerPos x numOfPos]
            % phiCh(output)  - [ch x numOfPos] 
            % phi(output)    - [numOfPos x 1]
            
            this.timeTable.FFT = tic;
            data = (2./this.samples.samplesPerPos)*fft(data,[],2);
            this.timeTable.FFT = toc(this.timeTable.FFT); 
            
            if this.uVars.exportRawData
                this.timeTable.copyFFT = tic;
                this.res.fftRes = gather(data);
                this.timeTable.copyFFT = toc(this.timeTable.copyFFT);
            end
            
            this.timeTable.PowerSpectrunm = tic;
            this.res.phiCh  = gather(squeeze(abs(data(:,this.freq.fSinIdx,:))));
            this.timeTable.PowerSpectrunm = toc(this.timeTable.PowerSpectrunm); 
            
            this.timeTable.ChAvg = tic;
            this.res.phi = permute(mean(this.res.phiCh,1), [2,1]);    
            this.timeTable.ChAvg = toc(this.timeTable.ChAvg);   
            
                       
            if this.graphics.gReq.validStruct.FFT
                this.graphics.obj.plotFFT('FFT', this.freq.frequencyBarShifted*1e-6, abs(data),...
                    this.freq.frequencyBar(this.freq.fSinIdx)*1e-6, abs(data(:, this.freq.fSinIdx, :)));
            end
            
            if this.graphics.gReq.validStruct.phiCh
                this.graphics.obj.plotMultipleSignals('phiCh', this.len.zVecUSRes*1e3, this.res.phiCh) 
            end
            
            if this.graphics.gReq.validStruct.phi
                this.graphics.obj.plotSingleSignal('phi', this.len.zVecUSRes*1e3, this.res.phi);
            end
            
            
%             if this.plotReq.FFT.getIsValid() && isgraphics(this.plotReq.FFT.getAx())
%                 [~, ch, pos, ~] = this.plotReq.FFT.getRequest();
%                 this.displayResults(this.plotReq.FFT, this.freq.frequencyBarShifted*1e6,...
%                                     abs(this.res.fftRes(ch, :, pos))',...
%                                     sprintf('FFT, ch: %d, pos idx: %d', ch, pos),...
%                                     'f[MHz]', 'Power Spectrum');
%             end

%             if this.plotReq.phiCh.getIsValid() && isgraphics(this.plotReq.phiCh.getAx())
%                legstr = cell(1,this.digitizer.channels);
%                for i=1:this.digitizer.channels
%                     legstr{i} = sprintf('Ch: %d', i);
%                end
%                this.displayResults(this.plotReq.phiCh, this.len.zVecUSRes*1e-3,...
%                                    this.res.phiCh', 'Phi Per Channel', 'z[mm]', 'Fluence Rate', legstr);
%             end
%            
%             if this.plotReq.phi.getIsValid() && isgraphics(this.plotReq.phi.getAx())
%                this.displayResults(this.plotReq.phi, this.len.zVecUSRes*1e-3, this.res.phi,...
%                                    'Phi', 'z[mm]', 'Fluence Rate');
%             end
           
        end
        
        function timeTable = getTimeTable(this)
            timeTable = this.timeTable;
        end
        
        function vars = getAlgoVars(this)
            vars.extClk  = this.extClk;
            vars.usSignal = this.usSignal;
            vars.geometry = this.geometry;
            vars.digitizer = this.digitizer;
            vars.samples = this.samples;
            vars.timing = this.timing;
            vars.len = this.len;
            vars.freq = this.freq;
        end
        
        function resetAlgoArrays(this)
            this.res=struct();
        end
        
        function setGraphicsStaticVars(this)
%             this.graphics.setGeneral(internal, valid, type);
%             this.graphics.setStrings(title, xlabel, ylabel, legend);
%             this.graphics.setValues(pos, ch, posDim, chDim, dataDim);
%             this.graphics.setMarkers(plotMark, markX, markY);
%             this.graphics.setLims(xlims, ylims, colors);
%             this.graphics.setFonts(type, titleSize, labelsSize, axisSize);
            
            % extClk
            this.graphics.obj.setType(this.graphics.graphicsNames{1}, 'stem');
            this.graphics.obj.setStrings(this.graphics.graphicsNames{1}, "External Clk (Sampling Clock)", "t[\mus]", "Amp", []);
            this.graphics.obj.setDims(this.graphics.graphicsNames{1}, [], [], 1);
            
            % usSignal
            this.graphics.obj.setType(this.graphics.graphicsNames{2}, 'stem');
            this.graphics.obj.setStrings(this.graphics.graphicsNames{2}, "Ultrasound Pulse", "t[\mus]", "Amp", []);
            this.graphics.obj.setDims(this.graphics.graphicsNames{2}, [], [], 1);
           
            % fullSignal [ch x samplesPerAcq]
            this.graphics.obj.setType(this.graphics.graphicsNames{3}, 'stem');
            this.graphics.obj.setStrings(this.graphics.graphicsNames{3}, "Full Signal Channel: %d", "t[\mus]", "Amp[V]", []);
            this.graphics.obj.setDims(this.graphics.graphicsNames{3}, [], 1, 2);

            % measSamples [ch x samplesPerMeas]
            this.graphics.obj.setType(this.graphics.graphicsNames{4}, 'stem');
            this.graphics.obj.setStrings(this.graphics.graphicsNames{4}, "Measured Signal Channel: %d", "t[\mus]", "Amp[V]", []);
            this.graphics.obj.setDims(this.graphics.graphicsNames{4}, [], 1, 2);
            
            % netSignal [ch x samplesPerSignal]
            this.graphics.obj.setType(this.graphics.graphicsNames{5}, 'stem');
            this.graphics.obj.setStrings(this.graphics.graphicsNames{5}, "Net Signal Channel: %d", "t[\mus]", "Amp[V]", []);
            this.graphics.obj.setDims(this.graphics.graphicsNames{5}, [], 1, 2);

            % reshapedSignal [ch x samplesPerPos x numOfPos]
            this.graphics.obj.setType(this.graphics.graphicsNames{6}, 'stem');
            this.graphics.obj.setStrings(this.graphics.graphicsNames{6}, "Reshaped Signal Channel: %d, Pos: %.2f[mm] (idx:(%d))", "t[\mu s]", "Amp[V]", []);
            this.graphics.obj.setDims(this.graphics.graphicsNames{6}, 3, 1, 2);
                        
            %FFT [ch x samplesPerPos x numOfPos]
            this.graphics.obj.setType(this.graphics.graphicsNames{7}, 'plot');
            this.graphics.obj.setStrings(this.graphics.graphicsNames{7}, "FFT Pos: %.2f[mm] (idx:(%d))", "f[MHz]", "Power Spectrum", "Ch %d");
            this.graphics.obj.setDims(this.graphics.graphicsNames{7}, 3, 1, 2);
            this.graphics.obj.setMarkersEnable(this.graphics.graphicsNames{7}, true);
            
            % phiCh [ch x numOfPos]
            this.graphics.obj.setType(this.graphics.graphicsNames{8}, 'stem');
            this.graphics.obj.setStrings(this.graphics.graphicsNames{8}, "\\phi_{ch}", "z[mm]", "\phi", "Ch %d");
            this.graphics.obj.setDims(this.graphics.graphicsNames{8}, 2, 1, 2);
            
            % phi [numOfPos]
            this.graphics.obj.setType(this.graphics.graphicsNames{9}, 'stem');
            this.graphics.obj.setStrings(this.graphics.graphicsNames{9}, "\\phi", "z[mm]", "\phi", []);
            this.graphics.obj.setDims(this.graphics.graphicsNames{1}, 1, [], 1);

        end
        
        function setGraphicsDynamicVars(this)
            this.graphics.obj.setChAndPos(this.graphics.gReq.ch, this.graphics.gReq.pos)
            
            z = this.len.zVecUSRes(this.graphics.gReq.ch)*1e3;
            ch = this.digitizer.channels;
            fSin = this.usSignal.fSin;
            
            
            this.graphics.obj.setTitleVariables(this.graphics.graphicsNames{1}, []);
            this.graphics.obj.setTitleVariables(this.graphics.graphicsNames{2}, []);
            this.graphics.obj.setTitleVariables(this.graphics.graphicsNames{3}, [this.graphics.gReq.ch]);              % fullSignal 
            this.graphics.obj.setTitleVariables(this.graphics.graphicsNames{4}, [this.graphics.gReq.ch]);              % measSamples
            this.graphics.obj.setTitleVariables(this.graphics.graphicsNames{5}, [this.graphics.gReq.ch]);              % netSignal
            this.graphics.obj.setTitleVariables(this.graphics.graphicsNames{6}, [this.graphics.gReq.ch, z, this.graphics.gReq.pos]); % reshapedSignal
            this.graphics.obj.setTitleVariables(this.graphics.graphicsNames{7}, [z, this.graphics.gReq.pos]);          % FFT
            this.graphics.obj.setTitleVariables(this.graphics.graphicsNames{8}, []);
            this.graphics.obj.setTitleVariables(this.graphics.graphicsNames{9}, []);
            
            this.graphics.obj.setLegendVariables(this.graphics.graphicsNames{7}, 1:ch);
            this.graphics.obj.setLegendVariables(this.graphics.graphicsNames{8}, 1:ch);
             
            this.graphics.obj.setLimits(this.graphics.graphicsNames{7}, (fSin + fSin*[-0.01, 0.01])*1e-6, []);
            
            this.graphics.obj.updateGraphicsConstruction()
        end
    end
end
