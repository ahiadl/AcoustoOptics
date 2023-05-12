classdef Algo < handle
    %ALGO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        uVars;
        
        general;
        geometry;
        
        sClk;
        usSignal;
        samples;
        digitizer;
        freq;
        timing;
        len;
        hadamard;

        timeTable;
        
        data;
        dataCal;
        emdData;
        measSamples;

        curRes;
        result;
        curSplit;

%         dataInSplit;
%         resSplit;
    end
 
    methods (Static)   
        function uVars = createUserVars()
            uVars.fSin              = [];              
            uVars.fSqnc             = [];
            uVars.cycPerPulse       = [];
            uVars.channels          = [];
            
            uVars.c                 = [];  
            uVars.distFromPhantom   = [];

            uVars.fs                = []; %fs
            
            uVars.fgClk             = [];
            uVars.sClkDcyc          = [];% [%]
            
            uVars.timeToSample      = [];
            uVars.frameTime         = [];

            uVars.envDC             = [];
            uVars.envUS             = [];
            uVars.useFrame          = false;

            uVars.useGPU              = false;
            uVars.useHadamard         = false;
            uVars.analyzeSingleCh     = false;
            uVars.calibrate           = false;
            uVars.acCoupling          = true;
            uVars.measTimeLimit       = [];

            uVars.cutArtfct           = false;
            uVars.artfctIdxVec        = [];
            
            uVars.calcMuEff     = false;
            uVars.muEffModel    = 'Uniform';
            uVars.muEffIdxs     = []; % After Cutting Artifacts

            uVars.extSNRDef    = false;
            uVars.extNoiseIdx  = [];
            uVars.extPeakIdx   = [];

            uVars.noiseStd      = []; % For virtual Data only

            uVars.export.meas        = false;
            uVars.export.signal      = false;
            uVars.export.deMul       = false;
            uVars.export.reshaped    = false;
            uVars.export.fft         = false;
            uVars.export.usCompCmplx = false;
        end
        
        function convertBaseChToSplits(baseCh)
            numReps = length(baseCh);
            numSplit =  length(baseCh.splitRes);
            
            
            dataInSplit = models.baseChSplitArr;

            for j=1:numReps
                for i = 1:numSplit
                    dataInSplit(j,:) = baseCh(j).splitRes(i);
                end
            end
        end
    
    end
    
    methods
        function this = Algo()
            this.initTimeTable();
            this.general.splitNum = 1;
            this.general.analysisReps = 1;
            models = this.getResultsModels();
            this.curRes = models.base;
            this.result = models.base;
        end

        %% Set/Get Functions
        function [vars] = setVars(this, user)
           this.uVars.fSin                  = user.fSin; %V
           this.uVars.fSqnc                 = user.fSqnc; %V
           this.uVars.cycPerPulse           = user.cycPerPulse; %V
           this.uVars.channels              = user.channels; %V
           this.uVars.timeToSample          = user.timeToSample;
           this.uVars.c                     = user.c;
           this.uVars.distFromPhantom       = user.distFromPhantom;
           this.uVars.fgClk                 = user.fgClk; %[S/s]
           this.uVars.fs                    = user.fs; %raw fs
           this.uVars.sClkDcyc              = user.sClkDcyc;
           this.uVars.envDC                 = user.envDC;
           this.uVars.envUS                 = user.envUS;
           this.uVars.useFrame              = user.useFrame;
           this.uVars.frameTime             = user.frameTime;
           this.uVars.useGPU                = user.useGPU;
           this.uVars.useHadamard           = user.useHadamard;
           this.uVars.analyzeSingleCh       = user.analyzeSingleCh;
           this.uVars.calibrate             = user.calibrate;
           this.uVars.acCoupling            = user.acCoupling;
           this.uVars.cutArtfct             = user.cutArtfct;
           this.uVars.artfctIdxVec          = user.artfctIdxVec;
           this.uVars.calcMuEff             = user.calcMuEff;
           this.uVars.muEffModel            = user.muEffModel;
           this.uVars.muEffIdxs             = user.muEffIdxs;
           this.uVars.extSNRDef             = user.extSNRDef;   
           this.uVars.extNoiseIdx           = user.extNoiseIdx;
           this.uVars.extPeakIdx            = user.extPeakIdx;    
           this.uVars.noiseStd              = user.noiseStd;
           this.uVars.export                = user.export;
           this.uVars.measTimeLimit         = user.measTimeLimit;

           this.calcParameters();
           vars = this.getVars();
        end
        
        function setRawData (this, measSamples)
            this.measSamples = measSamples;
            if this.general.analyzeSingleCh
                this.measSamples = gather(this.measSamples);
            end
        end
        
        function setCalibrationData(this,data)
            this.dataCal.raw = mean(data.frameAvgPowerFFT,3);
        end
        
        function vars = getVars(this)
            vars.uVars     = this.uVars;
            
            vars.general   = this.general;
            vars.geometry  = this.geometry;
            
            vars.sClk      = this.sClk;
            vars.usSignal  = this.usSignal;
            vars.samples   = this.samples;
            vars.digitizer = this.digitizer;
            vars.freq      = gather(this.freq);
            vars.timing    = this.timing;
            vars.len       = this.len;
            vars.hadamard  = this.hadamard;
        end
        
        function timeTable = getTimeTable(this)
            timeTable = this.timeTable;
        end
        
        %% Config Functions
        function calcParameters(this)
            this.calcGeneral();
            this.calcGeometry();
            
            this.calcSamplingClkDim();
            this.calcUsSignalDim();
            this.calcSamples();
            this.calcFreq();
            this.calcTiming();
            this.calcLength();

%             this.allocateReshapeMemory()
        end

        function calcGeneral(this)  
            channels            = this.uVars.channels; 

            useFrame        = this.uVars.useFrame;
            useGPU          = this.uVars.useGPU;
            useHadamard     = this.uVars.useHadamard;         
            analyzeSingleCh = this.uVars.analyzeSingleCh;
            calibrate       = this.uVars.calibrate;
            acCoupling      = this.uVars.acCoupling;
            cutArtfct       = this.uVars.cutArtfct;
            calcMuEff       = this.uVars.calcMuEff;
            muEffModel      = this.uVars.muEffModel;
            muEffIdxs       = this.uVars.muEffIdxs;
            extSNRDef       = this.uVars.extSNRDef ;   
            extNoiseIdx     = this.uVars.extNoiseIdx;
            extPeakIdx      = this.uVars.extPeakIdx;
            noiseStd        = this.uVars.noiseStd;

            %IndividualCh
            channelsToSample  = channels;
            channelsToAnalyze = channels;
            analysisReps      = 1;
            
            if analyzeSingleCh
                channelsToAnalyze = 1;
                analysisReps      = channels;
            end
            
            % Split Measurement
            timeToSampleTotal = this.uVars.timeToSample;
            measTimeLimit     = this.uVars.measTimeLimit;
            
            splitMeas = timeToSampleTotal > measTimeLimit;

            if splitMeas
                %Calculate time to sample
                splitTime = measTimeLimit+1;
                splitNum = 0;
                while splitTime > measTimeLimit
                    splitNum = splitNum+1;
                    splitTime = round(timeToSampleTotal/splitNum, 3);
                end
                timeToSample = splitTime;
                this.uVars.distFromPhantom = 0;
            else
                splitNum     = 1;
                timeToSample = timeToSampleTotal; 
            end
            
            % Export Data
            export.meas           = this.uVars.export.meas        && ~splitMeas;
            export.signal         = this.uVars.export.signal      && ~splitMeas;
            export.deMul          = this.uVars.export.deMul       && ~splitMeas && ~contHadamard;
            export.reshaped       = this.uVars.export.reshaped    && ~splitMeas && ~contSpeckleAnalysis;
            export.fft            = this.uVars.export.fft         && ~splitMeas;
            export.usCompCmplx    = this.uVars.export.usCompCmplx && ~splitMeas; %#ok<STRNU> 
            
            %--------------------
            % Collect parameters
            %--------------------
            clear ('measTimeLimit');
               
            vars = who();
            for i = 1:length(vars)
                if strcmp(vars{i}, 'this'); continue; end
                this.general.(vars{i}) = eval(vars{i});
            end
        end

        function calcGeometry(this)
            this.geometry.c = this.uVars.c;
            this.geometry.distFromPhantom = this.uVars.distFromPhantom;
        end

        function calcSamplingClkDim(this)
            % The sampling clock needs to by synchronized with the signal
            % fed into the US transducer. 
            % Therefore, we generate the sampling clk (sClk , fs) externally by the FG
            % rather than use the Digitizer internal sClk.
            % This function calculates the External sClock Params
            % to create with function generator.
            
            fs    = this.uVars.fs;
            fgClk = this.uVars.fgClk;
            
            % The data loaded to the AFG must be a multiplication of 16 in manners of data length.
            % creating sClk with 16 cycles meets this requirement no matter
            % the fFgClk used.
            sClkCycles  = 16;
            fsNaive     = fs;
            sClkDutyCyc = 50/100;
            
            fgClkSamplesPerSClkCyc = fgClk/fs;
            
            % If there is no round number of fgClk samples in sClk cycle, 
            % slow the sClk to the period where fgClk fits in.
            if mod(fgClkSamplesPerSClkCyc,1)~=0
                fgClkSamplesPerSClkCyc    =  ceil(fgClkSamplesPerSClkCyc);
                fs                   = fgClk / fgClkSamplesPerSClkCyc;
                fprintf("Notice: The sampling frequence you have chosen cannot be genrated by the AFG. \n the closest sampling frequency is: %d\n", fs);
            end

            sClkT                = 1/fs;
            fgClkSamplesPerFsSig = fgClkSamplesPerSClkCyc * sClkCycles;
            
            %Create sClk Data
            cycleData                   = ones(fgClkSamplesPerSClkCyc,1);
            dutyCycleIdx                = floor(fgClkSamplesPerSClkCyc*(1-sClkDutyCyc))+1;
            cycleData(dutyCycleIdx:end) = 0;
            clkData                     = repmat(cycleData, sClkCycles, 1);
            
            %--------------------
            % Collect parameters
            %--------------------
            vars = who();
            for i = 1:length(vars)
                if strcmp(vars{i}, 'this'); continue; end
                this.sClk.(vars{i}) = eval(vars{i});
            end
        end
        
        function calcUsSignalDim(this)
            % This function generates the US signal (single sequence) in manners of fgClk.
            % The US signal number of samples should be a multiply of:
            % Pulses, sClk, and 16 (fg demand).
            % This function is called after the function: calcSamplingClkDim
            % that makes sure that fgClk/sClk is a natural number.

            cycPerPulse = this.uVars.cycPerPulse;
            fSin        = this.uVars.fSin;
            fSqnc       = this.uVars.fSqnc;
            fgClk       = this.uVars.fgClk;
            
            % US Signal Vars
            fSqncNaive = fSqnc;
            
            fPulse = fSin/cycPerPulse; 
            
            fgClkSamplesPerSin   = fgClk / fSin;    % Hidden assumtion: fgClk/fSin is a Natural number.        
            fgClkSamplesPerPulse = fgClk / fPulse; 
            fgClkSamplesPerSqnc  = fgClk / fSqnc;   %this may vary later
            
            fgClkSamplesPerSClkCyc = this.sClk.fgClkSamplesPerSClkCyc;
            
            factor = lcm( lcm(fgClkSamplesPerPulse, fgClkSamplesPerSClkCyc), 16); %lcm = least common multiple
            
            fgClkSamplesPerSqnc = ceil(fgClkSamplesPerSqnc/factor)*factor;
            fSqnc               = fgClk / fgClkSamplesPerSqnc;
            
            pulsePerSqnc = fPulse / fSqnc;
            
            if this.uVars.useHadamard
                sMatrix = createSMatrix(pulsePerSqnc, 'QR'); 
                N       = size(sMatrix, 1);
                
                %createSMatrix returns an S-matrix with minimal order of N 
                %that satisfies: N>=pulsePerSqnc 
                %In case N~=pulsePerSqnc, fSqnc needs to be changed to fit N pulses.
                if N ~= pulsePerSqnc
                    fgClkSamplesPerSqnc = fgClkSamplesPerSqnc + (N-pulsePerSqnc)*fgClkSamplesPerPulse;
                    fSqnc               = fgClk / fgClkSamplesPerSqnc;
                    pulsePerSqnc        = N;
                    fprintf("NOTICE: Multiplexing matrix size is different than num of pos\n")
                end
            end
            
            % Create the signal
            samplesVec = 0: 1 : fgClkSamplesPerPulse-1;
            sinSamples = sin(2*pi * (fSin/fgClk) * samplesVec);
            
            if ~this.uVars.useHadamard
                sigVec = zeros(1, pulsePerSqnc);
                sigVec(1) = 1;
            else
                sigVec = sMatrix(:, 1)';
            end
            
            sigData = repmat(sigVec, fgClkSamplesPerPulse, 1);
            idxs = find(sigVec);
            sigData(:, idxs) = repmat(sinSamples', 1, length(idxs));
            sigData = sigData(:);

            %--------------------
            % Collect parameters
            %--------------------
            if this.uVars.useHadamard
                this.hadamard.N       = N;
                this.hadamard.sMat    = sMatrix;
%             this.hadamard.sMatInv = sMatInv;
            clear ('N', 'sMatrix', 'sMatInv');
            end
            
            clear('cycPerPulse', 'fSin', 'fgClk')

            vars = who();
            for i = 1:length(vars)
                if strcmp(vars{i}, 'this'); continue; end
                this.usSignal.(vars{i}) = eval(vars{i});
            end    
        end

        function initSamplesParameters(this)
            %--------------------------------------
            % Init Parameters (For Organization)
            %--------------------------------------
            this.samples.samplesPerSin          = [];
            this.samples.samplesPerPulse        = [];
            this.samples.samplesPerSqnc         = [];
            this.samples.samplesPerFrame        = [];
            this.samples.samplesPerSignal       = [];
            this.samples.samplesPerData         = [];
            this.samples.samplesPerMeas         = [];
            this.samples.samplesPerDataAllCh    = [];
            this.samples.samplesPerPosPerSqnc   = [];
            this.samples.samplesPerPosPerFrame  = [];
            this.samples.samplesPerPosPerSignal = [];
            
            this.samples.pulsePerSqnc           = [];
            this.samples.pulsePerFrame          = [];
            
            this.samples.sqncPerFrame           = [];
            this.samples.sqncPerSignal          = [];
            this.samples.sqncPerData            = [];
            this.samples.sqncPerDataAllCh       = [];
            
            this.samples.framesPerSignal        = [];
            
            this.samples.numOfPos               = [];
            this.samples.numOfPosUS             = [];
            
            this.samples.prePhantomSamples      = [];
            this.samples.inPhantomPropSamples   = [];
            this.samples.postSignalSamples      = [];
            this.samples.preSignalSamples       = [];
            
            %Raw
            this.samples.samplesPerFrameRaw        = [];
            this.samples.samplesPerSignalRaw       = [];
            this.samples.samplesPerPosPerFrameRaw  = [];
            this.samples.samplesPerPosPerSignalRaw = [];
            this.samples.samplesPerDataRaw         = [];
            this.samples.pulsePerFrameRaw          = [];
            this.samples.sqncPerFrameRaw           = [];
            this.samples.sqncPerSignalRaw          = [];
            this.samples.sqncPerSignalAllChRaw     = [];
            this.samples.sqncPerDataRaw            = [];
            this.samples.framesPerSignalRaw        = [];
            this.samples.framesPerSignalAllChRaw   = [];
            this.samples.numOfPosRaw               = [];

            % Hadamard
            this.samples.samplesPerSqncHad         = [];
            this.samples.samplesPerFrameHad        = [];
            this.samples.samplesPerSignalHad       = [];
            this.samples.samplesPerPosPerFrameHad  = [];
            this.samples.samplesPerPosPerSignalHad = [];
            
            this.samples.pulsePerFrameHad   = [];
            this.samples.sqncPerFrameHad    = [];
            this.samples.sqncPerSignalHad   = [];
            this.samples.framesPerSignalHad = [];
            
            this.samples.samplesPerDataInHad      = [];
            this.samples.samplesPerDataAllChInHad = [];
            this.samples.sqncPerDataInHad         = [];
            this.samples.sqncPerDataAllChInHad    = [];
            
            this.samples.sqncPerDataOutHad          = [];
            this.samples.sqncPerDataAllChOutHad     = [];
            this.samples.samplesPerDataOutHad       = [];
            this.samples.samplesPerDataAllChOutHad  = [];
            this.samples.samplesPerDataOutPreCutHad = [];

            % Continues Speckle
            this.samples.sqncPerDataInCS           = [];
            this.samples.sqncPerDatalAllChInCS     = [];
            this.samples.sqncPerFrameCS            = [];
            this.samples.samplesPerSqncCS          = [];
            this.samples.framesPerSignalCS         = [];
            this.samples.sqncPerSignalCS           = [];
            this.samples.samplesPerSignalCS        = [];
            this.samples.samplesPerPosPerSignalCS  = [];
            this.samples.samplesPerDataOutPreCutCS = [];
            this.samples.sqncPerDataOutCS          = [];
            this.samples.sqncPerDataAllChOutCS     = [];
            this.samples.samplesPerDataOutCS       = [];
            this.samples.samplesPerDataAllChOutCS  = [];
        end
        
        function calcSamples(this)
            % Terminology:
            % Sin   - single sinus
            % Pulse - few sinuses
            % sqnc  - sequence of pulses
            % frame - # of sqnces under SDT
            % Pos   - # of pulses of a region that fits in frame
            % sig   - # of frames in time to sample
            % data  - sig  + sampls to allow contSpeckle & contHadamard (excessive sqnc/frames)
            % meas  - data + propagation samples (pre phantom + in phantom)
            
            this.initSamplesParameters();
            
            %--------------------------------------
            % Collect Variables
            %--------------------------------------
            c               = this.geometry.c;
            distFromPhantom = this.geometry.distFromPhantom;
            channels        = this.general.channelsToAnalyze;
            
%             preTriggerSamples = this.uVars.preTriggerSamples;

            timeToSample = this.general.timeToSample;
            fs           = this.sClk.fs;
            fSin         = this.uVars.fSin;
            cycPerPulse  = this.uVars.cycPerPulse;
            fSqnc        = this.usSignal.fSqnc;
            pulsePerSqnc = this.usSignal.pulsePerSqnc;
            frameTime    = this.uVars.frameTime;
            if ~this.uVars.useFrame
                this.uVars.frameTime = this.general.timeToSample;
            end
            
            artfctIdxVec = this.uVars.artfctIdxVec;
            
            %--------------------------------------
            % ---- Raw Signal Calculation
            %--------------------------------------
            sqncPerFrameRaw  = ceil(frameTime * fSqnc);
            
            sqncPerSignalRaw = ceil(timeToSample * fSqnc);                             %Make sure that in timeToSample there is a whole number of sqnces
            sqncPerSignalRaw = ceil(sqncPerSignalRaw / sqncPerFrameRaw) * sqncPerFrameRaw; %Make sure that in timeToSample there is a whole number of frames
            
            sqncPerSignalAllChRaw   = sqncPerSignalRaw * channels;
            
            framesPerSignalRaw      = sqncPerSignalRaw   / sqncPerFrameRaw;
            framesPerSignalAllChRaw = framesPerSignalRaw * channels;
            
            samplesPerSin       = fs/fSin;
            samplesPerPulse     = samplesPerSin   * cycPerPulse;
            samplesPerSqnc      = samplesPerPulse * pulsePerSqnc;
            samplesPerFrameRaw  = samplesPerSqnc  * sqncPerFrameRaw;
            samplesPerSignalRaw = samplesPerFrameRaw * framesPerSignalRaw;
            
            samplesPerPosPerSqnc       = samplesPerPulse;
            samplesPerPosPerFrameRaw   = samplesPerPosPerSqnc * sqncPerFrameRaw;
            samplesPerPosPerSignalRaw  = samplesPerPosPerFrameRaw * framesPerSignalRaw; %not very useful info
            
            pulsePerFrameRaw = pulsePerSqnc * sqncPerFrameRaw;
            
            %TODO: reduce this methodology
            numOfPosRaw   = pulsePerSqnc;
            numOfPosUS    = pulsePerSqnc;
            numOfPos      = numOfPosRaw;
            numOfPosAlgo  = numOfPosRaw;

            %--------------------------------------
            % ---- Meas Signal & Data Size Calculation
            %--------------------------------------
            prePhantomSamples     = ceil( (distFromPhantom / c) * fs ); % the propagation of the first ultrasound pulse from transducer to medium surface
            inPhantomPropSamples  = 0;
            postSignalSamples     = 0;
            
            if this.uVars.useHadamard
                inPhantomPropSamples  = 0;
%                 inPhantomPropSamples  = samplesPerSqnc - samplesPerPulse;
%                 postSignalSamples     = samplesPerSqnc;
            end

            samplesPerDataRaw = samplesPerSignalRaw  + postSignalSamples;
            preSignalSamples  = prePhantomSamples    + inPhantomPropSamples;
            samplesPerMeas    = preSignalSamples     + samplesPerDataRaw;
            
            sqncPerDataRaw    = samplesPerDataRaw / samplesPerSqnc;
            
            %------------------------------
            % Collect Raw Parameters
            %------------------------------
            sqncPerFrame           = sqncPerFrameRaw;
            sqncPerSignal          = sqncPerSignalRaw;
            
            framesPerSignal        = framesPerSignalRaw;
            
            samplesPerFrame        = samplesPerFrameRaw;
            samplesPerSignal       = samplesPerSignalRaw;
            samplesPerPosPerFrame  = samplesPerPosPerFrameRaw;
            samplesPerPosPerSignal = samplesPerPosPerSignalRaw;

            pulsePerFrame          = pulsePerFrameRaw;
            
            %----------------------------------------------------
            % Calculate Raw Data Size
            samplesPerData      = samplesPerDataRaw;
            samplesPerDataAllCh = samplesPerData*channels;
            sqncPerData         = sqncPerDataRaw;
            sqncPerDataAllCh    = sqncPerData*channels;
 
            %--------------------------------------
            % ---- Hadamard Calculation
            %--------------------------------------
            % Creating the sMatrix is hapeening inside calcUsSignalDim when
            % constructing the US transmission signal
            
            if this.uVars.useHadamard
                sMatrix = this.hadamard.sMat;
                sMatInv = inv(sMatrix);
                N       = this.hadamard.N;
                
                %-------------------------------------
                % Collect Input Data Size
                samplesPerDataInHad      = samplesPerData;
                samplesPerDataAllChInHad = samplesPerDataInHad * channels;
                sqncPerDataInHad         = sqncPerData;
                sqncPerDataAllChInHad    = sqncPerDataInHad * channels;
                
                %--------------------------------------
                % ---- Naive Hadamard
                %--------------------------------------
                sVecInv = sMatInv(1,:);

                sVecInvSqnc      = zeros(samplesPerPulse, N);
                sVecInvSqnc(1,:) = sVecInv;
                sVecInvSqnc      = sVecInvSqnc(:);

                sMatInvSqnc = zeros(samplesPerSqnc);
                for i=1:samplesPerSqnc
                    sMatInvSqnc(i,:) = circshift(sVecInvSqnc, (i-1));
                end

%                 for i = 1:N
%                     sMatInvSingleSqnc( (i-1)*samplesPerPulse+1:i*samplesPerPulse, ...
%                                         1:samplesPerPulse:end)         = ...
%                                         repmat(sMatInv(i,:),  samplesPerPulse, 1);
%                 end
%                 
%                 for i = 1:samplesPerPulse
%                     sMatInvSingleSqnc(i:samplesPerPulse:end, :) = circshift(sMatInvSingleSqnc(i:samplesPerPulse:end, :), (i-1), 2);
%                 end
                
%                 sMatInvSingleSqnc = flip(sMatInvSingleSqnc,1);
                
                
                %----------------------------------------------------------
                % Collect Naive-Hadamard Calculation to Hadamard Parameters
%                 sMatInvSingleSqnc         = sMatInvSingleSqnc;
                
                sqncPerFrameHad           = sqncPerFrame;
                framesPerSignalHad        = framesPerSignal;
                sqncPerSignalHad          = sqncPerSignal;

                samplesPerSqncHad         = samplesPerSqnc;
                samplesPerFrameHad        = samplesPerFrame;
                samplesPerSignalHad       = samplesPerSignal;

                samplesPerPosPerFrameHad  = samplesPerPosPerFrame;
                samplesPerPosPerSignalHad = samplesPerPosPerSignal;

                pulsePerFrameHad          = pulsePerFrame;

                %----------------------------------------------------
                % Calculate Naive-Hadamard Output Data Size
                sqncPerDataOutHad          = sqncPerDataInHad;
                sqncPerDataAllChOutHad     = sqncPerDataOutHad    * channels;
                samplesPerDataOutHad       = sqncPerDataOutHad    * samplesPerSqnc;
                samplesPerDataAllChOutHad  = samplesPerDataOutHad * channels;
                samplesPerDataOutPreCutHad = sqncPerDataInHad     * samplesPerSqnc;

%                 %---------------------------------
%                 % Continues Hadamard
%                 %---------------------------------
%                 if this.general.contHadamard
%                     sMatInvSqncCont = zeros(samplesPerSqnc*samplesPerSqnc, 2*samplesPerSqnc);
%                     for i=1:samplesPerSqnc
%                         tmp = circshift(sMatInvSingleSqnc, -(i-1), 2);
%                         sMatInvSqncCont( (i-1)*samplesPerSqnc+1 : (i*samplesPerSqnc), i:i+samplesPerSqnc-1) = tmp;
%                     end
% 
%                     sMatInvSqnc = sMatInvSqncCont;
%                     
%                     %----------------------------------------------------------
%                     % Collect Naive-Hadamard Calculation to Hadamard Parameters
%                     sqncPerFrameHad          = sqncPerFrameRaw * samplesPerSqnc;
%                     framesPerSignalHad        = framesPerSignal;
%                     sqncPerSignalHad          = framesPerSignalHad * sqncPerFrameHad;
% 
%                     samplesPerFrameHad        = sqncPerFrameHad   * samplesPerSqnc;
%                     samplesPerSignalHad       = samplesPerFrameHad * framesPerSignalHad;
% 
%                     samplesPerPosPerFrameHad  = samplesPerPulse * sqncPerFrameHad;
%                     samplesPerPosPerSignalHad = samplesPerPulse * samplesPerSignalHad;
% 
%                     pulsePerFrameHad          = pulsePerSqnc * sqncPerFrameHad;
%                     
%                     %----------------------------------------------------
%                     % Calculate Data Size
%                     sqncPerDataOutHad          = (sqncPerDataInHad-1) * samplesPerSqnc;
%                     sqncPerDataAllChOutHad     = sqncPerDataOutHad    * channels;
%                     samplesPerDataOutHad       = sqncPerDataOutHad    * samplesPerSqnc;
%                     samplesPerDataAllChOutHad  = samplesPerDataOutHad * channels;
%                     samplesPerDataOutPreCutHad = sqncPerDataInHad * samplesPerSqnc;
%                 end
                
                %---------------------------------------------
                % Collect Had Calculation To Global Parameters
                sqncPerFrame           = sqncPerFrameHad;
                sqncPerSignal          = sqncPerSignalHad;

                framesPerSignal        = framesPerSignalHad;

                samplesPerFrame        = samplesPerFrameHad;
                samplesPerSignal       = samplesPerSignalHad;
                samplesPerPosPerFrame  = samplesPerPosPerFrameHad;
                samplesPerPosPerSignal = samplesPerPosPerSignalHad;

                pulsePerFrame          = pulsePerFrameHad;
                
                samplesPerData = samplesPerDataOutHad;
                sqncPerData    = sqncPerDataOutHad;
            end
            
%             %--------------------------------------
%             % ---- Continues Speckle Calculation
%             %--------------------------------------
%             if this.general.contSpeckleAnalysis
%                 sqncPerDataInCS        = sqncPerData;
%                 sqncPerDatalAllChInCS  = channels * sqncPerData;
%                 
%                 sqncPerFrameCS           = sqncPerFrame;
%                 samplesPerSqncCS         = samplesPerSqnc;
%                 framesPerSignalCS        = framesPerSignal   * sqncPerFrameCS;
%                 sqncPerSignalCS          = framesPerSignalCS * sqncPerFrameCS;
%                 samplesPerSignalCS       = samplesPerFrame   * framesPerSignalCS;
%                 samplesPerPosPerSignalCS = samplesPerPulse * sqncPerFrame * framesPerSignalCS;
% 
%                 %----------------------------------------------------
%                 % Collect Cont Parameters
%                 sqncPerSignal          = sqncPerSignalCS;
%                 framesPerSignal        = framesPerSignalCS;
%                 samplesPerSignal       = samplesPerSignalCS;
%                 samplesPerPosPerSignal = samplesPerPosPerSignalCS;
%                 
%                 %----------------------------------------------------
%                 % calculate Data Size
%                 samplesPerDataOutPreCutCS = sqncPerDataInCS * sqncPerFrameCS * samplesPerSqncCS;
%                 sqncPerDataOutCS          = (sqncPerDataInCS - sqncPerFrameCS) * sqncPerFrameCS;
%                 sqncPerDataAllChOutCS     = sqncPerDataOutCS     * channels;
%                 samplesPerDataOutCS       = sqncPerDataOutCS     * samplesPerSqncCS;
%                 samplesPerDataAllChOutCS  = samplesPerDataOutCS  * channels;
%                 
%                 samplesPerData = samplesPerDataOutCS;
%                 sqncPerData    = sqncPerDataOutCS;
%             end
              
            %--------------------------------------
            % ---- Cutting Artifacts & MuEff
            %--------------------------------------
            if this.general.cutArtfct
                numOfPos = numOfPosAlgo - length(artfctIdxVec);
            end

            %--------------------
            % Collect All parameters
            %--------------------
            if  this.general.useHadamard
                this.hadamard.sMatInv = sMatInv;
                this.hadamard.sMatInvSqnc       = sMatInvSqnc;
            end
            
            clear ('c', 'channels', 'fSin', 'fSqnc', 'fs', 'sinPerPulse', ... 
                   'frameTime', 'timeToSample', 'preTriggerSamples', ...
                   'distFromPhantom', 'sMatInvSqnc', 'sMatInvSqnc', 'sMatInvSqncCont');
               
            vars = who();
            for i = 1:length(vars)
                if strcmp(vars{i}, 'this'); continue; end
                this.samples.(vars{i}) = eval(vars{i});
            end
        end

        function calcFreq(this)
            
            envDC  =  this.uVars.envDC;
            envUS  =  this.uVars.envUS;
            
            N = this.samples.samplesPerPosPerFrame;
            
            fUS = this.uVars.fSin;
            fs  = this.sClk.fs;
            
            df               = fs / N;
            samplesPerUsHar  = floor(fUS / df);
            shiftFactor      = floor(N/2);
            
            fBarRaw = (0 : 1 : N-1) * (fs/N);
            fBar    = (fs/N) *  ( (-N/2) : 1 : (N/2)-1 );

            fUsIdxPos = floor(fUS /df) +1;
            fUsIdxNeg = floor( (fs - fUS) /df) + 1;
            fDcIdx     = 1;
            
            fUsIdxPosShift = fUsIdxPos + shiftFactor;
            fUsIdxNegShift = fUsIdxNeg - shiftFactor;
            fDcIdxShift    = fDcIdx    + shiftFactor;
            
            fUSIdx = fUsIdxPosShift;
            
            harFreqs = fs/2 : -fUS : (-fs/2);
            harIdxs  = 1 : samplesPerUsHar : N+1;
            harIdxs(harIdxs == fDcIdxShift) = [];
            numOfHar = length(harIdxs);
            
            % Calc indeces of Harmonics
            envHarLen     = ceil(envUS / df);
            envVec        = -envHarLen : 1 : envHarLen;
            envMat        = repmat(envVec', 1, numOfHar);
            envHarFullLen = length(envVec);

            fUSIdxs = fUSIdx-envHarLen : fUSIdx+envHarLen;
            fUSIdxsNoise = fUSIdx-3*envHarLen : fUSIdx+3*envHarLen;
            fUSIdxsNoise (ismember(fUSIdxsNoise, fUSIdxs)) = [];

            harEnvIdxMat = repmat(harIdxs, envHarFullLen, 1);
            harEnvIdxs   = envMat + harEnvIdxMat;
            harEnvIdxs   = harEnvIdxs(:)';
            harEnvIdxs(harEnvIdxs <= 0 | harEnvIdxs > N) = [];
            
            if envUS == 0
                harEnvIdxs = [fUsIdxNegShift, fUsIdxPosShift];
            end
            % Calc indeces of DC
            envDCLen  = floor(envDC / df);
            envVec    = -envDCLen : 1 : envDCLen;
            dcEnvIdxs = fDcIdxShift + envVec;
            
            cutIdxs = [harEnvIdxs, dcEnvIdxs];
            fitIdxShift = 1 : 1 : N;
            fitIdxShift (cutIdxs ) = [];

            if this.uVars.useGPU
                fBarRaw        = gpuArray(fBarRaw);
                fBar           = gpuArray(fBar);
                fitIdxShift    = gpuArray(fitIdxShift);
            end
            
            %--------------------
            % Collect parameters
            %--------------------
            vars = who();
            for i = 1:length(vars)
                if strcmp(vars{i}, 'this'); continue; end
                this.freq.(vars{i}) = eval(vars{i});
            end
        end
        
        function calcTiming(this)
            fs    = this.sClk.fs;
            fSin  = this.uVars.fSin;
            fSqnc = this.usSignal.fSqnc;
            fgClk = this.sClk.fgClk;
            cycPerPulse = this.uVars.cycPerPulse;
            
            samplesPerPulse       = this.samples.samplesPerPulse;
            samplesPerSqnc        = this.samples.samplesPerSqnc;
            samplesPerFrame       = this.samples.samplesPerFrame;
            samplesPerPosPerFrame = this.samples.samplesPerPosPerFrame;
            samplesPerSignal      = this.samples.samplesPerSignalRaw;
            samplesPerMeas        = this.samples.samplesPerMeas;
            
            fgClkSamplesPerFsSig = this.sClk.fgClkSamplesPerFsSig;
            fgClkSamplesPerSqnc  = this.usSignal.fgClkSamplesPerSqnc;
            
            % Timing Vars
            dts      = 1/fs;     %[s], duration of a single sample
            dtFgClk  = 1/fgClk;  %[s], duration of a single function generator clk
            
            sinT    = 1/fSin;              %[s], duration of single period of the sin
            pulseT  = cycPerPulse * sinT;  %[s], duration of single pulse
            sqncT   = 1/fSqnc;             %[s], duration of single sqnce
            frameT  = samplesPerFrame /fs; %[s], duration of single Frame
            signalT = samplesPerSignal /fs; %#ok<*NASGU>

            tVecSampleClk = ( 0 : 1 : (fgClkSamplesPerFsSig - 1) ) * dtFgClk;
            tVecUS        = ( 0 : 1 : (fgClkSamplesPerSqnc - 1) )  * dtFgClk;
            
            tVecPulse  = ( 0 : 1 : (samplesPerPulse   - 1) ) * dts;
            tVecSqnc   = ( 0 : 1 : (samplesPerSqnc   - 1) ) * dts;
            tVecFrame  = ( 0 : 1 : (samplesPerFrame   - 1) ) * dts;
            tVecSig    = ( 0 : 1 : (samplesPerSignal - 1) ) * dts;

            tVecPosPerFrame  = ( 0 : 1 : (samplesPerPosPerFrame - 1) ) * dts;
            if ~this.general.splitMeas 
                tVecMeas   = ( 0 : 1 : (samplesPerMeas    - 1) ) * dts;
            else
%                 tVecPulse        = [];
%                 tVecSqnc         = [];
%                 tVecFrame        = [];
%                 tVecSig          = [];
%                 tVecPosPerFrame  = [];
                tVecMeas         = [];
            end
            
            %--------------------
            % Collect parameters
            %--------------------
            clear ('fs', 'fSin', 'fSqnc', 'fgClk', 'cycPerPulse', ... 
                   'samplesPerPulse', 'samplesPerSqnc', 'samplesPerFrame', ...
                   'samplesPerPosPerFrame', 'samplesPerSignal', 'samplesPerMeas', ... 
                   'fgClkSamplesPerFsSig', 'fgClkSamplesPerSqnc');
            
            vars = who();
            for i = 1:length(vars)
                if strcmp(vars{i}, 'this'); continue; end
                this.timing.(vars{i}) = eval(vars{i});
            end
            
        end

        function calcLength(this)
            c      = this.geometry.c;
            
            dts    = this.timing.dts; 
            sinT   = this.timing.sinT;
            pulseT = this.timing.pulseT;
            frameT = this.timing.frameT;
            
%             numOfPosHR = this.samples.numOfPosHR;
            numOfPosAlgo = this.samples.numOfPosAlgo;
            numOfPos     = this.samples.numOfPos;
            
            % Length Vars
            sinLen     = c * sinT;
            pulseLen   = c * pulseT;
            frameLen   = c * frameT;
            
            dDepth = c * pulseT;
            depthLen = frameLen;
            
            depthVecAlgo    = ( 0 : 1 : (numOfPosAlgo - 1) ) * dDepth;
            depthIdxAlgo    =  1 : 1 : numOfPosAlgo;
            depthIdxLenAlgo = numOfPosAlgo;
            
            depthVec    = ( 0 : 1 : (numOfPos - 1) ) * dDepth;
            depthIdxLen = numOfPos;

            depthNorm   = depthVec;
            depthCntr   = depthNorm - mean(depthNorm);
            depthZero   = depthNorm - min(depthNorm);
            depthIdx    = 1 : 1 : depthIdxLen;
            
            depthVecNormAlgo = 1 : 1 : numOfPosAlgo;
            %--------------------
            % Collect parameters
            %--------------------
            clear ('c', 'dts', 'tSin', 'tPulse', 'tFrame', ... 
                   'numOfPosHR', 'numOfPosUS');
            vars = who();
            for i = 1:length(vars)
                if strcmp(vars{i}, 'this'); continue; end
                this.len.(vars{i}) = eval(vars{i});
            end
        end

        %% Algorithm Functions
        function reconAlgo(this)
            % measSamples(input)  - [ch x samplesPerMeas]
            % res(output) - struct:
            %               measSamples - [ch x samplesPerMeas]
            %               signal      - [ch x samplesPerSignal]
            %               reshaped    - [ch x frames x samplesPerPos x numOfPos]
            %               fft      - [ch x samplesPerPos x numOfPos]
            %               phiCh       - [ch x 1 x numOfPos] 
            %               phi         - [1  x 1 x numOfPos]
            
            this.timeTable.FullAnalysis = tic;

            % Export Raw Data
            this.exportData('meas');
            
            % Extraxt Net signal 
            this.timeTable.overallNetSignal = tic;
            this.extractSignal();
            this.timeTable.overallNetSignal = toc(this.timeTable.overallNetSignal);
            
            % DeMultiplex
            this.timeTable.overallDemultiplex = tic;
            if this.general.useHadamard
                this.demultiplexSignal();
            end
            this.timeTable.overallDemultiplex = toc(this.timeTable.overallDemultiplex);
            
            this.timeTable.overallACCoupling = tic;
            this.acCoupling();
            this.timeTable.overallACCoupling = toc(this.timeTable.overallACCoupling);
            
            % Reshape Data
            this.timeTable.overallReshape = tic;
            this.reshapeSignal();
            this.timeTable.overallReshape = toc(this.timeTable.overallReshape);

            % Signal Processing
            this.timeTable.overallSignalProcessing = tic;
            this.fourierTransform();
            this.timeTable.overallSignalProcessing = toc(this.timeTable.overallSignalProcessing);

            if ~this.general.splitMeas
                this.extractPhiFromFFT();
            else
                this.gatherCurrentSplit();
            end
                        
            this.timeTable.FullAnalysis = toc(this.timeTable.FullAnalysis);
        end
        
        function res = reconExtFFTData(this, extData)
            this.result = [];
            models = this.getResultsModels();

            for i=1:length(extData)
                this.curRes = models.base;
%                 if isfield(extData(i), 'splitRes')
%                     this.curRes = models.baseSplit;
%                 else
%                     this.curRes = models.base;
%                 end
                this.curRes.frameAvgPowerFFT = extData(i).frameAvgPowerFFT;
                
                this.extractPhiFromFFT();
                if i==1
                    this.result = this.curRes;
                else
                    this.result(i) = this.curRes;
                end
            end
            
            res = this.result;
        end
        
        function extractSignal(this)
            % measSamples(input)  - [ch x samplesPerAcq]
            % measSamples(output) - [ch x samplesPerSignal]
            
            % Collect Variables
            preSignalSamples    = this.samples.preSignalSamples;
            samplesPerMeas      = this.samples.samplesPerMeas;
            
            % Extract netSignal samples
            this.timeTable.extractNetSignal = tic;
            this.data = this.data(:, ( (preSignalSamples + 1) : samplesPerMeas ) );
            this.timeTable.extractNetSignal = toc(this.timeTable.extractNetSignal);
            
            this.exportData ('signal')
        end
        
        function demultiplexSignal(this)
            % NetSignal(output) - [ch x samplesPerSignal]
            % deMulNetSignal - [ch x samplesPerSignal]
            
            % Collect Variables
            channels                   = this.general.channelsToAnalyze;
            samplesPerSqncHad          = this.samples.samplesPerSqncHad;
            sqncPerDataAllChInHad      = this.samples.sqncPerDataAllChInHad;
            sMatInvSqnc                = this.hadamard.sMatInvSqnc;
            samplesPerDataOutHad       = this.samples.samplesPerDataOutHad;
            samplesPerDataOutPreCutHad = this.samples.samplesPerDataOutPreCutHad;

            % Separate each sequence
            this.timeTable.deMulReshape1 = tic;
            this.data = reshape(this.data', samplesPerSqncHad, sqncPerDataAllChInHad);
            this.timeTable.deMulReshape1 = toc(this.timeTable.deMulReshape1);

            sMatInvSqnc = gpuArray(sMatInvSqnc);

            this.timeTable.deMultiplex = tic;
            this.data   = sMatInvSqnc * this.data;
            this.timeTable.deMultiplex = toc(this.timeTable.deMultiplex);

            this.timeTable.deMulReshape2 = tic;
            this.data = reshape(this.data, samplesPerDataOutPreCutHad, channels)';
            this.timeTable.deMulReshape2 = toc(this.timeTable.deMulReshape2);

            this.data(:, samplesPerDataOutHad+1:end) = [];
            
            this.exportData ('deMul')
        end
        
        function acCoupling(this)
            if this.general.acCoupling
                this.data = this.data - mean(this.data, 2);
            end
        end
        
        function reshapeSignal(this)
            % signal(input)  - [ch x samplesPerSignal]
            % signal(output) - [frame x ch x samplesPerPos x numOfPos]                                

            % Collect Variables
            channels              = this.general.channelsToAnalyze;
            numOfPos              = this.samples.numOfPosAlgo;
            samplesPerPulse       = this.samples.samplesPerPulse;
            pulsePerSqnc          = this.samples.pulsePerSqnc;
            samplesPerSqnc        = this.samples.samplesPerSqnc;
            framesPerSignal       = this.samples.framesPerSignal;
            
            samplesPerFrame       = this.samples.samplesPerFrame;
            sqncPerFrame          = this.samples.sqncPerFrame;
            pulsePerFrame         = this.samples.pulsePerFrame;
            samplesPerPosPerFrame = this.samples.samplesPerPosPerFrame;
         
            
            % Separate signal to Frames
            % [ch x samples x frame]
            this.timeTable.reshape1 = tic;
            this.data = reshape(this.data, channels,...
                                           samplesPerFrame, ...
                                           framesPerSignal);
            this.timeTable.reshape1 = toc(this.timeTable.reshape1);          

            this.data = this.data - mean(this.data, 2);

            % Bring frame-dim to be first (averaging convinient)
            this.timeTable.permute1 = tic;
            this.data = permute(this.data, [1,3,2]);
            this.timeTable.permute1 = toc(this.timeTable.permute1);

            % Break all the pulses in each frame to pulses from same sqnc
            % [frames x ch x samples x pulsePerSqnc x sqncPerFrame]
            this.timeTable.reshape2 = tic;
            this.data = reshape(this.data, channels,... 
                                           framesPerSignal,...
                                           samplesPerPulse,...
                                           pulsePerSqnc,...
                                           sqncPerFrame);
            this.timeTable.reshape2 = toc(this.timeTable.reshape2); 

            % Switch between sqnc and pulse dims
            % [frames x ch x samples x sqncPerFrame x pulsePerSqnc]
            this.timeTable.permute2 = tic;
            this.data = permute(this.data, [1,2,3,5,4]);
            this.timeTable.permute2 = toc(this.timeTable.permute2); 

            % Merge all samples of a single Poistion from all the sqnces
            % within the frame
            % [frames x ch x samples x numOfPos]
            this.timeTable.reshape3 = tic;
            this.data = reshape(this.data, channels,...
                                           framesPerSignal,...
                                           samplesPerPosPerFrame,...
                                           numOfPos);
            this.timeTable.reshape3 = toc(this.timeTable.reshape3);                       
            
            this.data = this.data - mean(this.data, 3);
            
            this.exportData ('reshaped');
        end
        
        function fourierTransform(this)
            % data(input)               - [ch x frames x samples x numOfPos]]
            % frameAvgPowerFFT(output)  - [ch x samplesPerPos x numOfPos]
            
            samplesPerPosPerFrame = this.samples.samplesPerPosPerFrame;

            % Perform Fourier Transform
            this.timeTable.FFT = tic;
            this.data = (2./samplesPerPosPerFrame)*fft(this.data, [], 3);
            this.timeTable.FFT = toc(this.timeTable.FFT);

            this.exportData('FFT')
            this.exportData('usComp')

            % Absolute and square
            this.timeTable.absAndPower = tic;
            this.data                  = abs(this.data).^2;
            this.timeTable.absAndPower = toc(this.timeTable.absAndPower);
            
            % Average over all frames
            % NOTE: shift is here to perform on a lower dimension data.
            this.timeTable.meanFrameFFT = tic;
            this.curRes.frameAvgPowerFFT = permute(mean(this.data, 2), [1,3,4,2]);
            this.timeTable.meanFrameFFT = toc(this.timeTable.meanFrameFFT);
            
            this.timeTable.Shift = tic;
            this.curRes.frameAvgPowerFFT = fftshift(this.curRes.frameAvgPowerFFT, 2);
            this.timeTable.Shift = toc(this.timeTable.Shift);

        end
        
        %% Post Processing:
        function extractRawPhi(this)
            %-----------------------------------------------------
            % Raw Phi Calculation (no fit - averaging over all ch)
            %-----------------------------------------------------
            channels = this.general.channelsToAnalyze;
            fUsIdx = this.freq.fUSIdx;

            % Average over all channels and then SQRT(!)
            this.timeTable.meanChFFTRaw = tic;
            this.curRes.rawChAvgFFT     = sqrt(permute(mean(this.curRes.frameAvgPowerFFTCut, 1), [2,3,1]));
            this.timeTable.meanChFFTRaw = toc(this.timeTable.meanChFFTRaw);
            
            this.timeTable.rawPhiExtract = tic;
            this.curRes.rawPhi = this.curRes.rawChAvgFFT(fUsIdx, :);
            this.timeTable.rawPhiExtract = toc(this.timeTable.rawPhiExtract);

            this.timeTable.AnalyseRawPhi = tic;
            this.curRes.analysis.rawPhi = this.analyseRecon(this.curRes.rawPhi);
            this.timeTable.AnalyseRawPhi = toc(this.timeTable.AnalyseRawPhi);
        end

        function extractPhiFromFFT(this)
            % Input:
            % frameAvgPowerFFT - [ch x samplesPerPos x numOfPos]
            % Outputs:
            % rawChAvgFFT    - [samplesPerPos x numOfPos]
            % rawChAvgFFTStd - [samplesPerPos x numOfPos]
            % rawPhi            - [1 x numOfPos]       
            % posAvgPowerFFT       - [ch x samplesPerPos]
            % fitModelMat       - [ch x samplesPerPos]
            % fittedQAvgFFT     - [ch x samplesPerPos x numOfPos]
            % fittedChAvgFFT    - [samplesPerPos x numOfPos]
            % fittedFFT         - [samplesPerPos x numOfPos]
            % fittedFFTNorm     - [samplesPerPos x numOfPos]
            % phi               - [1 x numOfPos]      
            
            % Collect Variables
            channels = this.general.channelsToAnalyze;
            fUsIdx = this.freq.fUSIdx;
            fUSIdxs = this.freq.fUSIdxs;

            %Cut Artifacts
            this.cutArtifactIndexes();
            
            % Extract Raw Phi - Naive Signal Processing
            this.extractRawPhi();
            
            %-----------------------------------------------------
            % Phi Calculation (fit each channel and then average)
            %-----------------------------------------------------
            
            %Averag over all positions to create avg channel response before fit
            this.timeTable.meanPos = tic;
            this.curRes.posAvgPowerFFT    = mean(this.curRes.frameAvgPowerFFTCut(:,:,this.curRes.analysis.rawPhi.bkgIdx), 3);
            this.timeTable.meanPos = toc(this.timeTable.meanPos);
            
            % Calculate fit model to FFT
            this.timeTable.calculateFit = tic;
            if this.general.calibrate
                this.curRes.fitModelMat = this.getCalibration();
            else
                this.curRes.fitModelMat = this.calcFittedFFT();
            end
            this.timeTable.calculateFit = toc(this.timeTable.calculateFit);
            
            % Apply Fit
            this.timeTable.applyFit = tic;
            this.curRes.fittedPowerFFT = this.curRes.frameAvgPowerFFTCut ./ this.curRes.fitModelMat;
%             this.curRes.fittedPowerFFT = this.curRes.frameAvgPowerFFTCut - this.curRes.fitModelMat;
            this.timeTable.applyFit    = toc(this.timeTable.applyFit);

            % Average over all channels
            this.timeTable.meanChFFT = tic;
            this.curRes.fittedChAvgPowerFFT = permute(mean(this.curRes.fittedPowerFFT, 1), [2,3,1]);
            this.timeTable.meanChFFT = toc(this.timeTable.meanChFFT);

            % Sqrt FFT
            this.timeTable.SqrtAbsFittedFFT = tic;
%             this.curRes.fittedFFT = sqrt(abs(this.curRes.fittedChAvgPowerFFT - 1));
            this.curRes.fittedFFT = abs(sqrt(abs(this.curRes.fittedChAvgPowerFFT)) - 1);
%             this.curRes.fittedFFT = sqrt(abs(this.curRes.fittedChAvgPowerFFT));
            this.timeTable.SqrtAbsFittedFFT = toc(this.timeTable.SqrtAbsFittedFFT);

            % Extract Ultrasound Component
            this.timeTable.phiExtract = tic;
            if false
                this.curRes.phiPreCut = this.curRes.fittedFFT(fUSIdxs, :);
            else
                this.curRes.phiPreCut = max(this.curRes.fittedFFT(fUsIdx, :), [], 1);
            end
            this.timeTable.phiExtract = toc(this.timeTable.phiExtract);
            
            this.curRes.phi = this.curRes.phiPreCut;
            
%             figure();
%             subplot(3,2,1)
%             stem(this.curRes.rawPhi)
% %             ylim([0, 2e-7])
%             title("Raw")
%             subplot(3,2,2)
%             stem(sqrt(abs(squeeze(this.curRes.frameAvgPowerFFT(1, 91,:)))))
%             title("Channel Before Fit")
%             subplot(3,2,3)
%             stem(sqrt(abs(squeeze(this.curRes.fittedPowerFFT(1, 91,:)))))
%             title("Channel After Fit")
%             subplot(3,2,4)
%             stem(sqrt(abs(squeeze(this.curRes.fittedChAvgPowerFFT(91,:)))))
%             title("Phi After Ch Averaging")
%             subplot(3,2,5)
%             stem(this.curRes.fittedFFT(91,:))
%             title("Phi After Sqrt")
%             subplot(3,2,6)
%             stem(this.curRes.phi)
%             title("Final Phi")
            
%             figure();
%             subplot(2,2,1)
%             plot(squeeze(this.curRes.frameAvgPowerFFTCut(1,:,103)))
%             subplot(2,2,2)
%             plot(squeeze(this.curRes.fitModelMat))
%             subplot(2,2,3)
%             plot(this.curRes.fittedChAvgPowerFFT(:, 103))

            % Mark Signal and Background indices
            this.timeTable.AnalysePhi = tic;
            this.curRes.analysis.phi = this.analyseRecon(this.curRes.phi);
            this.timeTable.AnalysePhi = toc(this.timeTable.AnalysePhi);
            
            %Selected Normalization:
            this.curRes.phiNorm = this.curRes.analysis.phi.normTypes.phiNorm2;
            this.curRes.phiLog  = this.curRes.analysis.phi.normTypes.phiLog2;
            
            % Calc MuEff
            this.calcMuEff();
            
            % Gather from GPU
            this.timeTable.gather = tic;
            
            this.curRes.frameAvgPowerFFT = gather(this.curRes.frameAvgPowerFFT);
            
            this.curRes.frameAvgPowerFFTCut = gather(this.curRes.frameAvgPowerFFTCut);
            
            this.curRes.rawChAvgFFT    = gather(this.curRes.rawChAvgFFT);
            this.curRes.rawPhi            = gather(this.curRes.rawPhi);
            
            this.curRes.posAvgPowerFFT       = gather(this.curRes.posAvgPowerFFT);
            this.curRes.fitModelMat          = gather(this.curRes.fitModelMat);
            this.curRes.fittedPowerFFT       = gather(this.curRes.fittedPowerFFT);
            this.curRes.fittedChAvgPowerFFT  = gather(this.curRes.fittedChAvgPowerFFT);
            this.curRes.fittedFFT            = gather(this.curRes.fittedFFT);

            this.curRes.phiPreCut      = gather(this.curRes.phiPreCut);
            this.curRes.phi            = gather(this.curRes.phi);
            this.curRes.phiLog         = gather(this.curRes.phiLog);
            this.curRes.phiNorm        = gather(this.curRes.phiNorm);
            
            this.timeTable.gather        = toc(this.timeTable.gather);
        end
        
        function cutArtifactIndexes(this)
            this.curRes.frameAvgPowerFFTCut = this.curRes.frameAvgPowerFFT;
            
            this.timeTable.CutArtifacts = tic;
            if this.general.cutArtfct
                this.curRes.frameAvgPowerFFTCut(:,:,this.samples.artfctIdxVec) = [];
            end
            this.timeTable.CutArtifacts = toc(this.timeTable.CutArtifacts);
        end
        
        function cal = getCalibration(this)
            cal = this.dataCal.raw;
        end
        
        function fitModelMat = calcFittedFFT(this)
            fitModelMat = ones(size(this.curRes.posAvgPowerFFT));
            if this.general.useGPU
                fitModelMat = gpuArray(fitModelMat);
            end

            if length(this.freq.fitIdxShift) <11
                return;
            end
            
            % Collect Variables
            channels  = this.general.channelsToAnalyze;
            numOfPos  = this.samples.numOfPos;
            fBar      = this.freq.fBar;
            idxs      = this.freq.fitIdxShift;

            % Calculate Fit
            for i = 1:channels
                fitVec = this.curRes.posAvgPowerFFT(i, idxs);
                fitModelMat(i,:) = interp1(fBar(idxs), fitVec, fBar, 'linear', 'extrap');
%                 fitModelMat(i,:) = interp1(fBar(idxs), fitVec, fBar, 'spline', 'extrap');
            end

%             fitModelMat = repmat(fitModelMat, 1, 1, numOfPos);
        end
        
        function analysis = analyseRecon(this, phi)
            phi = gather(phi);

            bkgIdxUser  = this.general.extNoiseIdx;
            peakIdxUser = this.general.extPeakIdx;

            fittedFFT = gather(this.curRes.fittedFFT);
            fUSIdx    = this.freq.fUSIdx;
            usEnvIdx  = this.freq.fUSIdxsNoise;

            %-------------------------
            % Separate Bakground and Signal:
            %-------------------------
            idxVec = 1:length(phi);

            minVal = min(phi);
            
            if this.general.extSNRDef
                peakIdx = peakIdxUser;
                maxVal  = phi(peakIdx);
                spanAbs =  maxVal - minVal;

                tmpStd      = std(phi(bkgIdxUser));
                tmpBkgLevel = mean(phi(bkgIdxUser));
                bkgIdxMask  = phi < (tmpBkgLevel+3*tmpStd);
            else
                [maxVal, peakIdx] = max(phi);
                spanAbs = maxVal - minVal;

                noiseLevel = 0.2;
                
                tmpStd      = std(phi( phi<(minVal+4*std(phi)) ) );
                tmpBkgLevel = mean(phi( phi<(minVal+4*std(phi)) ));
                bkgIdxMask  = phi < (minVal + spanAbs*noiseLevel);              
            end
            
            bkgIdx = idxVec(bkgIdxMask);
            bkg    = phi(bkgIdx);
            
            signalIdx = idxVec(~bkgIdxMask);
            signal    = phi(signalIdx);

            depthVecBkg    = this.len.depthVec(bkgIdx)*1e3;
            depthVecSignal = this.len.depthVec(signalIdx)*1e3;
            
            %-------------------------
            % Calculate SNR & Statisics:
            %-------------------------
            %Spatial SNR:
            stdNoise  = std(bkg);
            avgNoise  = mean(bkg);
            spanNoise = maxVal - avgNoise;
            
            SNRAbs = spanAbs/stdNoise;
            SNR = spanNoise/stdNoise;
            
            %Spectral SNR:
            if ~isempty(fittedFFT)
                spectSamples = fittedFFT(:, peakIdx);
    
                spectSignal    = spectSamples(fUSIdx);
                spectNoiseMean = mean(spectSamples(usEnvIdx));
                spectNoiseStd  = std(spectSamples(usEnvIdx));
                
                SNRSpectral = (spectSignal - spectNoiseMean) / spectNoiseStd;
            else
                SNRSpectral = [];
            end
            %-------------------------
            % Norm and Log:
            %-------------------------
            % ----- Phi:
            % Type 1 - Minimum:
            normTypes.phiNorm1 = (phi - minVal) ./ spanAbs;
            normTypes.phiLog1  = log(normTypes.phiNorm1);

            % Type 2 - Average Noise:
            normTypes.phiNorm2 = (phi - avgNoise) ./ spanNoise;
            normTypes.phiLog2  = log(abs(normTypes.phiNorm2));
            
            % Type3 - Maximum:
            normTypes.phiNorm3 = phi ./ maxVal;
            normTypes.phiLog3  = log(abs(normTypes.phiNorm3)); %#ok<STRNU> 
            
            % Processing:
            dX = this.len.dDepth*1e3;
            phiSmooth  = imfilter(phi, fspecial("average", [1,3]), 'circular', 'same');
            phiRoot    = abs(sqrt(phi));
            phiLaplace = sqrt(abs((gradient(gradient(phiRoot)*dX)*dX) ./ phiRoot));

            %--------------------
            % Collect parameters
            %--------------------
            clear ('phi', 'bkgIdxUser' ,'peakIdxUser', 'fittedFFT', 'fUSIdx', "usEnvIdx");
            clear ('idxVec', 'tmpStd', 'tmpBkgLevel', 'bkgIdxMask');
            clear ('spectSamples');

            vars = who();
            for i = 1:length(vars)
                if strcmp(vars{i}, 'this'); continue; end
                analysis.(vars{i}) = eval(vars{i});
            end
            
        end

        function calcMuEff(this)
            this.timeTable.CalcMuEff = tic;

            if this.general.calcMuEff
                phiCut = this.curRes.phiLog(this.general.muEffIdxs);
                phiRawCut = this.curRes.analysis.rawPhi.normTypes.phiLog2(this.general.muEffIdxs);
                xVecCut = this.len.depthVec(this.general.muEffIdxs)*1e3;
                dxInt   = (this.len.dDepth*1e3)/10;
                xVecInt = (xVecCut(1)):dxInt:(xVecCut(end));
                xVecIntIdx = (this.general.muEffIdxs(1) : 0.1 : this.general.muEffIdxs(end));
                phiCutInt = interp1(xVecCut, phiCut', xVecInt, 'pchip')';
                phiRawCutInt = interp1(xVecCut, phiRawCut', xVecInt, 'pchip')';

                switch this.general.muEffModel
                    case 'Uniform'
                       gradVals = gradient(phiCutInt /2);
                       muEffVal = mean(abs(gradVals(2:end-1))) /dxInt;
                       fitRes      = fit(xVecInt', cast(phiCutInt/2, 'double'), 'poly1');
                       muEffFitVal = fitRes.p1;

                       gradValsRaw = gradient(phiRawCutInt/2);
                       muEffRawVal = mean(abs(gradValsRaw(2:end-1)))  /dxInt;
                       fitRawRes      = fit(xVecInt', cast(phiRawCutInt/2, 'double'), 'poly1');
                       muEffRawFitVal = fitRawRes.p1;
                    case 'Point'

                    case 'TwoFibers'
                end

                vars = who();
                for i = 1:length(vars)
                    if strcmp(vars{i}, 'this'); continue; end
                    this.curRes.analysis.muEff.(vars{i}) = eval(vars{i});
                end
            end

            this.timeTable.CalcMuEff        = toc(this.timeTable.CalcMuEff);
        end
        
        %% Reconstruction types:
        function res = reconstruct(this) 
            this.resetAlgoRes();
            if this.general.analyzeSingleCh
                for i = 1:this.general.channelsToSample
                    this.data = this.measSamples (1, :);
                    this.measSamples (1, :) = [];
                    if this.general.useGPU
                        this.data = gpuArray(this.data);
                    end
                    this.reconAlgo();
                    this.result(i) = this.curRes;
                end
            else
                % In regular all-ch gpu based analysis measSamples will be
                % upload to GPU externally
                this.data = this.measSamples;
                this.measSamples = [];
                this.reconAlgo();
                this.result = this.curRes;
            end
            
            res = this.result;

            % All data in res is CPU
            this.data = [];
        end

        function res = reconSplittedData(this, dataInSplit)
            this.timeTable.OverallSplitAnalysis = tic;
            fprintf("ALGO: Analyzing splitted data...");
            this.resetSplitRes();
            splitNum = this.general.splitNum;
            
            splitFrameAvgPowerFFT = zeros(this.general.channelsToAnalyze, this.samples.samplesPerPosPerFrame);
            
            for j = 1:this.general.analysisReps

                this.timeTable.SplitsMeanFrameFFT  = tic;
                for i = 1:splitNum
                    splitFrameAvgPowerFFT = splitFrameAvgPowerFFT + dataInSplit(j,i).frameAvgPowerFFT;
                end
                splitFrameAvgPowerFFT = splitFrameAvgPowerFFT / splitNum;
                this.curRes.frameAvgPowerFFT = gpuArray(mean(splitFrameAvgPowerFFT, 4));
                this.timeTable.SplitsMeanFrameFFT = toc(this.timeTable.SplitsMeanFrameFFT);
                
                this.extractPhiFromFFT();
                this.copyStruct(this.curRes, 'result', j);
                this.result(j).splitRes = dataInSplit(j,:);
            end
            this.timeTable.OverallSplitAnalysis = toc(this.timeTable.OverallSplitAnalysis);
            res = this.result;
            fprintf("Done!\n");
        end

        function addCurrentSplit(this, curSplitData, idx)
            this.timeTable.AddSingleSplit = tic;
            for j = 1:this.general.analysisReps
                if idx==1
                    this.curSplit(j).frameAvgPowerFFT = curSplitData(j).frameAvgPowerFFT;
                else
                    this.curSplit(j).frameAvgPowerFFT = (1/idx)*(this.curSplit(j).frameAvgPowerFFT * (idx-1) + curSplitData(j).frameAvgPowerFFT);
                end
            end
            this.timeTable.AddSingleSplit = toc(this.timeTable.AddSingleSplit);
        end
        
        function resCur = reconCurSplit(this)
            for j = 1:this.general.analysisReps
                this.curRes.frameAvgPowerFFT   = gpuArray(this.curSplit(j).frameAvgPowerFFT);

                this.extractPhiFromFFT();
                this.copyStruct(this.curRes, 'result', j);
            end
            resCur = this.result;
        end
        
        function resNew = reconFromFFT(this, res)
            for j = 1:this.general.analysisReps
                this.curRes.frameAvgPowerFFT   = gpuArray(res(j).frameAvgPowerFFT);

                this.extractPhiFromFFT();
                this.copyStruct(this.curRes, 'result', j);
            end
            resNew = this.result;
        end

        function gatherCurrentSplit(this)
            this.curRes.frameAvgPowerFFT = gather(this.curRes.frameAvgPowerFFT);
        end

        %% Misc
        function bufferDataOut = createVirtualData(this, std)
            
            fus              = this.usSignal.fPulse;
            samplesPerPulse  = this.samples.samplesPerPulse;
            samplesPerSqnc   = this.samples.samplesPerSqnc;
            sqncPerSignal    = this.samples.sqncPerSignal;
            samplesPerSignal = this.samples.samplesPerSignal;
            samplesPerMeas   = this.samples.samplesPerMeas;
            pulsePerSqnc     = this.samples.pulsePerSqnc;
            
            tPulse   = this.timing.tVecSig(1:samplesPerPulse);
            sig      = sin(2*pi*fus*tPulse);
            sigSqnc = zeros(1, samplesPerSqnc);
            
            pointWithFluence = ceil(pulsePerSqnc/2);
            pointWithFluence = this.samples.numOfPos;
            factor = 150/pointWithFluence;
            factor = this.len.depthVec*1e3;
            muEff = 0.074;
%             expVec = (0:(pointWithFluence-1))*factor;
            expVec = factor;
            
            for i=1:pointWithFluence
                sigSqnc(samplesPerPulse*(i-1)+1:i*samplesPerPulse) = sig*(exp(-2*muEff*expVec(i)));
            end
            
            prePhantom  = this.samples.prePhantomSamples;
            postSamples = samplesPerMeas - (prePhantom + samplesPerSignal);
            
            preSig  = zeros(1, prePhantom);
            postSig = zeros(1, postSamples);

            tmpSig = repmat(sigSqnc, 1, sqncPerSignal);
            bufferDataOut = repmat([preSig, tmpSig, postSig], this.general.channels,1) + ...
                            std*rand(this.general.channels, samplesPerMeas);
            
            if this.general.useGPU
                bufferDataOut = gpuArray(bufferDataOut);
            end
        
            this.timeTable.moveData = tic;
        end
        
        function initTimeTable(this)
            this.timeTable = [];
            timeTableNames = ["FullAnalysis", ...
                              "overallNetSignal", "overallDemultiplex", "overallACCoupling", "overallReshape", "overallSignalProcessing",...
                              "extractNetSignal",...
                              "deMulReshape1","deMultiplex","deMulReshape2",...
                              "reshape1","permute1","reshape2","permute2","reshape3",...
                              "FFT","Shift","absAndPower","meanFrameFFT",...
                              "CutArtifacts",...
                              "meanChFFTRaw","rawPhiExtract","AnalyseRawPhi",...
                              "meanPos","calculateFit","applyFit","meanChFFT","SqrtAbsFittedFFT",...
                              "phiExtract",...
                              "AnalysePhi","CalcMuEff","gather",...
                              "exportMeas","exportSignal","exportDeMultiplexed","exportReshaped","exportFFT","exportUsComp",...
                              "AddSingleSplit", "SplitsMeanFrameFFT", "OverallSplitAnalysis"];
            

            for i=1:length(timeTableNames)
                this.timeTable.(timeTableNames{i}) = 0;
            end
        end
        
        function exportData(this, dataName)
            switch dataName
                case 'meas'
                    this.timeTable.exportMeas = tic;
                    if this.general.export.meas
                        this.curRes.export.meas = gather(this.data);
                    end
                    this.timeTable.exportMeas = toc(this.timeTable.exportMeas);
                case 'signal'
                    this.timeTable.exportSignal = tic;
                    if this.general.export.signal
                        this.curRes.export.signal = gather(this.data);
                    end
                    this.timeTable.exportSignal = toc(this.timeTable.exportSignal);
                case 'deMul'
                    this.timeTable.exportDeMultiplexed = tic;
                    if this.general.export.deMul
                        this.curRes.export.deMul = gather(this.data);
                    end
                    this.timeTable.exportDeMultiplexed = toc(this.timeTable.exportDeMultiplexed);
                case 'reshaped'
                    this.timeTable.exportReshaped = tic;
                    if this.general.export.reshaped
                        this.curRes.export.reshaped = gather(this.data);
                    end
                    this.timeTable.exportReshaped = toc(this.timeTable.exportReshaped);
                case 'FFT'
                    this.timeTable.exportFFT = tic;
                    if this.general.export.fft
                        this.curRes.export.fft = gather(this.data);
                    end
                    this.timeTable.exportFFT = toc(this.timeTable.exportFFT);
                case 'usComp'
                     this.timeTable.exportUsComp = tic;
                    if this.general.export.usCompCmplx
                        this.curRes.export.usCompCmplx = squeeze(this.data(:,:, this.freq.fUSIdx, :)); %quite large
                        this.curRes.export.usComp      = abs(this.curRes.export.usCompCmplx);
                    end
                    this.timeTable.exportUsComp = toc(this.timeTable.exportUsComp);
            end       
        end
        
        function resetAlgoRes(this)
             models = this.getResultsModels();
             this.curRes      = models.base;
             this.result      = models.baseChArr; 
        end

        function resetSplitRes(this)
            models      = this.getResultsModels();
            this.curRes = models.base;
            this.result = models.baseCh; 
        end

        function resetCurSplit(this)
            models        = this.getResultsModels();
            this.curSplit = models.baseChArr;
        end
        
        function models = getResultsModels(this)
            %Generic result Array
            base = struct();

            base.frameAvgPowerFFT = [];
            
            base.frameAvgPowerFFTCut = [];
            
            base.rawChAvgFFT    = [];
            base.rawChAvgFFTStd = [];
            base.rawPhi            = [];

            base.posAvgPowerFFT      = [];
            base.fitModelMat         = [];
            base.fittedPowerFFT      = [];
            base.fittedChAvgPowerFFT = [];
            base.fittedFFT           = [];

            base.phiPreCut = [];
            base.phi       = [];
            base.phiLog    = [];
            base.phiNorm   = [];
            
            base.analysis = [];           
            
            base.export.meas           = [];
            base.export.signal         = [];
            base.export.deMul          = [];
            base.export.reshaped       = [];
            base.export.fft            = [];
            base.export.usCompCmplx    = [];
            base.export.usComp         = [];

            %Split Result Struct and Array
            baseSplitArr(this.general.splitNum)  = base;
            baseChArr(this.general.analysisReps) = base;
            
            baseSplit = base;
            baseSplit.splitRes = baseSplitArr;
            
            % Separate Ch Analysis Extension
            baseChSplitArr(this.general.analysisReps, :) = baseSplitArr;
            baseCh(this.general.analysisReps)            = baseSplit;

            % Collect models:
            models.base         = base;
            models.baseSplitArr = baseSplitArr;
            models.baseChArr    = baseChArr;
            models.baseSplit    = baseSplit;
            models.baseCh       = baseCh;

            models.baseChSplitArr = baseChSplitArr;
        end

        function copyStruct(this, source, dest, dstIdx)
            names = fieldnames(source);
            for i = 1:length(names)
               this.(dest)(dstIdx).(names{i}) = source.(names{i});                
            end
        end
    end
end
