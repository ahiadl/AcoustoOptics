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

        measLimit
        timeTable;
        %         reshapedData;
        
        data;
        measSamples;

        curRes;
        result;
        dataInSplit;
        resSplit;
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
            uVars.useFrames         = false;

            uVars.useGPU              = false;
            uVars.useHadamard         = false;
            uVars.contHadamard        = false;
            uVars.highResAO           = false;
            uVars.analyzeSingleCh     = false;
            uVars.contSpeckleAnalysis = false;
            
            uVars.cutArtfct           = false;
            uVars.artfctIdxVec        = [];
            
            uVars.export.meas           = false;
            uVars.export.signal         = false;
            uVars.export.deMul          = false;
            uVars.export.reshaped       = false;
            uVars.export.fft            = false;
            uVars.export.usCompCmplx    = false;
        end

    end
    
    methods
        function this = Algo()
            this.initTimeTable();  
            this.measLimit.time = 8;
        end

        % Set/Get Functions
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
           this.uVars.contHadamard          = user.contHadamard;
           this.uVars.highResAO             = user.highResAO;
           this.uVars.analyzeSingleCh       = user.analyzeSingleCh;
           this.uVars.contSpeckleAnalysis   = user.contSpeckleAnalysis;
           this.uVars.cutArtfct             = user.cutArtfct;
           this.uVars.artfctIdxVec          = user.artfctIdxVec;
           this.uVars.export                = user.export;
           this.calcDimensions();
           vars = this.getVars();
        end

        function setMeasLimit(this, limit)
           this.measLimit.time = limit; 
        end
        
        function setRawData (this, measSamples)
            this.measSamples = measSamples;
            if this.general.analyzeSingleCh
                this.measSamples = gather(this.measSamples);
            end
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
        
        %Config Functions
        function calcDimensions(this)
            this.calcGeneral();
            this.calcGeometry();
            
            this.calcSamplingClkDim();
            this.calcUsSignalDim();
            this.calcSamples();
%             this.calcDigitizer();
            this.calcFreq();
            this.calcTiming();
            this.calcLength();

%             this.allocateReshapeMemory()
        end

        function calcGeneral(this)  
            useGPU              = this.uVars.useGPU;
            useHadamard         = this.uVars.useHadamard;
            contHadamard        = this.uVars.contHadamard;
            highResAO           = this.uVars.highResAO;
            useFrame            = this.uVars.useFrame;
            analyzeSingleCh     = this.uVars.analyzeSingleCh;
            channels            = this.uVars.channels;
            contSpeckleAnalysis = this.uVars.contSpeckleAnalysis;
            cutArtfct           = this.uVars.cutArtfct;
            
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
            acqLimitT         = this.measLimit.time;
            
            splitMeas = timeToSampleTotal > acqLimitT;

            if splitMeas
                %Calculate time to sample
                splitTime = acqLimitT+1;
                splitNum = 0;
                while splitTime > acqLimitT
                    splitNum = splitNum+1;
                    splitTime = round(timeToSampleTotal/splitNum, 3);
                end
                timeToSample = splitTime;
            else
                splitNum     = 1;
                timeToSample = timeToSampleTotal; 
            end
            
            % Export Data
            export.meas           = this.uVars.export.meas ;
            export.signal         = this.uVars.export.signal      && ~splitMeas;
            export.deMul          = this.uVars.export.deMul       && ~splitMeas && ~contHadamard;
            export.reshaped       = this.uVars.export.reshaped    && ~splitMeas && ~contSpeckleAnalysis;
            export.fft            = this.uVars.export.fft         && ~splitMeas;
            export.usCompCmplx    = this.uVars.export.usCompCmplx && ~splitMeas; 
            
            %--------------------
            % Collect parameters
            %--------------------
            clear ('acqLimitT');
               
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
            
            fs     = this.uVars.fs;
            fgClk  = this.uVars.fgClk;
            
            % The data loaded to the AFG must be a multiplication of 16 in manners of data length.
            % creating sClk with 16 cycles meets this requirement no matter
            % the fFgClk used.
            sClkCycles  = 16;
            fsNaive     = fs;
            sClkDutyCyc = 50/100;
            
            fgClkSamplesPerSClkCyc           = fgClk/fs;
            
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
            sinSamples = sin(2*pi * (fPulse/fgClk) * samplesVec);
            
            if ~this.uVars.useHadamard
                sigVec = zeros(1, pulsePerSqnc);
                sigVec(1) = 1;
            else
                sigVec = sMatrix(1, :);
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
                inPhantomPropSamples  = samplesPerSqnc - samplesPerPulse;
%                 postSignalSamples     = samplesPerSqnc;
            end
            
            if this.uVars.contSpeckleAnalysis
                postSignalSamples = postSignalSamples + samplesPerFrameRaw;
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
                sMatInvSingleSqnc = zeros(samplesPerSqnc);
                for i = 1:N
                    sMatInvSingleSqnc( (i-1)*samplesPerPulse+1:i*samplesPerPulse, ...
                                        1:samplesPerPulse:end)         = ...
                                        repmat(sMatInv(i,:),  samplesPerPulse, 1);
                end
                
                for i = 1:samplesPerPulse
                    sMatInvSingleSqnc(i:samplesPerPulse:end, :) = circshift(sMatInvSingleSqnc(i:samplesPerPulse:end, :), (i-1), 2);
                end
                
                sMatInvSingleSqnc = flip(sMatInvSingleSqnc,1);
                
%                 for i = 1:N
%                     sMatInvSingleSqnc( (i-1)*samplesPerPulse+1:i*samplesPerPulse, ...
%                                         samplesPerPulse:samplesPerPulse:end)         = ...
%                                         repmat(sMatInv(i,:),  samplesPerPulse, 1);
%                 end
%                 
%                 for i = 1:samplesPerPulse
%                     sMatInvSingleSqnc(i:samplesPerPulse:end, :) = circshift(sMatInvSingleSqnc(i:samplesPerPulse:end, :), -(i-1), 2);
%                 end
                %----------------------------------------------------------
                % Collect Naive-Hadamard Calculation to Hadamard Parameters
                sMatInvSqnc               = sMatInvSingleSqnc;
                
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
                %---------------------------------
                % Continues Hadamard
                %---------------------------------
                if this.general.contHadamard
                    sMatInvSqncCont = zeros(samplesPerSqnc*samplesPerSqnc, 2*samplesPerSqnc);
                    for i=1:samplesPerSqnc
                        tmp = circshift(sMatInvSingleSqnc, -(i-1), 2);
                        sMatInvSqncCont( (i-1)*samplesPerSqnc+1 : (i*samplesPerSqnc), i:i+samplesPerSqnc-1) = tmp;
                    end

                    sMatInvSqnc = sMatInvSqncCont;
                    
                    %----------------------------------------------------------
                    % Collect Naive-Hadamard Calculation to Hadamard Parameters
                    sqncPerFrameHad          = sqncPerFrameRaw * samplesPerSqnc;
                    framesPerSignalHad        = framesPerSignal;
                    sqncPerSignalHad          = framesPerSignalHad * sqncPerFrameHad;

                    samplesPerFrameHad        = sqncPerFrameHad   * samplesPerSqnc;
                    samplesPerSignalHad       = samplesPerFrameHad * framesPerSignalHad;

                    samplesPerPosPerFrameHad  = samplesPerPulse * sqncPerFrameHad;
                    samplesPerPosPerSignalHad = samplesPerPulse * samplesPerSignalHad;

                    pulsePerFrameHad          = pulsePerSqnc * sqncPerFrameHad;
                    
                    %----------------------------------------------------
                    % Calculate Data Size
                    sqncPerDataOutHad          = (sqncPerDataInHad-1) * samplesPerSqnc;
                    sqncPerDataAllChOutHad     = sqncPerDataOutHad    * channels;
                    samplesPerDataOutHad       = sqncPerDataOutHad    * samplesPerSqnc;
                    samplesPerDataAllChOutHad  = samplesPerDataOutHad * channels;
                    samplesPerDataOutPreCutHad = sqncPerDataInHad * samplesPerSqnc;
                end
                
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
            
            %--------------------------------------
            % ---- Continues Speckle Calculation
            %--------------------------------------
            
            if this.general.contSpeckleAnalysis
                sqncPerDataInCS        = sqncPerData;
                sqncPerDatalAllChInCS  = channels * sqncPerData;
                
                sqncPerFrameCS           = sqncPerFrame;
                samplesPerSqncCS         = samplesPerSqnc;
                framesPerSignalCS        = framesPerSignal   * sqncPerFrameCS;
                sqncPerSignalCS          = framesPerSignalCS * sqncPerFrameCS;
                samplesPerSignalCS       = samplesPerFrame   * framesPerSignalCS;
                samplesPerPosPerSignalCS = samplesPerPulse * sqncPerFrame * framesPerSignalCS;

                %----------------------------------------------------
                % Collect Cont Parameters
                sqncPerSignal          = sqncPerSignalCS;
                framesPerSignal        = framesPerSignalCS;
                samplesPerSignal       = samplesPerSignalCS;
                samplesPerPosPerSignal = samplesPerPosPerSignalCS;
                
                %----------------------------------------------------
                % calculate Data Size
                samplesPerDataOutPreCutCS = sqncPerDataInCS * sqncPerFrameCS * samplesPerSqncCS;
                sqncPerDataOutCS          = (sqncPerDataInCS - sqncPerFrameCS) * sqncPerFrameCS;
                sqncPerDataAllChOutCS     = sqncPerDataOutCS     * channels;
                samplesPerDataOutCS       = sqncPerDataOutCS     * samplesPerSqncCS;
                samplesPerDataAllChOutCS  = samplesPerDataOutCS  * channels;
                
                samplesPerData = samplesPerDataOutCS;
                sqncPerData    = sqncPerDataOutCS;
            end
                        
            if this.general.cutArtfct
                numOfPos = numOfPosAlgo - length(artfctIdxVec);
            end
            
            %--------------------
            % Collect All parameters
            %--------------------
            if  this.general.useHadamard
                this.hadamard.sMatInv = sMatInv;
                this.hadamard.sMatInvSingleSqnc = sMatInvSingleSqnc;
                this.hadamard.sMatInvSqnc       = sMatInvSqnc;
                if this.general.contHadamard
                    this.hadamard.sMatInvSqncCont   = sMatInvSqncCont;
                end
            end
            
            clear ('c', 'channels', 'fSin', 'fSqnc', 'fs', 'sinPerPulse', ... 
                   'frameTime', 'timeToSample', 'preTriggerSamples', ...
                   'distFromPhantom', 'sMatInvSqnc', 'sMatInvSingleSqnc', 'sMatInvSqncCont');
               
            vars = who();
            for i = 1:length(vars)
                if strcmp(vars{i}, 'this'); continue; end
                this.samples.(vars{i}) = eval(vars{i});
            end
        end

        function calcDigitizer(this)
            channels          = this.general.channelsToSample;
            bufferSizeBytes   = this.uVars.bufferSizeBytes;
            bytesPerSample    = this.uVars.bytesPerSample;
            
            samplesPerMeas = this.samples.samplesPerMeas;
            
            samplesPerBufferAllCh = bufferSizeBytes / bytesPerSample;
            samplesPerBuffer      = samplesPerBufferAllCh / channels;
            numOfBuffers          = ceil(samplesPerMeas / samplesPerBuffer);
            samplesPerAcq         = samplesPerBuffer * numOfBuffers;
            samplesPerAcqAllCh    = samplesPerAcq * channels;
            
            %--------------------
            % Collect parameters
            %--------------------            
            vars = who();
            for i = 1:length(vars)
                if strcmp(vars{i}, 'this'); continue; end
                this.digitizer.(vars{i}) = eval(vars{i});
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
            harIdxs  = [1 : samplesPerUsHar : N+1];
            harIdxs(harIdxs == fDcIdxShift) = [];
            numOfHar = length(harIdxs);
            
            % Calc indeces of Harmonics
            envHarLen     = ceil(envUS / df);
            envVec        = -envHarLen : 1 : envHarLen;
            envMat        = repmat(envVec', 1, numOfHar);
            envHarFullLen = length(envVec);
            
            harEnvIdxMat = repmat(harIdxs, envHarFullLen, 1);
            harEnvIdxs   = envMat + harEnvIdxMat;
            harEnvIdxs   = harEnvIdxs(:)';
            harEnvIdxs(harEnvIdxs <= 0 | harEnvIdxs > N) = [];
            
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
            
            if ~this.general.splitMeas
                tVecPulse  = ( 0 : 1 : (samplesPerPulse   - 1) ) * dts;
                tVecSqnc   = ( 0 : 1 : (samplesPerSqnc   - 1) ) * dts;
                tVecFrame  = ( 0 : 1 : (samplesPerFrame   - 1) ) * dts;
                tVecSig    = ( 0 : 1 : (samplesPerSignal - 1) ) * dts;
                
                tVecPosPerFrame  = ( 0 : 1 : (samplesPerPosPerFrame - 1) ) * dts;
                
                tVecMeas   = ( 0 : 1 : (samplesPerMeas    - 1) ) * dts;
            else
                tVecPulse        = [];
                tVecSqnc         = [];
                tVecFrame        = [];
                tVecSig          = [];
                tVecPosPerFrame  = [];
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
            depthIdx    = 1 : 1 : numOfPosAlgo;
            depthIdxLen = numOfPos;

            depthNorm   = depthVec;
            depthCntr   = depthNorm - mean(depthNorm);
            depthZero   = depthNorm - min(depthNorm);
            depthIdx    = 1 : 1 : depthIdxLen;
            
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
        
        % Algorithm Functions
        function res = analyse(this) 
            this.resetAlgoRes();
            if this.general.analyzeSingleCh
                for i = 1:this.general.channelsToSample
                    this.data = this.measSamples (i, :);
                    if this.general.useGPU
                        this.data = gpuArray(this.data);
                    end
                    this.analyseAlgo();
                    this.result(i) = this.curRes;
                end
            else
                % In regular all-ch gpu based analysis measSamples will be
                % upload to GPU externally
                this.data = this.measSamples;
                this.analyseAlgo();
                this.result = this.curRes;
            end
            
            res = this.result;

            % All data in res is CPU
            this.data = [];
            if this.general.splitMeas
                this.curRes    = [];
                this.result = [];
            end
        end
        
        function analyseAlgo(this)
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

            % Global Statistics
            this.timeTable.overallStatistics = tic;
            this.collectGlobalStatistics();
            this.timeTable.overallStatistics = toc(this.timeTable.overallStatistics);
            
            % DeMultiplex
            this.timeTable.overallDemultiplex = tic;
            if this.general.useHadamard
                this.demultiplexSignal();
            end
            this.timeTable.overallDemultiplex = toc(this.timeTable.overallDemultiplex);
            
            % Continues Speckle
            if this.general.contSpeckleAnalysis
                this.contSpeckle();                
            end
            
            % Reshape Data
            this.timeTable.overallReshape = tic;
            this.reshapeSignal();
            this.timeTable.overallReshape = toc(this.timeTable.overallReshape);

            % Signal Processing
            this.timeTable.overallSignalProcessing = tic;
            this.signalProcessing();
            this.timeTable.overallSignalProcessing = toc(this.timeTable.overallSignalProcessing);

            this.timeTable.FullAnalysis = toc(this.timeTable.FullAnalysis);
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
        
        function collectGlobalStatistics(this)
            this.curRes.std       = gather(std(this.data, 0, 2));
            this.curRes.avg       = gather(mean(this.data, 2));
            this.curRes.normNoise = mean(this.curRes.std ./ this.curRes.avg);
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
            
            if this.general.contHadamard
                this.data = repmat(this.data, 2, 1);
                this.data(samplesPerSqncHad:end, :) = circshift(this.data(samplesPerSqncHad:end, :), -1, 2);
            end
            
            this.data = sMatInvSqnc * this.data;
            
            this.timeTable.deMulReshape2 = tic;
            this.data = reshape(this.data, samplesPerDataOutPreCutHad, channels)';
            this.timeTable.deMulReshape2 = toc(this.timeTable.deMulReshape2);

            this.data(:, samplesPerDataOutHad+1:end) = [];
            
            this.exportData ('deMul')
        end
        
        function contSpeckle(this)
            % Collect Variables
            channels                  = this.general.channelsToAnalyze;
            samplesPerSqncCS          = this.samples.samplesPerSqncCS;
            sqncPerDatalAllChInCS     = this.samples.sqncPerDatalAllChInCS;
            sqncPerFrameCS            = this.samples.sqncPerFrameCS;
            samplesPerDataOutPreCutCS = this.samples.samplesPerDataOutPreCutCS;
            samplesPerDataOutCS       = this.samples.samplesPerDataOutCS;
            
            % Create Continues Frames
            this.timeTable.csReshape1 = tic;
            this.data = reshape(this.data', samplesPerSqncCS, sqncPerDatalAllChInCS);
            this.timeTable.csReshape1 = toc(this.timeTable.csReshape1);
            
            this.data = repmat(this.data, sqncPerFrameCS, 1);
            for i = 1 : sqncPerFrameCS
                this.data((i-1)*samplesPerSqncCS+1:i*samplesPerSqncCS, :) =...
                    circshift(this.data((i-1)*samplesPerSqncCS+1:i*samplesPerSqncCS, :), -(i-1), 2);
            end
            
            this.timeTable.csReshape2 = tic;
            this.data = reshape(this.data, samplesPerDataOutPreCutCS, channels)';
            this.timeTable.csReshape2 = toc(this.timeTable.csReshape2);
            
            this.data(:, samplesPerDataOutCS+1:end) = [];
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
            this.timeTable.reshape4 = tic;
            this.data = reshape(this.data, channels,...
                                           framesPerSignal,...
                                           samplesPerPosPerFrame,...
                                           numOfPos);
            this.timeTable.reshape4 = toc(this.timeTable.reshape4);                       
            
            % Bring ch dim to first
            % [ch x frames x samples x numOfPos]
%             this.timeTable.permute3 = tic;
%             this.data = permute(this.data, [2,1,3,4]);
%             this.timeTable.permute3 = toc(this.timeTable.permute3); 
            
            this.exportData ('reshaped');
        end
        
        function signalProcessing(this)
            % data(input)              - [ch x frames x samples x numOfPos]]
            % qAvgChFFT(output)        - [ch x samplesPerPos x numOfPos]
            
            samplesPerPosPerFrame = this.samples.samplesPerPosPerFrame;
            
            % Perform Fourier Transform
            this.timeTable.FFTandShift = tic;
            this.data = (2./samplesPerPosPerFrame) * fftshift( fft(this.data, [], 3) , 3);
            this.timeTable.FFTandShift = toc(this.timeTable.FFTandShift);

            this.exportData('FFT')
            this.exportData('usComp')

            % Absolute and square
            this.timeTable.absAndPower = tic;
            this.data                  = abs(this.data).^2;
            this.timeTable.absAndPower = toc(this.timeTable.absAndPower);
            
            % Average over all frames
            this.timeTable.meanFrameFFT = tic;
            this.curRes.qAvgChFFT          = permute(mean(this.data, 2), [1,3,4,2]);
            this.curRes.qStdChFFT          = permute(std(this.data, 0, 2), [1,3,4,2]);
            this.timeTable.meanFrameFFT = toc(this.timeTable.meanFrameFFT);
           
            if ~this.general.splitMeas
                this.extractPhiFromFFT();
            else
                this.curRes.qAvgChFFT = gather(this.curRes.qAvgChFFT);
                this.curRes.qStdChFFT = gather(this.curRes.qStdChFFT);
            end
        end
        
        function resSplit = analyseSplittedData(this, dataInSplit)
            % usCompCmplx(output)      - [quants x ch x numOfPos]
            % qAvgFFT(output)          - [ch x samplesPerPos x numOfPos]
            % unFittedFFT(output)      - [samplesPerPos x numOfPos]
            % unFittedFFTShift(output) - [samplesPerPos x numOfPos]
            % fitModel(output)         - [samplesPerPos x numOfPos]
            % fittedFFT(output)        - [samplesPerPos x numOfPos]
            % phi(output)              - [1 x numOfPos]
            
            this.timeTable.splitAnalysis = tic;
            fprintf("ALGO: Analyzing splitted data.\n");
            this.resetAlgoRes();
            this.dataInSplit = dataInSplit;
            splitNum = this.general.splitNum;
            
            splitQAvgFFT = zeros(this.general.channelsToAnalyze, this.samples.samplesPerPosPerFrame, this.samples.numOfPosAlgo, splitNum);
            splitQStdFFt = zeros(this.general.channelsToAnalyze, this.samples.samplesPerPosPerFrame, this.samples.numOfPosAlgo, splitNum);
            
            for j = 1:this.general.analysisReps
                this.timeTable.splitAnalysis = tic;
                this.timeTable.meanQuantFFT  = tic;
                for i = 1:splitNum
                    splitQAvgFFT(:,:,:,i) = gather(dataInSplit(j,i).qAvgChFFT);
%                     splitQStdFFt(:,:,:,i) = gather(dataInSplit(j,i).qStdChFFT);
                end
                this.curRes.qAvgChFFT   = mean(splitQAvgFFT, 4);
                this.curRes.qStdChFFT   = (1/splitNum)*sqrt(sum(splitQStdFFt.^2 ,4));
                this.timeTable.meanQuantFFT = toc(this.timeTable.meanQuantFFT);

                this.extractPhiFromFFT();
                this.setStruct(this.curRes, 'resSplit', j);
                this.resSplit(j).totalRes = this.curRes;
                this.timeTable.splitAnalysis = toc(this.timeTable.splitAnalysis);
            end
            resSplit = this.resSplit;
        end
        
        function extractPhiFromFFT(this)
            % Input:
            % qAvgChFFT - [ch x samplesPerPos x numOfPos]
            % Outputs:
            % unFittedAvgFFT    - [samplesPerPos x numOfPos]
            % unFittedAvgFFTStd - [samplesPerPos x numOfPos]
            % rawPhi            - [1 x numOfPos]       
            % posAvgChFFT       - [ch x samplesPerPos]
            % fitModelMat       - [ch x samplesPerPos]
            % fittedQAvgFFT     - [ch x samplesPerPos x numOfPos]
            % fittedChAvgFFT    - [samplesPerPos x numOfPos]
            % fittedFFT         - [samplesPerPos x numOfPos]
            % fittedFFTNorm     - [samplesPerPos x numOfPos]
            % phi               - [1 x numOfPos]      
            
            % Collect Variables
            channels = this.general.channelsToAnalyze;
            fUsIdx = this.freq.fUSIdx;
            
            %-----------------------------------------------------
            % Raw Phi Calculation (no fit - averaging over all ch)
            %-----------------------------------------------------
            % Average over all channels and SQRT(!)
            this.timeTable.meanChFFTRaw = tic;
            this.curRes.unFittedAvgFFT     = sqrt(permute(mean(this.curRes.qAvgChFFT, 1), [2,3,1]));
            this.curRes.unFittedAvgFFTStd  = sqrt(permute((1/channels)*sqrt(sum(this.curRes.qStdChFFT.^2,1)),[2,3,1]));
            this.timeTable.meanChFFTRaw = toc(this.timeTable.meanChFFTRaw);
            
            this.curRes.rawPhi = this.curRes.unFittedAvgFFT(fUsIdx, :);

            %-----------------------------------------------------
            % Phi Calculation (fit each channel and then averag)
            %-----------------------------------------------------
            %Averag over all pos to create avg channel response before fit
            this.timeTable.meanPos = tic;
            this.curRes.posAvgChFFT    = mean(this.curRes.qAvgChFFT, 3);
            this.timeTable.meanPos = toc(this.timeTable.meanPos);
            
            % Calculate fit model to FFT
            this.timeTable.calculateFit = tic;
            this.curRes.fitModelMat = this.calcFittedFFT();
            this.timeTable.calculateFit = toc(this.timeTable.calculateFit);

            % Apply Fit
            this.timeTable.applyFit = tic;
            this.curRes.fittedQAvgFFT = this.curRes.qAvgChFFT ./ this.curRes.fitModelMat;
            this.timeTable.applyFit = toc(this.timeTable.applyFit);

            % Average over all channels
            this.timeTable.meanChFFT = tic;
            this.curRes.fittedChAvgFFT = permute(mean(this.curRes.fittedQAvgFFT, 1), [2,3,1]);
            this.timeTable.meanChFFT = toc(this.timeTable.meanChFFT);

            % Sqrt Average FFT
            this.timeTable.sqrtFittedFFT = tic;
            this.curRes.fittedFFT = sqrt(this.curRes.fittedChAvgFFT);
            this.timeTable.sqrtFittedFFT = toc(this.timeTable.sqrtFittedFFT);

            %Normalize FFT
            this.timeTable.normalizeFFT = tic;
            this.curRes.fittedFFTNorm = this.curRes.fittedFFT - 1;
            this.timeTable.normalizeFFT = toc(this.timeTable.normalizeFFT);    

            % Extract Ultrasound Component
            this.timeTable.phiExtract = tic;
            this.curRes.phiPreCut     = this.curRes.fittedFFTNorm(fUsIdx, :);
            this.timeTable.phiExtract = toc(this.timeTable.phiExtract);

            this.curRes.phi = this.curRes.phiPreCut;
            
            %Cut Artifacts
            this.cutArtifactIndexes();
            
            %Normalize and Log
            this.normAndLog();
            
            % Gather from GPU
            this.timeTable.gather = tic;
            
            this.curRes.qAvgChFFT = gather(this.curRes.qAvgChFFT);
            this.curRes.qStdChFFT = gather(this.curRes.qStdChFFT); 
            
            this.curRes.unFittedAvgFFT    = gather(this.curRes.unFittedAvgFFT);
            this.curRes.unFittedAvgFFTStd = gather(this.curRes.unFittedAvgFFTStd);
            this.curRes.rawPhi            = gather(this.curRes.rawPhi);
            
            this.curRes.posAvgChFFT    = gather(this.curRes.posAvgChFFT);
            this.curRes.fittedQAvgFFT  = gather(this.curRes.fittedQAvgFFT);
            this.curRes.fittedChAvgFFT = gather(this.curRes.fittedChAvgFFT);
            this.curRes.fitModelMat    = gather(this.curRes.fitModelMat);
            this.curRes.fittedFFT      = gather(this.curRes.fittedFFT);
            this.curRes.fittedFFTNorm  = gather(this.curRes.fittedFFTNorm);
            this.curRes.phiPreCut      = gather(this.curRes.phiPreCut);
            this.curRes.phi            = gather(this.curRes.phi);
            this.curRes.phiNorm        = gather(this.curRes.phiNorm);
            this.curRes.phiLog         = gather(this.curRes.phiLog);
            
            this.timeTable.gather        = toc(this.timeTable.gather);

            % Calculate SNR
            this.timeTable.calcSNR = tic;
            this.curRes.SNR = this.calcSNR();
            this.timeTable.calcSNR = toc(this.timeTable.calcSNR);
        end

        function fitModelMat = calcFittedFFT(this)
            fitModelMat = ones(size(this.curRes.posAvgChFFT));
            if this.general.useGPU
                fitModelMat = gpuArray(fitModelMat);
            end
            if length(this.freq.fitIdxShift) <11
                return;
            end
            
            % Collect Variables
            channels  = this.general.channelsToAnalyze;
            numOfPos  = this.samples.numOfPosAlgo;
            fBar      = this.freq.fBar;
            idxs      = this.freq.fitIdxShift;
            
            % Calculate Fit
            for i = 1:channels
                fitVec = this.curRes.posAvgChFFT(i, idxs);
                fitModelMat(i,:) = interp1(fBar(idxs), fitVec, fBar);
            end
            
%             figure(); 
%             stem(fBar(idxs), ones(1,length(fBar(idxs)))); hold on
%             yyaxis right
%             plot(fBar, squeeze(this.curRes.qAvgChFFT(1,:,30)))

            fitModelMat = repmat(fitModelMat, 1, 1, numOfPos);
        end
        
        function [SNR] = calcSNR(this)
            [maxVal, maxIdx] = max(this.curRes.phi);
            avgVal = mean(this.curRes.phi);
            
            bkgPts = find (this.curRes.phi < (avgVal - 0.01*avgVal));
            tailIdx =  bkgPts(bkgPts > maxIdx);
            if isempty(tailIdx)
                tailIdx = maxIdx + 2;
            end
            
            SNR.tailIdx  = tailIdx(1);
            SNR.peak     = maxVal;
            SNR.tailAvg  = mean(this.curRes.phi(tailIdx(1):end));
            SNR.noiseStd = std(this.curRes.phi(tailIdx(1):end));
            SNR.val      = maxVal/std(this.curRes.phi(tailIdx(1):end));
            
            [maxVal, maxIdx] = max(this.curRes.rawPhi);
            avgVal = mean(this.curRes.rawPhi);
            
            bkgPts = find (this.curRes.rawPhi < (avgVal - 0.01*avgVal));
            tailIdx =  bkgPts(bkgPts > maxIdx);
            if isempty(tailIdx)
                tailIdx = maxIdx + 2;
            end
            
            noiseAvg   =  mean(this.curRes.rawPhi(tailIdx(1):end));
            normRawPhi =  this.curRes.rawPhi - noiseAvg;
            maxVal = max(normRawPhi);
            
            SNR.tailIdxRaw = tailIdx(1);
            SNR.peakRaw = maxVal;
            SNR.noiseStdRaw = std(normRawPhi(tailIdx(1):end));
            SNR.valRaw = maxVal/std(normRawPhi(tailIdx:end));
            SNR.tailAvgRaw = noiseAvg;
        end
        
        function cutArtifactIndexes(this)
            if this.general.cutArtfct
                this.curRes.phi(this.samples.artfctIdxVec) = [];
            end
        end
        
        function normAndLog(this)
             minVal = min(this.curRes.phi);
             maxVal = max(this.curRes.phi);
             span = maxVal - minVal;
             
             this.curRes.phiNorm = (this.curRes.phi - minVal) ./ span;
             this.curRes.phiLog = log(this.curRes.phiNorm);         
        end
        
        % Misc
        function initTimeTable(this)
            this.timeTable.FullAnalysis             = 0;
            
            this.timeTable.overallNetSignal         = 0;
            this.timeTable.overallStatistics        = 0;
            this.timeTable.overallDemultiplex       = 0;
            this.timeTable.overallReshape           = 0;
            this.timeTable.overallSignalProcessing  = 0;
            
            this.timeTable.extractNetSignal         = 0;
            
            this.timeTable.deMulReshape1            = 0;
            this.timeTable.deMultiplex              = 0; 
            this.timeTable.deMulReshape2            = 0;
            
            this.timeTable.reshape1                 = 0;
            this.timeTable.permute1                 = 0;
            this.timeTable.reshape2                 = 0;
            this.timeTable.reshape3                 = 0;
            this.timeTable.permute2                 = 0;
            this.timeTable.reshape4                 = 0;
            this.timeTable.permute3                 = 0;
            
            this.timeTable.FFTandShift              = 0;
            this.timeTable.absAndPower              = 0;
            this.timeTable.meanFrameFFT             = 0;
            this.timeTable.meanChFFTRaw             = 0;
            this.timeTable.meanPos                  = 0;
            this.timeTable.calculateFit             = 0;
            this.timeTable.applyFit                 = 0;
            this.timeTable.meanChFFT                = 0;
            this.timeTable.sqrtFittedFFT            = 0;
            this.timeTable.normalizeFFT             = 0;
            this.timeTable.phiExtract               = 0;
            
            this.timeTable.gather                   = 0;
            
            this.timeTable.calcSNR                  = 0;
            
            this.timeTable.exportRawData            = 0;
            this.timeTable.exportNetSignal          = 0;
            this.timeTable.exportDeMultiplexed      = 0;
            this.timeTable.exportReshaped           = 0;
            this.timeTable.exportFFT                = 0;
            this.timeTable.exportUsComp             = 0;
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
             [this.curRes, this.result, this.dataInSplit, this.resSplit]= this.getResultsStruct;
        end
        
        function [res, result, dataInSplit, resSplit] = getResultsStruct(this)
            %Generic result Array
            tmpRes = struct();
            
            tmpRes.std       = [];
            tmpRes.avg       = [];
            tmpRes.normNoise = [];
            
            tmpRes.qAvgChFFT = [];
            tmpRes.qStdChFFT = [];
            
            tmpRes.unFittedAvgFFT    = [];
            tmpRes.unFittedAvgFFTStd = [];
            tmpRes.rawPhi            = [];

            tmpRes.posAvgChFFT    = [];
            tmpRes.fitModelMat    = [];
            tmpRes.fittedQAvgFFT  = [];
            tmpRes.fittedChAvgFFT = [];
            tmpRes.fittedFFT      = [];
            tmpRes.fittedFFTNorm  = [];
            tmpRes.phiPreCut      = [];
            tmpRes.phi            = [];
            tmpRes.phiNorm        = [];
            tmpRes.phiLog         = [];

            tmpRes.SNR = [];            
            
            tmpRes.export.meas           = [];
            tmpRes.export.signal         = [];
            tmpRes.export.deMul          = [];
            tmpRes.export.reshaped       = [];
            tmpRes.export.fft            = [];
            tmpRes.export.usCompCmplx    = [];
            tmpRes.export.usComp         = [];
            
            res = tmpRes;
            
            %Split Result Struct and Array
            tmpSplitArr(this.general.splitNum) = tmpRes; 
            tmpSplitRes = tmpRes;
            tmpSplitRes.splitRes = tmpSplitArr;
            
            % Single Ch analysis extension 
            result(this.general.analysisReps)         = tmpRes;
            dataInSplit(this.general.analysisReps, :) = tmpSplitArr;
            resSplit(this.general.analysisReps)       = tmpSplitRes;
        end
        
        function setStruct(this, source, dest, dstIdx)
            names = fieldnames(source);
            for i = 1:length(names)
               this.(dest)(dstIdx).(names{i}) = source.(names{i});                
            end
        end
    end
end