classdef PDExp < handle
    %PDEXP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        app
        stagesT
        ao
        vars
        data
    end
    
    methods (Static)
        
        function uVars = createUserVars()
            %AO
            uVars.quant           = 0.002;
            uVars.timeToSample    = 2;
            uVars.distFromPhantom = 8.5e-2;
            uVars.fs              = 5e6;
            uVars.useGPU          = true;
            uVars.exportRawData   = true;
            uVars.idxToPlot       = 20;
            uVars.useHadamard     = true;
            
            uVars.fwl = 0.5;
            
            uVars.usePhiPP = true;
            uVars.refSNR = 100;
            
            % Extract Net Signal
            uVars.peakEnv = 20;
            uVars.peakIdx = 20;
            uVars.cutIdxs = [];
            
            % Power
            uVars.netSignal = 200e-9;
            uVars.minPwr    = 4e-6 ;
            uVars.maxPwr    = 5e-3;
            
            % Measurement
            uVars.reps       = 30;
            uVars.stabTime   = 1;
            uVars.subSetSize = 10;
            uVars.measType   = 'PMT';

            uVars.minTheta  = 58.8;
            uVars.thetaVec  = [];
            uVars.refFactor = [];
            uVars.voltage   = []; 
            uVars.rotate    = true;
            
            uVars.updateAOGui = false;
            uVars.logScale    = false;
            
            % File System
            uVars.measName = 'pdExp';
            uVars.dir      = 'D:\Photodiode';
            uVars.saveFlag = false;
            uVars.saveN    = 10;
            
            uVars.tgUpdate = false;
        end
        
        function in = dialogBox(var)
            beep()
            prompt = {sprintf('Change the Gain of the PMT to %d.\n Enter its DC level(mV):', var)};
            dlgtitle = 'Gain Change';
            definput = {'0'};
            dims = [1 40];
        %     opts.Interpreter = 'tex';
            opts = struct();
            answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
            in = str2double(answer);
        end
        
    end
    
    methods
        function this = PDExp(ao, stagesIn)           
            if nargin == 2
                this.ao = ao;
                this.stagesT = stagesIn;
            else
                app = AOGUI();
                this.ao = app.ao;
                this.stagesT = stages('Zaber', 'COM22');
                this.stagesT.connect();
                this.stagesT.home('X');
            end
            this.data = struct();

            this.vars.figs.single.hFig      = gobjects(1, 1);
            this.vars.figs.single.axRaw     = [];
            this.vars.figs.single.axPhiPost = [];
            this.vars.figs.single.axStats   = [];
            this.vars.figs.single.axTrend   = [];
            this.vars.figs.single.axErr     = [];
            this.vars.figs.single.axErrLog  = [];
            this.vars.figs.single.axRStd    = [];
            this.vars.figs.single.axFFT     = [];

            this.vars.figs.dual.hFig         = gobjects(1, 3);
            this.vars.figs.dual.axRawPD      = [];
            this.vars.figs.dual.axRawPMT     = [];
            this.vars.figs.dual.axPostCutPD  = [];
            this.vars.figs.dual.axPostCutPMT = [];
            this.vars.figs.dual.axFFTPD      = [];
            this.vars.figs.dual.axFFTPMT     = [];
            this.vars.figs.dual.axStatsPD    = [];
            this.vars.figs.dual.axStatsPMT   = [];
            this.vars.figs.dual.axTrendPD    = [];
            this.vars.figs.dual.axTrendPMT   = [];
            this.vars.figs.dual.axStd        = [];
            this.vars.figs.dual.axGain       = [];
            this.vars.figs.dual.axGainLog    = [];
        end
        
        function setVars(this, uVars)
            this.vars.uVars = uVars;
                        
            %AO Vars
            this.vars.ao = this.ao.getVars();
            this.vars.uVarsAO = this.vars.ao.uVars;
            
            this.vars.uVarsAO.ao.frametTime         = uVars.frame;
            this.vars.uVarsAO.ao.timeToSample       = uVars.timeToSample;
            this.vars.uVarsAO.ao.distFromPhantom    = uVars.distFromPhantom;
            this.vars.uVarsAO.ao.fs                 = uVars.fs;
            this.vars.uVarsAO.ao.fSqnc              = uVars.fSqnc;
            this.vars.uVarsAO.ao.useGPU             = uVars.useGPU;
            this.vars.uVarsAO.ao.exportData.meas    = uVars.exportMeas;
            this.vars.uVarsAO.ao.useHadamard        = uVars.useHadamard;
            
            this.vars.uVarsAO.ao.channels        = 1;
            this.vars.uVarsAO.ao.analyzeSingleCh = false;
            this.vars.dualOP = false;
            
             if strcmp(uVars.measType, 'Dual')
                this.vars.uVarsAO.ao.channels        = 2;
                this.vars.uVarsAO.ao.analyzeSingleCh = true;
                this.vars.dualOP                     = true;
             end

            this.vars.uVarsAO.ao.cutArtfct = false;
            this.vars.uVarsAO.ao.artfctIdxVec = [];
             
            % AO Figures
            this.vars.uVarsAO.figs.depthIdx = uVars.idxToPlot;
            this.vars.updateAOGui           = uVars.updateAOGui;
            
            this.vars.uVarsAO.figs.validStruct.usSignal            = false;
            this.vars.uVarsAO.figs.validStruct.sClk                = false;
            
            this.vars.uVarsAO.figs.validStruct.meas                = false;
            this.vars.uVarsAO.figs.validStruct.signal              = false;
            this.vars.uVarsAO.figs.validStruct.deMul               = false;
            this.vars.uVarsAO.figs.validStruct.reshaped            = false;

            this.vars.uVarsAO.figs.validStruct.unFittedFFTAllCh    = false;
            this.vars.uVarsAO.figs.validStruct.unFittedFFTSingleCh = false;
            this.vars.uVarsAO.figs.validStruct.fittedFFTAllCh      = false;
            this.vars.uVarsAO.figs.validStruct.fittedFFTSingleCh   = false;
            this.vars.uVarsAO.figs.validStruct.avgFFT              = false;
            this.vars.uVarsAO.figs.validStruct.normFFT             = false;
            
            this.vars.uVarsAO.figs.validStruct.phi                 = false;
            this.vars.uVarsAO.figs.validStruct.rawPhi              = false;
            
            if this.vars.updateAOGui
                this.vars.uVarsAO.figs.validStruct.unFittedFFTAllCh    = this.ao.graphics.figs.validStruct.unFittedFFTAllCh;
                this.vars.uVarsAO.figs.validStruct.unFittedFFTSingleCh = this.ao.graphics.figs.validStruct.unFittedFFTSingleCh;
                this.vars.uVarsAO.figs.validStruct.fittedFFTAllCh      = this.ao.graphics.figs.validStruct.fittedFFTAllCh;
                this.vars.uVarsAO.figs.validStruct.fittedFFTSingleCh   = this.ao.graphics.figs.validStruct.fittedFFTSingleCh;
                this.vars.uVarsAO.figs.validStruct.avgFFT              = this.ao.graphics.figs.validStruct.avgFFT;
                this.vars.uVarsAO.figs.validStruct.normFFT             = this.ao.graphics.figs.validStruct.normFFT;
                this.vars.uVarsAO.figs.validStruct.phi                 = this.ao.graphics.figs.validStruct.phi;
                this.vars.uVarsAO.figs.validStruct.rawPhi              = this.ao.graphics.figs.validStruct.rawPhi;
            end
            
            this.vars.fwl = uVars.fwl;
            
            this.ao.setVars(this.vars.uVarsAO);
            this.ao.configPeripherals();
            this.vars.ao = this.ao.getVars();
            
            % PD Exp - Analysis Vars
            this.vars.phiLen = this.vars.ao.measVars.algo.samples.numOfPos;
            this.vars.fBar   = gather(this.vars.ao.measVars.algo.freq.fBar);
            this.vars.fIdx   = this.vars.ao.measVars.algo.freq.fUSIdx;
            this.vars.fftLen = this.vars.ao.measVars.algo.samples.samplesPerPosPerFrame;
            
            this.vars.usePhiPP = uVars.usePhiPP;
            this.vars.refSNR   = uVars.refSNR;
            
            % Extract Net Signal, Peak and Background
            this.vars.peakEnv      = uVars.peakEnv;
            this.vars.peakIdx      = uVars.peakIdx;
            this.vars.bkgIdx       = uVars.bkgIdx;
            this.vars.cutIdxs      = uVars.cutIdxs;

            keepIdxs                    = 1:this.vars.phiLen; 
            keepIdxs(this.vars.cutIdxs) = [];
            this.vars.keepIdxs          = keepIdxs;
            this.vars.cutLen            = length(this.vars.cutIdxs);
            
            this.vars.phiLenPostCut  = length(this.vars.keepIdxs);
            this.vars.sigIdxs        = (this.vars.peakIdx-this.vars.peakEnv):(this.vars.peakIdx+this.vars.peakEnv);
            this.vars.sigIdxs(this.vars.sigIdxs<1) = [];
            this.vars.sigIdxs(this.vars.sigIdxs>this.vars.phiLenPostCut) = [];
            
            this.vars.preSigIdx      = sum( this.vars.cutIdxs < min(this.vars.sigIdxs));
            this.vars.postSigIdx     = sum( this.vars.cutIdxs > max(this.vars.sigIdxs));
            this.vars.peakIdxPostCut = this.vars.peakIdx - this.vars.preSigIdx;
            this.vars.sigIdxPostCut  = (this.vars.peakIdxPostCut-this.vars.peakEnv):(this.vars.peakIdxPostCut+this.vars.peakEnv);
            this.vars.sigIdxPostCut(this.vars.sigIdxPostCut<1) = [];
            this.vars.sigIdxPostCut(this.vars.sigIdxPostCut>this.vars.phiLenPostCut) = [];
            this.vars.bkgIdxPostCut  = [1:this.vars.sigIdxPostCut(1)-1, this.vars.sigIdxPostCut(end)+1:this.vars.phiLenPostCut];
            this.vars.bkgLen         = length(this.vars.bkgIdxPostCut);
            this.vars.sigLen         = length(this.vars.sigIdxPostCut);
            
            bkgIdxs                = 1:this.vars.phiLen; 
            bkgIdxs([this.vars.cutIdxs, this.vars.sigIdxs]) = [];
            this.vars.bkgIdxPreCut = bkgIdxs;
            
            % Power Parameters
            this.vars.pwr.netSignal = uVars.netSignal;
            this.vars.pwr.minPwr    = uVars.minPwr;
            this.vars.pwr.maxPwr    = uVars.maxPwr;
            
            this.vars.pwr.span = this.vars.pwr.maxPwr - this.vars.pwr.minPwr;

            % Scan Parameters (Gain/Reference)
            this.vars.reps       = uVars.reps;
            this.vars.stabTime   = uVars.stabTime;
            this.vars.subSetSize = uVars.subSetSize;
            this.vars.subSetAx   = this.vars.subSetSize:1:this.vars.reps;
            this.vars.subSetLen  = length(this.vars.subSetAx);            
            
            this.vars.measType = uVars.measType;

            switch this.vars.measType
                case {'PD'; 'Dual'}
                    this.vars.minTheta      = uVars.minTheta;
                    this.vars.thetaTrend    = uVars.thetaTrend;
                    this.vars.thetaVec      = uVars.thetaVec;
                    this.vars.pwr.refFactor = uVars.refFactor;
                    if length(this.vars.thetaVec) > 1
                        this.vars.thetaVecNorm  = abs(this.vars.thetaVec - uVars.minTheta); 
%                         this.vars.thetaVecNorm  = this.vars.thetaVecNorm - this.vars.thetaVecNorm(1);
                        this.vars.pwr.Ivec      = ((sin(deg2rad(this.vars.thetaVecNorm*2))).^2)*this.vars.pwr.span+this.vars.pwr.minPwr;
                        this.vars.pwr.refFactor = this.vars.pwr.Ivec/this.vars.pwr.netSignal;
                    else
                        this.vars.pwr.Ivec     =  this.vars.pwr.refFactor .* this.vars.pwr.netSignal;
                        this.vars.thetaVecNorm = rad2deg( asin(  sqrt( (this.vars.pwr.Ivec - this.vars.pwr.minPwr)/this.vars.pwr.span) )/2);
                        switch this.vars.thetaTrend 
                            case 'Positive'
                                this.vars.thetaVec     = this.vars.minTheta + this.vars.thetaVecNorm;
                            case 'Negative'
                                this.vars.thetaVec     = this.vars.minTheta - this.vars.thetaVecNorm;
                        end
                    end
                    this.vars.rotate  = true;
                    this.vars.loopVar = this.vars.thetaVec;
                case 'PMT'
                    this.vars.voltage = uVars.voltage;
                    m              = 0.007;
                    n              = 0.4;
                    log10Gain      = this.vars.voltage*m+n;
                    this.vars.gain = 10.^log10Gain;
                    
                    this.vars.rotate  = false;
                    this.vars.loopVar = this.vars.voltage;
                case {'Point', 'Dual-Point'}
                    this.vars.thetaVec = uVars.thetaVec;
                    if ~isempty(this.vars.thetaVec)
                        this.vars.thetaVecNorm  = abs(this.vars.thetaVec - uVars.minTheta); 
%                         this.vars.thetaVecNorm  = this.vars.thetaVecNorm - this.vars.thetaVecNorm(1);
                        this.vars.pwr.Ivec      = ((sin(deg2rad(this.vars.thetaVecNorm*2))).^2)*this.vars.pwr.span+this.vars.pwr.minPwr;
                        this.vars.pwr.refFactor = this.vars.pwr.Ivec/this.vars.pwr.netSignal;
                    else
                        this.vars.pwr.Ivec     =  this.vars.pwr.refFactor .* this.vars.pwr.netSignal;
                        this.vars.thetaVecNorm = rad2deg( asin(  sqrt( (this.vars.pwr.Ivec - this.vars.pwr.minPwr)/this.vars.pwr.span) )/2);
                        switch this.vars.thetaTrend 
                            case 'Positive'
                                this.vars.thetaVec     = this.vars.minTheta + this.vars.thetaVecNorm;
                            case 'Negative'
                                this.vars.thetaVec     = this.vars.minTheta - this.vars.thetaVecNorm;
                        end
                    end
                    this.vars.rotate   = uVars.rotate;
                    this.vars.loopVar  = 1;
            end 
            this.vars.minLoop = min(this.vars.loopVar);
            this.vars.maxLoop = max(this.vars.loopVar);
            if length(this.vars.loopVar)>1
                this.vars.loopVarDiff = [abs(this.vars.loopVar(1) - this.vars.loopVar(2)),...
                                         abs(this.vars.loopVar(end) - this.vars.loopVar(end-1))];
            else
                this.vars.loopVarDiff = ones(1,2);
            end
            this.vars.loopVarMargins = [this.vars.loopVar(1) - 0.1*this.vars.loopVarDiff(1), this.vars.loopVar(end) + 0.1*this.vars.loopVarDiff(2)];
            
            this.vars.measLen = length(this.vars.loopVar);

            % Init Graphics
            this.vars.figs.logScale = uVars.logScale;

            if this.vars.measLen > length(unique(this.vars.loopVar))
                this.vars.errBarAx = 1:this.vars.measLen;
                this.vars.figs.lims = [0, this.vars.measLen+1];
            else
                this.vars.errBarAx = this.vars.loopVar;
                this.vars.figs.lims = [this.vars.minLoop - this.vars.loopVarDiff(1), this.vars.maxLoop + this.vars.loopVarDiff(2)];
            end

            % File System
            this.vars.measName = uVars.measName;
            this.vars.dir      = uVars.dir;
            this.vars.saveFlag = uVars.saveFlag;
            this.vars.saveN    = uVars.saveN;
                       
            this.vars.measDir = sprintf('%s/%s', this.vars.dir, this.vars.measName);
            
            this.vars.tgUpdate = uVars.tgUpdate;
        end
        
        function vars = getVars(this)
            vars = this.vars;
        end
        
        function initFS(this)
            if ~exist(this.vars.measDir, 'dir')
                mkdir(this.vars.measDir)
            end
        end
        
        function initFigsSingle(this)
            try 
                if ~isgraphics(this.vars.figs.single.hFig)      || ~isgraphics(this.vars.figs.single.axRaw)   || ...
                   ~isgraphics(this.vars.figs.single.axPhiPost) || ~isgraphics(this.vars.figs.single.axStats) || ...
                   ~isgraphics(this.vars.figs.single.axTrend)   || ~isgraphics(this.vars.figs.single.axErr)   || ...
                   ~isgraphics(this.vars.figs.single.axErrLog)  || ~isgraphics(this.vars.figs.single.axRStd)  || ...
                   ~isgraphics(this.vars.figs.single.axFFT)
                
                    if ~isfield(this.vars.figs, 'hFig')      || ~isgraphics(this.vars.figs.hFig)
                        this.vars.figs.hFig = figure();
                        clf
                    end

                    % Current Raw Signal
                    this.vars.figs.axRaw     = subplot(3,3,1);
                    this.vars.figs.hNet      = stem(this.vars.figs.axRaw, this.vars.keepIdxs, zeros(1, this.vars.phiLenPostCut)); hold on
                    if ~isempty(this.vars.cutIdxs)
                        this.vars.figs.hCut = stem(this.vars.figs.axRaw, this.vars.cutIdxs,  zeros(1, this.vars.cutLen), 'r'); hold on
                    else
                        this.vars.figs.hCut = 0;
                    end
                    hold off
                    title(this.vars.figs.axRaw, "Raw Phi")
                    
                    %Current Post-Cutting Signal
                    this.vars.figs.axPhiPost = subplot(3,3,2);
                    this.vars.figs.hSig      = stem(this.vars.figs.axPhiPost, this.vars.sigIdxPostCut, zeros(1,this.vars.sigLen), 'g'); hold on
                    this.vars.figs.hBkg      = stem(this.vars.figs.axPhiPost, this.vars.bkgIdxPostCut, zeros(1,this.vars.bkgLen), 'r');
                    hold off
                    title(this.vars.figs.axPhiPost, "Post Cutting Phi")
                    
                    % Statistics Relations
                    this.vars.figs.axStats      = subplot(3,3,6);
                    this.vars.figs.hSNRStats    = plot(this.vars.figs.axStats, 1:this.vars.reps, zeros(1,this.vars.reps)); hold on
                    this.vars.figs.hMaxStats    = plot(this.vars.figs.axStats, 1:this.vars.reps, zeros(1,this.vars.reps));
                    this.vars.figs.hBkgStdStats = plot(this.vars.figs.axStats, 1:this.vars.reps, zeros(1,this.vars.reps));
                    legend (this.vars.figs.axStats, "SNR", "Sig", "Bkg std");
                    hold off
                    title(this.vars.figs.axStats, "Normalized SNR Components")
                    
                    % SNR Loop Trend
                    this.vars.figs.axTrend    = subplot(3,3,3);
                    this.vars.figs.hSNR       = plot(this.vars.figs.axTrend, 1:this.vars.reps, zeros(1,this.vars.reps));hold on
                    this.vars.figs.hSNRSubset = plot(this.vars.figs.axTrend, this.vars.subSetAx, zeros(1,this.vars.subSetLen));
                    this.vars.figs.hSNRTotal  = plot(this.vars.figs.axTrend, 1:this.vars.reps, zeros(1,this.vars.reps)); 
                    legend(this.vars.figs.axTrend, "SNR", "SubSet SNR", "Avg SNR")
                    hold off
                    title(this.vars.figs.axTrend, "SNR Trends")
                    
                    % Measurement SNR- Log Scaled
                    this.vars.figs.axErrLog = subplot(3,3,4);
                    yyaxis(this.vars.figs.axErrLog, 'left')
                    this.vars.figs.hErrLog  = errorbar(this.vars.figs.axErrLog, this.vars.pwr.refFactor, zeros(1, this.vars.measLen), zeros(1, this.vars.measLen)); hold on
                    this.vars.figs.hErrNLog = errorbar(this.vars.figs.axErrLog, this.vars.pwr.refFactor, zeros(1, this.vars.measLen), zeros(1, this.vars.measLen));
                    yyaxis(this.vars.figs.axErrLog, 'right')
                    this.vars.figs.hErrGLog = errorbar(this.vars.figs.axErrLog, this.vars.pwr.refFactor, zeros(1, this.vars.measLen), zeros(1, this.vars.measLen));
                    hold off
                    legend(this.vars.figs.axErrLog, "SNR", "Norm. SNR", 'SNR Gain')
                    set(this.vars.figs.axErrLog, 'XScale', 'log')
                    title(this.vars.figs.axErrLog, "Meas. SNR trend and Gain - Log")
                    
                    % Measurement SNR 
                    this.vars.figs.axErr = subplot(3,3,7);
                    yyaxis(this.vars.figs.axErr, 'left')
                    this.vars.figs.hErr  = errorbar(this.vars.figs.axErr, this.vars.errBarAx, zeros(1, this.vars.measLen), zeros(1, this.vars.measLen)); hold on
                    this.vars.figs.hErrN = errorbar(this.vars.figs.axErr, this.vars.errBarAx, zeros(1, this.vars.measLen), zeros(1, this.vars.measLen));
                    yyaxis(this.vars.figs.axErr, 'right')
                    this.vars.figs.hErrG = errorbar(this.vars.figs.axErr, this.vars.errBarAx, zeros(1, this.vars.measLen), zeros(1, this.vars.measLen));
                    hold off
                    legend(this.vars.figs.axErr, "SNR", "Norm. SNR", 'SNR Gain')
                    title(this.vars.figs.axErr, "Meas. SNR trend and Gain")
                    
                    % RawStd
                    this.vars.figs.axRStd   = subplot(3,3,9);
                    this.vars.figs.hRawStd = plot(this.vars.figs.axRStd, 1:this.vars.reps, zeros(1,this.vars.reps));
                    title(this.vars.figs.axRStd, "Raw Data Std")
                    
                    % FFT
                    this.vars.figs.axFFT = subplot(3,3,5);
                    this.vars.figs.hFftSig = plot(this.vars.figs.axFFT, this.vars.fBar*1e-6, zeros(1,this.vars.fftLen)); hold on
                    this.vars.figs.hFftBkg = plot(this.vars.figs.axFFT, this.vars.fBar*1e-6, zeros(1,this.vars.fftLen)); hold off
                    title(this.vars.figs.axFFT, "Sig and Bkg FFT")
                else
                    % Current Raw Signal
                    set(this.vars.figs.hNet, 'XData', this.vars.keepIdxs, 'YData', zeros(1, this.vars.phiLenPostCut));
                    if ~isempty(this.vars.cutIdxs)
                        if ~isfield(this.vars.figs, 'hCut') || ~isgraphics(this.vars.figs.hCut)
                            hold(this.vars.figs.axRaw, 'on');
                            this.vars.figs.hCut = stem(this.vars.figs.axRaw, this.vars.cutIdxs,  zeros(1, this.vars.cutLen), 'r'); hold on
                            hold(this.vars.figs.axRaw, 'off');
                        else
                            set(this.vars.figs.hCut, 'XData', this.vars.cutIdxs,  'YData', zeros(1, this.vars.cutLen));
                        end
                    else
                         this.vars.figs.hCut = 0;
                    end
                    
                    % Current Post-Cutting Signal
                    set(this.vars.figs.hSig, 'XData', this.vars.sigIdxPostCut, 'YData', zeros(1,this.vars.sigLen)); 
                    set(this.vars.figs.hBkg, 'XData', this.vars.bkgIdxPostCut, 'YData', zeros(1,this.vars.bkgLen));

                    % Statistics Relations
                    set(this.vars.figs.hSNRStats,    'XData', 1:this.vars.reps, 'YData', zeros(1,this.vars.reps));
                    set(this.vars.figs.hMaxStats,    'XData', 1:this.vars.reps, 'YData', zeros(1,this.vars.reps));
                    set(this.vars.figs.hBkgStdStats, 'XData', 1:this.vars.reps, 'YData', zeros(1,this.vars.reps)); 

                    % SNR Loop Trend
                    set(this.vars.figs.hSNR,       'XData', 1:this.vars.reps,   'YData', zeros(1,this.vars.reps));
                    set(this.vars.figs.hSNRSubset, 'XData', this.vars.subSetAx, 'YData', zeros(1,this.vars.subSetLen));
                    set(this.vars.figs.hSNRTotal,  'XData', 1:this.vars.reps,   'YData', zeros(1,this.vars.reps));

                    % Measurement SNR - Log
                    set(this.vars.figs.hErrLog,  'XData', this.vars.pwr.refFactor, 'YData', zeros(1, this.vars.measLen), 'YNegativeDelta', zeros(1, this.vars.measLen), 'YPositiveDelta', zeros(1, this.vars.measLen));
                    set(this.vars.figs.hErrNLog, 'XData', this.vars.pwr.refFactor, 'YData', zeros(1, this.vars.measLen), 'YNegativeDelta', zeros(1, this.vars.measLen), 'YPositiveDelta', zeros(1, this.vars.measLen));
                    set(this.vars.figs.hErrGLog, 'XData', this.vars.pwr.refFactor, 'YData', zeros(1, this.vars.measLen), 'YNegativeDelta', zeros(1, this.vars.measLen), 'YPositiveDelta', zeros(1, this.vars.measLen));
                    
                    % Measurement SNR
                    set(this.vars.figs.hErr,  'XData', this.vars.errBarAx, 'YData', zeros(1, this.vars.measLen), 'YNegativeDelta', zeros(1, this.vars.measLen), 'YPositiveDelta', zeros(1, this.vars.measLen));
                    set(this.vars.figs.hErrN, 'XData', this.vars.errBarAx, 'YData', zeros(1, this.vars.measLen), 'YNegativeDelta', zeros(1, this.vars.measLen), 'YPositiveDelta', zeros(1, this.vars.measLen));
                    set(this.vars.figs.hErrG, 'XData', this.vars.errBarAx, 'YData', zeros(1, this.vars.measLen), 'YNegativeDelta', zeros(1, this.vars.measLen), 'YPositiveDelta', zeros(1, this.vars.measLen));
                    
                    % RawStd
                    set(this.vars.figs.hRawStd, 'XData', 1:this.vars.reps, 'YData', zeros(1,this.vars.reps));
                    
                    
                    % FFT
                    set(this.vars.figs.hFftSig, 'XData', this.vars.fBar*1e-6, 'YData', zeros(1,this.vars.fftLen));
                    set(this.vars.figs.hFftBkg, 'XData', this.vars.fBar*1e-6, 'YData', zeros(1,this.vars.fftLen));
                    
                end

                xlim(this.vars.figs.axErr, this.vars.figs.lims)
                
            catch
                fprintf("Error initing figures\n");
                return 
            end
        end
        
        function initFigsDual(this)
            try 
                if ~isgraphics(this.vars.figs.dual.hFig(1))      || ~isgraphics(this.vars.figs.dual.hFig(2))      || ~isgraphics(this.vars.figs.dual.hFig(3)) ||...
                   ~isgraphics(this.vars.figs.dual.axRawPD)      || ~isgraphics(this.vars.figs.dual.axRawPMT)     ||...
                   ~isgraphics(this.vars.figs.dual.axPostCutPD)  || ~isgraphics(this.vars.figs.dual.axPostCutPMT) ||...
                   ~isgraphics(this.vars.figs.dual.axFFTPD)      || ~isgraphics(this.vars.figs.dual.axFFTPMT)     ||...
                   ~isgraphics(this.vars.figs.dual.axStatsPD)    || ~isgraphics(this.vars.figs.dual.axStatsPMT)   ||...
                   ~isgraphics(this.vars.figs.dual.axTrendPD)    || ~isgraphics(this.vars.figs.dual.axTrendPMT)   ||...
                   ~isgraphics(this.vars.figs.dual.axStd)        || ~isgraphics(this.vars.figs.dual.axGain)       ||...
                   ~isgraphics(this.vars.figs.dual.axGainLog)
                
                    if ~isgraphics(handle(this.vars.figs.dual.hFig))
                        desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
                        myGroup = desktop.addGroup('myGroup');
                        desktop.setGroupDocked('myGroup', 0);
                        myDim   = java.awt.Dimension(4, 2);                         
                        % 4 columns, 2 rows
                        % 1: Maximized, 2: Tiled, 3: Floating
                        desktop.setDocumentArrangement('myGroup', 2, myDim)
                        bakWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
                        for i = 1:3
                           this.vars.figs.dual.hFig(i) = figure('WindowStyle', 'docked');
                           if i == 1
                               this.vars.figs.dual.axRawPD      = subplot(2,3,1);
                               this.vars.figs.dual.axRawPMT     = subplot(2,3,4);
                               this.vars.figs.dual.axPostCutPD  = subplot(2,3,2);
                               this.vars.figs.dual.axPostCutPMT = subplot(2,3,5);
                               this.vars.figs.dual.axFFTPD      = subplot(2,3,3);
                               this.vars.figs.dual.axFFTPMT     = subplot(2,3,6);
                           elseif i == 2
                               this.vars.figs.dual.axStatsPD  = subplot(2,3,1);
                               this.vars.figs.dual.axStatsPMT = subplot(2,3,4);
                               this.vars.figs.dual.axTrendPD  = subplot(2,3,2);
                               this.vars.figs.dual.axTrendPMT = subplot(2,3,5);
                               this.vars.figs.dual.axStd      = subplot(1,3,3);
                           else                               
                               this.vars.figs.dual.axGain     = subplot(1,2,1);
                               this.vars.figs.dual.axGainLog  = subplot(1,2,2);
                           end
                           drawnow; pause(0.02);  % Magic, reduces rendering errors
                           set(get(handle(this.vars.figs.dual.hFig(i)), 'javaframe'), 'GroupName', 'myGroup');
                        end
                        warning(bakWarn);
                    end

                    % Current Raw Signal - PD
                    hold(this.vars.figs.dual.axRawPD, 'on');
                    this.vars.figs.dual.hNetPD      = stem(this.vars.figs.dual.axRawPD, this.vars.keepIdxs, zeros(1, this.vars.phiLenPostCut)); 
                    if ~isempty(this.vars.cutIdxs)
                        this.vars.figs.dual.hCutPD = stem(this.vars.figs.dual.axRawPD, this.vars.cutIdxs,  zeros(1, this.vars.cutLen), 'r');
                    else
                        this.vars.figs.dual.hCutPD = 0;
                    end
                    hold(this.vars.figs.dual.axRawPD, 'off');
                    title(this.vars.figs.dual.axRawPD, "Raw Phi - PD")
                    
                    % Current Raw Signal - PMT
                    hold(this.vars.figs.dual.axRawPMT, 'on');
                    this.vars.figs.dual.hNetPMT      = stem(this.vars.figs.dual.axRawPMT, this.vars.keepIdxs, zeros(1, this.vars.phiLenPostCut));
                    if ~isempty(this.vars.cutIdxs)
                        this.vars.figs.dual.hCutPMT = stem(this.vars.figs.dual.axRawPMT, this.vars.cutIdxs,  zeros(1, this.vars.cutLen), 'r');
                    else
                        this.vars.figs.dual.hCutPMT = 0;
                    end
                    hold(this.vars.figs.dual.axRawPMT, 'off');
                    title(this.vars.figs.dual.axRawPMT, "Raw Phi - PMT")

                    %Current Post-Cutting Signal - PD
                    hold(this.vars.figs.dual.axPostCutPD, 'on');
                    this.vars.figs.dual.hSigPD = stem(this.vars.figs.dual.axPostCutPD, this.vars.sigIdxPostCut, zeros(1,this.vars.sigLen), 'g');
                    this.vars.figs.dual.hBkgPD = stem(this.vars.figs.dual.axPostCutPD, this.vars.bkgIdxPostCut, zeros(1,this.vars.bkgLen), 'r');
                    hold(this.vars.figs.dual.axPostCutPD, 'off');
                    title(this.vars.figs.dual.axPostCutPD, "Post Cutting Phi - PD")
                    
                    %Current Post-Cutting Signal - PMT
                    hold(this.vars.figs.dual.axPostCutPMT, 'on');
                    this.vars.figs.dual.hSigPMT      = stem(this.vars.figs.dual.axPostCutPMT, this.vars.sigIdxPostCut, zeros(1,this.vars.sigLen), 'g'); 
                    this.vars.figs.dual.hBkgPMT      = stem(this.vars.figs.dual.axPostCutPMT, this.vars.bkgIdxPostCut, zeros(1,this.vars.bkgLen), 'r');
                    hold(this.vars.figs.dual.axPostCutPMT, 'off');
                    title(this.vars.figs.dual.axPostCutPMT, "Post Cutting Phi - PMT")
                    
                    % FFT - PD
                    hold(this.vars.figs.dual.axFFTPD, 'on');
                    this.vars.figs.dual.hFftSigPD = plot(this.vars.figs.dual.axFFTPD, this.vars.fBar*1e-6, zeros(1,this.vars.fftLen)); 
                    this.vars.figs.dual.hFftBkgPD = plot(this.vars.figs.dual.axFFTPD, this.vars.fBar*1e-6, zeros(1,this.vars.fftLen)); 
                    hold(this.vars.figs.dual.axFFTPD, 'off');
                    title(this.vars.figs.dual.axFFTPD, "Sig and Bkg FFT - PD")
                    
                    % FFT - PMT
                    hold(this.vars.figs.dual.axFFTPMT, 'on');
                    this.vars.figs.dual.hFftSigPMT = plot(this.vars.figs.dual.axFFTPMT, this.vars.fBar*1e-6, zeros(1,this.vars.fftLen)); 
                    this.vars.figs.dual.hFftBkgPMT = plot(this.vars.figs.dual.axFFTPMT, this.vars.fBar*1e-6, zeros(1,this.vars.fftLen)); 
                    hold(this.vars.figs.dual.axFFTPMT, 'off');
                    title(this.vars.figs.dual.axFFTPMT, "Sig and Bkg FFT - PMT")
                    
                    % Statistics Relations - PD
                    hold(this.vars.figs.dual.axStatsPD, 'on');
                    this.vars.figs.dual.hSNRStatsPD    = plot(this.vars.figs.dual.axStatsPD, 1:this.vars.reps, zeros(1,this.vars.reps)); 
                    this.vars.figs.dual.hMaxStatsPD    = plot(this.vars.figs.dual.axStatsPD, 1:this.vars.reps, zeros(1,this.vars.reps));
                    this.vars.figs.dual.hBkgStdStatsPD = plot(this.vars.figs.dual.axStatsPD, 1:this.vars.reps, zeros(1,this.vars.reps));
                    hold(this.vars.figs.dual.axStatsPD, 'off');
                    legend (this.vars.figs.dual.axStatsPD, "SNR", "Sig", "Bkg std");
                    title(this.vars.figs.dual.axStatsPD, "Norm. SNR Components - PD")
                    
                    % Statistics Relations - PMT
                    hold(this.vars.figs.dual.axStatsPMT, 'on');
                    this.vars.figs.dual.hSNRStatsPMT    = plot(this.vars.figs.dual.axStatsPMT, 1:this.vars.reps, zeros(1,this.vars.reps)); 
                    this.vars.figs.dual.hMaxStatsPMT    = plot(this.vars.figs.dual.axStatsPMT, 1:this.vars.reps, zeros(1,this.vars.reps));
                    this.vars.figs.dual.hBkgStdStatsPMT = plot(this.vars.figs.dual.axStatsPMT, 1:this.vars.reps, zeros(1,this.vars.reps));
                    hold(this.vars.figs.dual.axStatsPMT, 'off');
                    legend (this.vars.figs.dual.axStatsPMT, "SNR", "Sig", "Bkg std");
                    title(this.vars.figs.dual.axStatsPMT, "Norm. SNR Components - PMT")
                    
                    % SNR Loop Trend - PD
                    hold(this.vars.figs.dual.axTrendPD, 'on');
                    this.vars.figs.dual.hSNRPD       = plot(this.vars.figs.dual.axTrendPD, 1:this.vars.reps, zeros(1,this.vars.reps));
                    this.vars.figs.dual.hSNRSubsetPD = plot(this.vars.figs.dual.axTrendPD, this.vars.subSetAx, zeros(1,this.vars.subSetLen));
                    this.vars.figs.dual.hSNRTotalPD  = plot(this.vars.figs.dual.axTrendPD, 1:this.vars.reps, zeros(1,this.vars.reps)); 
                    hold(this.vars.figs.dual.axTrendPD, 'off');
                    legend(this.vars.figs.dual.axTrendPD, "SNR", "SubSet SNR", "Avg SNR")
                    title(this.vars.figs.dual.axTrendPD, "SNR Trends - PD")
                    
                    % SNR Loop Trend - PMT
                    hold(this.vars.figs.dual.axTrendPMT, 'on');
                    this.vars.figs.dual.hSNRPMT       = plot(this.vars.figs.dual.axTrendPMT, 1:this.vars.reps, zeros(1,this.vars.reps));
                    this.vars.figs.dual.hSNRSubsetPMT = plot(this.vars.figs.dual.axTrendPMT, this.vars.subSetAx, zeros(1,this.vars.subSetLen));
                    this.vars.figs.dual.hSNRTotalPMT  = plot(this.vars.figs.dual.axTrendPMT, 1:this.vars.reps, zeros(1,this.vars.reps)); 
                    hold(this.vars.figs.dual.axTrendPMT, 'off');
                    legend(this.vars.figs.dual.axTrendPMT, "SNR", "SubSet SNR", "Avg SNR")
                    title(this.vars.figs.dual.axTrendPMT, "SNR Trends - PMT")
                    
                    % RawStd
                    hold(this.vars.figs.dual.axStd, 'on');
                    yyaxis(this.vars.figs.dual.axStd, 'left')
                    this.vars.figs.dual.hRawStdPD  = plot(this.vars.figs.dual.axStd, 1:this.vars.reps, zeros(1,this.vars.reps)); 
                    yyaxis(this.vars.figs.dual.axStd, 'right')
                    this.vars.figs.dual.hRawStdPMT = plot(this.vars.figs.dual.axStd, 1:this.vars.reps, zeros(1,this.vars.reps));
                    hold(this.vars.figs.dual.axStd, 'off');
                    legend(this.vars.figs.dual.axStd, "PD", "PMT")
                    title(this.vars.figs.dual.axStd, "Raw Data Std")

                    % Measurement SNR Gain
                    hold(this.vars.figs.dual.axGain, 'on');
                    yyaxis(this.vars.figs.dual.axGain, 'left')
                    this.vars.figs.dual.hSNRMeasPD  = errorbar(this.vars.figs.dual.axGain, this.vars.errBarAx, zeros(1, this.vars.measLen), zeros(1, this.vars.measLen)); 
                    this.vars.figs.dual.hSNRMeasPMT = errorbar(this.vars.figs.dual.axGain, this.vars.errBarAx, zeros(1, this.vars.measLen), zeros(1, this.vars.measLen));
                    yyaxis(this.vars.figs.dual.axGain, 'right')
                    this.vars.figs.dual.hSNRGain = errorbar(this.vars.figs.dual.axGain, this.vars.errBarAx, zeros(1, this.vars.measLen), zeros(1, this.vars.measLen));
                    hold(this.vars.figs.dual.axGain, 'off');
                    legend(this.vars.figs.dual.axGain, "SNR-PD", "SNR-PMT", 'SNR Gain')
                    title(this.vars.figs.dual.axGain, "Meas. SNR trend and Gain")
                    
                    % Measurement SNR Gain - Log Scaled
                    hold(this.vars.figs.dual.axGainLog, 'on');
                    yyaxis(this.vars.figs.dual.axGainLog, 'left')
                    this.vars.figs.dual.hSNRMeasLogPD  = errorbar(this.vars.figs.dual.axGainLog, this.vars.pwr.refFactor, zeros(1, this.vars.measLen), zeros(1, this.vars.measLen)); 
                    this.vars.figs.dual.hSNRMeasLogPMT = errorbar(this.vars.figs.dual.axGainLog, this.vars.pwr.refFactor, zeros(1, this.vars.measLen), zeros(1, this.vars.measLen));
                    yyaxis(this.vars.figs.dual.axGainLog, 'right')
                    this.vars.figs.dual.hSNRGainLog = errorbar(this.vars.figs.dual.axGainLog, this.vars.pwr.refFactor, zeros(1, this.vars.measLen), zeros(1, this.vars.measLen));
                    hold(this.vars.figs.dual.axGainLog, 'off');
                    legend(this.vars.figs.dual.axGainLog, "SNR", "SNR-PMT", "SNR Gain")
                    set(this.vars.figs.dual.axGainLog, 'XScale', 'log')
                    title(this.vars.figs.dual.axGainLog, "Meas. SNR trend and Gain - Log") 
                else
                    % Current Raw Signal - PD
                    set(this.vars.figs.dual.hNetPD, 'XData', this.vars.keepIdxs, 'YData', zeros(1, this.vars.phiLenPostCut));
                    if ~isempty(this.vars.cutIdxs)
                        if ~isfield(this.vars.figs.dual, 'hCutPD') || ~isgraphics(this.vars.figs.dual.hCutPD)
                            hold(this.vars.figs.dual.axRawPD, 'on');
                            this.vars.figs.dual.hCutPD = stem(this.vars.figs.dual.axRawPD, this.vars.cutIdxs,  zeros(1, this.vars.cutLen), 'r'); hold on
                            hold(this.vars.figs.dual.axRawPD, 'off');
                        else
                            set(this.vars.figs.dual.hCutPD, 'XData', this.vars.cutIdxs,  'YData', zeros(1, this.vars.cutLen));
                        end
                    else
                         this.vars.figs.dual.hCutPD = 0;
                    end
                    
                    % Current Raw Signal - PMT
                    set(this.vars.figs.dual.hNetPMT, 'XData', this.vars.keepIdxs, 'YData', zeros(1, this.vars.phiLenPostCut));
                    if ~isempty(this.vars.cutIdxs)
                        if ~isfield(this.vars.figs.dual, 'hCutPD') || ~isgraphics(this.vars.figs.dual.hCutPMT)
                            hold(this.vars.figs.dual.axRawPMT, 'on');
                            this.vars.figs.dual.hCutPMT = stem(this.vars.figs.dual.axRawPMT, this.vars.cutIdxs,  zeros(1, this.vars.cutLen), 'r'); hold on
                            hold(this.vars.figs.dual.axRawPMT, 'off');
                        else
                            set(this.vars.figs.dual.hCutPMT, 'XData', this.vars.cutIdxs,  'YData', zeros(1, this.vars.cutLen));
                        end
                    else
                         this.vars.figs.dual.hCutPMT = 0;
                    end

                    % Current Post-Cutting Signal
                    set(this.vars.figs.dual.hSigPD,  'XData', this.vars.sigIdxPostCut, 'YData', zeros(1,this.vars.sigLen));
                    set(this.vars.figs.dual.hBkgPD,  'XData', this.vars.bkgIdxPostCut, 'YData', zeros(1,this.vars.bkgLen));
                    set(this.vars.figs.dual.hSigPMT, 'XData', this.vars.sigIdxPostCut, 'YData', zeros(1,this.vars.sigLen));
                    set(this.vars.figs.dual.hBkgPMT, 'XData', this.vars.bkgIdxPostCut, 'YData', zeros(1,this.vars.bkgLen));

                    % FFT
                    set(this.vars.figs.dual.hFftSigPD,  'XData', this.vars.fBar*1e-6, 'YData', zeros(1,this.vars.fftLen));
                    set(this.vars.figs.dual.hFftBkgPD,  'XData', this.vars.fBar*1e-6, 'YData', zeros(1,this.vars.fftLen));
                    set(this.vars.figs.dual.hFftSigPMT, 'XData', this.vars.fBar*1e-6, 'YData', zeros(1,this.vars.fftLen));
                    set(this.vars.figs.dual.hFftBkgPMT, 'XData', this.vars.fBar*1e-6, 'YData', zeros(1,this.vars.fftLen));

                    % Statistics Relations
                    set(this.vars.figs.dual.hSNRStatsPD,     'XData', 1:this.vars.reps, 'YData', zeros(1,this.vars.reps));
                    set(this.vars.figs.dual.hMaxStatsPD,     'XData', 1:this.vars.reps, 'YData', zeros(1,this.vars.reps));
                    set(this.vars.figs.dual.hBkgStdStatsPD,  'XData', 1:this.vars.reps, 'YData', zeros(1,this.vars.reps));
                    set(this.vars.figs.dual.hSNRStatsPMT,    'XData', 1:this.vars.reps, 'YData', zeros(1,this.vars.reps));
                    set(this.vars.figs.dual.hMaxStatsPMT,    'XData', 1:this.vars.reps, 'YData', zeros(1,this.vars.reps));
                    set(this.vars.figs.dual.hBkgStdStatsPMT, 'XData', 1:this.vars.reps, 'YData', zeros(1,this.vars.reps));

                    % SNR Loop Trend
                    set(this.vars.figs.dual.hSNRPD,        'XData', 1:this.vars.reps,   'YData', zeros(1,this.vars.reps));
                    set(this.vars.figs.dual.hSNRSubsetPD,  'XData', this.vars.subSetAx, 'YData', zeros(1,this.vars.subSetLen));
                    set(this.vars.figs.dual.hSNRTotalPD,   'XData', 1:this.vars.reps,   'YData', zeros(1,this.vars.reps));
                    set(this.vars.figs.dual.hSNRPMT,       'XData', 1:this.vars.reps,   'YData', zeros(1,this.vars.reps));
                    set(this.vars.figs.dual.hSNRSubsetPMT, 'XData', this.vars.subSetAx, 'YData', zeros(1,this.vars.subSetLen));
                    set(this.vars.figs.dual.hSNRTotalPMT,  'XData', 1:this.vars.reps,   'YData', zeros(1,this.vars.reps));

                    % RawStd
                    set(this.vars.figs.dual.hRawStdPD,  'XData', 1:this.vars.reps, 'YData', zeros(1,this.vars.reps));
                    set(this.vars.figs.dual.hRawStdPMT, 'XData', 1:this.vars.reps, 'YData', zeros(1,this.vars.reps));

                    % Measurement SNR
                    set(this.vars.figs.dual.hSNRMeasPD,  'XData', this.vars.errBarAx, 'YData', zeros(1, this.vars.measLen), 'YNegativeDelta', zeros(1, this.vars.measLen), 'YPositiveDelta', zeros(1, this.vars.measLen));
                    set(this.vars.figs.dual.hSNRMeasPMT, 'XData', this.vars.errBarAx, 'YData', zeros(1, this.vars.measLen), 'YNegativeDelta', zeros(1, this.vars.measLen), 'YPositiveDelta', zeros(1, this.vars.measLen));
                    set(this.vars.figs.dual.hSNRGain,    'XData', this.vars.errBarAx, 'YData', zeros(1, this.vars.measLen), 'YNegativeDelta', zeros(1, this.vars.measLen), 'YPositiveDelta', zeros(1, this.vars.measLen));

                    % Measurement SNR - Log
                    set(this.vars.figs.dual.hSNRMeasLogPD,  'XData', this.vars.pwr.refFactor, 'YData', zeros(1, this.vars.measLen), 'YNegativeDelta', zeros(1, this.vars.measLen), 'YPositiveDelta', zeros(1, this.vars.measLen));
                    set(this.vars.figs.dual.hSNRMeasLogPMT, 'XData', this.vars.pwr.refFactor, 'YData', zeros(1, this.vars.measLen), 'YNegativeDelta', zeros(1, this.vars.measLen), 'YPositiveDelta', zeros(1, this.vars.measLen));
                    set(this.vars.figs.dual.hSNRGainLog,    'XData', this.vars.pwr.refFactor, 'YData', zeros(1, this.vars.measLen), 'YNegativeDelta', zeros(1, this.vars.measLen), 'YPositiveDelta', zeros(1, this.vars.measLen));
                end

                xlim(this.vars.figs.dual.axGain,    this.vars.figs.lims)
%                 xlim(this.vars.figs.dual.axGainLog, this.vars.figs.lims)
            catch
                fprintf("Error initing figures\n");
                return 
            end
        end
        
        function plotResultsSingle(this)
            j = this.vars.idxs(1);
            i = this.vars.idxs(2);
            
            % Current Raw Signal
            set(this.vars.figs.hNet, 'YData', this.data.phiNet(j,i,:));
            if ~isempty(this.vars.cutIdxs)
                set(this.vars.figs.hCut, 'YData', this.data.phiCut(j,i,:));
            end
            
            % Current Post-Cutting Signal
            set(this.vars.figs.hSig, 'YData', this.data.sig(j,i,:));
            set(this.vars.figs.hBkg, 'YData', this.data.bkg(j,i,:));

            % Statistics Relations
            set(this.vars.figs.hSNRStats, 'YData',this.data.normSNR(j,:));
            set(this.vars.figs.hMaxStats, 'YData', this.data.normMax(j,:));
            set(this.vars.figs.hBkgStdStats, 'YData', this.data.normStd(j,:));
            
            % SNR Loop Trend
            set(this.vars.figs.hSNR,     'YData', this.data.SNR(j,:));
            set(this.vars.figs.hSNRTotal, 'YData', this.data.avgSNR(j,:));
            set(this.vars.figs.hSNRSubset, 'YData', this.data.SNRsubset(j,:));
            
            % Measurement SNR - Log
            set(this.vars.figs.hErrLog,   'YData',  this.data.SNRAvgMeas, 'YNegativeDelta', this.data.SNRstdMeas, 'YPositiveDelta', this.data.SNRstdMeas);
            set(this.vars.figs.hErrNLog,  'YData',  this.data.SNRAvgNorm, 'YNegativeDelta', this.data.SNRStdNorm, 'YPositiveDelta', this.data.SNRStdNorm);
            set(this.vars.figs.hErrGLog,  'YData',  this.data.SNRAvgGain, 'YNegativeDelta', this.data.SNRStdGain, 'YPositiveDelta', this.data.SNRStdGain);
            
            % Measurement SNR
            set(this.vars.figs.hErr,   'YData',  this.data.SNRAvgMeas, 'YNegativeDelta', this.data.SNRstdMeas, 'YPositiveDelta', this.data.SNRstdMeas);
            set(this.vars.figs.hErrN,  'YData',  this.data.SNRAvgNorm, 'YNegativeDelta', this.data.SNRStdNorm, 'YPositiveDelta', this.data.SNRStdNorm);
            set(this.vars.figs.hErrG,  'YData',  this.data.SNRAvgGain, 'YNegativeDelta', this.data.SNRStdGain, 'YPositiveDelta', this.data.SNRStdGain);
            
            % RawStd
            set(this.vars.figs.hRawStd, 'YData', this.data.rawStd(j, :));
            minstd = min(this.data.rawStd(j, 1:i));
            maxstd = max(this.data.rawStd(j, 1:i));
            avgstd = mean(this.data.rawStd(j, 1:i));
            
            ylim(this.vars.figs.axRStd, [minstd-0.1*avgstd, maxstd+0.1*avgstd]);
            xlim(this.vars.figs.axRStd, [1,this.vars.reps]);
            
            %FFT
            set(this.vars.figs.hFftSig, 'YData', this.data.fftResSig(j,i,:)); 
            set(this.vars.figs.hFftBkg, 'YData', this.data.fftResBkg(j,i,:));
            
            drawnow;  
        end
        
        function plotResultsDual(this)
            j = this.vars.idxs(1);
            i = this.vars.idxs(2);
            
            % Current Raw Signal
            set(this.vars.figs.dual.hNetPD, 'YData', this.data.phiNet(j,i,:,1));
            if ~isempty(this.vars.cutIdxs)
                set(this.vars.figs.dual.hCutPD, 'YData', this.data.phiCut(j,i,:,1));
            end
            set(this.vars.figs.dual.hNetPMT, 'YData', this.data.phiNet(j,i,:,2));
            if ~isempty(this.vars.cutIdxs)
                set(this.vars.figs.dual.hCutPMT, 'YData', this.data.phiCut(j,i,:,2));
            end
            
            % Current Post-Cutting Signal
            set(this.vars.figs.dual.hSigPD,  'YData', this.data.sig(j,i,:,1));
            set(this.vars.figs.dual.hBkgPD,  'YData', this.data.bkg(j,i,:,1));
            set(this.vars.figs.dual.hSigPMT, 'YData', this.data.sig(j,i,:,2));
            set(this.vars.figs.dual.hBkgPMT, 'YData', this.data.bkg(j,i,:,2));

            %FFT
            set(this.vars.figs.dual.hFftSigPD,  'YData', this.data.fftResSig(j,i,:,1)); 
            set(this.vars.figs.dual.hFftBkgPD,  'YData', this.data.fftResBkg(j,i,:,1));
            set(this.vars.figs.dual.hFftSigPMT, 'YData', this.data.fftResSig(j,i,:,2)); 
            set(this.vars.figs.dual.hFftBkgPMT, 'YData', this.data.fftResBkg(j,i,:,2));
            
            % Statistics Relations
            set(this.vars.figs.dual.hSNRStatsPD,     'YData', this.data.normSNR(j,:,1));
            set(this.vars.figs.dual.hMaxStatsPD,     'YData', this.data.normMax(j,:,1));
            set(this.vars.figs.dual.hBkgStdStatsPD,  'YData', this.data.normStd(j,:,1));
            set(this.vars.figs.dual.hSNRStatsPMT,    'YData', this.data.normSNR(j,:,2));
            set(this.vars.figs.dual.hMaxStatsPMT,    'YData', this.data.normMax(j,:,2));
            set(this.vars.figs.dual.hBkgStdStatsPMT, 'YData', this.data.normStd(j,:,2));
            
            % SNR Loop Trend
            set(this.vars.figs.dual.hSNRPD,        'YData', this.data.SNR(j,:,1));
            set(this.vars.figs.dual.hSNRTotalPD,   'YData', this.data.avgSNR(j,:,1));
            set(this.vars.figs.dual.hSNRSubsetPD,  'YData', this.data.SNRsubset(j,:,1));
            set(this.vars.figs.dual.hSNRPMT,       'YData', this.data.SNR(j,:,2));
            set(this.vars.figs.dual.hSNRTotalPMT,  'YData', this.data.avgSNR(j,:,2));
            set(this.vars.figs.dual.hSNRSubsetPMT, 'YData', this.data.SNRsubset(j,:,2));
            
            % RawStd
            set(this.vars.figs.dual.hRawStdPD, 'YData', this.data.rawStd(j, :, 1));
            minStdPD = min(this.data.rawStd (j, 1:i, 1));
            maxStdPD = max(this.data.rawStd (j, 1:i, 1));
            avgStdPD = mean(this.data.rawStd(j, 1:i, 1));
            minLimPD = minStdPD-0.1*avgStdPD;
            maxLimPD = maxStdPD+0.1*avgStdPD;
            
            set(this.vars.figs.dual.hRawStdPMT, 'YData', this.data.rawStd(j, :, 2));
            minStdPMT = min(this.data.rawStd (j, 1:i, 2));
            maxStdPMT = max(this.data.rawStd (j, 1:i, 2));
            avgStdPMT = mean(this.data.rawStd(j, 1:i, 2));
            
            minLimPMT = minStdPMT-0.1*avgStdPMT;
            maxLimPMT = maxStdPMT+0.1*avgStdPMT;
            
            lims = [min(minLimPD, minLimPMT), max(maxLimPD, maxLimPMT)];
            ylim(this.vars.figs.dual.axStd, lims);
            xlim(this.vars.figs.dual.axStd, [1,this.vars.reps]);
            
            % Measurement SNR
            set(this.vars.figs.dual.hSNRMeasPD,  'YData',  this.data.SNRAvgMeas(1,:,1), 'YNegativeDelta', this.data.SNRstdMeas(1,:,1), 'YPositiveDelta', this.data.SNRstdMeas(1,:,1));
            set(this.vars.figs.dual.hSNRMeasPMT, 'YData',  this.data.SNRAvgMeas(1,:,2), 'YNegativeDelta', this.data.SNRstdMeas(1,:,2), 'YPositiveDelta', this.data.SNRstdMeas(1,:,2));
            set(this.vars.figs.dual.hSNRGain,    'YData',  this.data.SNRAvgGain, 'YNegativeDelta', this.data.SNRStdGain, 'YPositiveDelta', this.data.SNRStdGain);
            
            % Measurement SNR - Log
            set(this.vars.figs.dual.hSNRMeasLogPD,  'YData',  this.data.SNRAvgMeas(1,:,1), 'YNegativeDelta', this.data.SNRstdMeas(1,:,1), 'YPositiveDelta', this.data.SNRstdMeas(1,:,1));
            set(this.vars.figs.dual.hSNRMeasLogPMT, 'YData',  this.data.SNRAvgMeas(1,:,2), 'YNegativeDelta', this.data.SNRstdMeas(1,:,2), 'YPositiveDelta', this.data.SNRstdMeas(1,:,2));
            set(this.vars.figs.dual.hSNRGainLog,    'YData',  this.data.SNRAvgGain, 'YNegativeDelta', this.data.SNRStdGain, 'YPositiveDelta', this.data.SNRStdGain);

            drawnow;
        end
        
        function initDataStructs(this)
            % Init database variables
            if this.vars.uVarsAO.ao.exportData.meas && this.vars.saveFlag && (this.vars.saveN>=1)
                this.data.rawData = zeros(this.vars.saveN, this.vars.ao.measVars.algo.digitizer.samplesPerAcqAllCh);
            end
            
            ch = 1;
            if this.vars.dualOP
                ch = 2;
            end
            
            this.data.rawPhi = zeros(this.vars.measLen, this.vars.reps, this.vars.phiLen, ch);
            
            this.data.phiCut = zeros(this.vars.measLen, this.vars.reps, this.vars.cutLen, ch);
            this.data.phiNet = zeros(this.vars.measLen, this.vars.reps, this.vars.phiLenPostCut, ch);
            this.data.bkg    = zeros(this.vars.measLen, this.vars.reps, this.vars.bkgLen, ch);
            this.data.sig    = zeros(this.vars.measLen, this.vars.reps, this.vars.sigLen, ch);
            
            this.data.fftResSig    = zeros(this.vars.measLen, this.vars.reps, this.vars.fftLen, ch);
            this.data.fftResBkg    = zeros(this.vars.measLen, this.vars.reps, this.vars.fftLen, ch);
            this.data.fftResBkgAvg = zeros(this.vars.measLen, this.vars.reps, this.vars.fftLen, ch);
            this.data.fftResBkgStd = zeros(this.vars.measLen, this.vars.reps, this.vars.fftLen, ch);
            
            this.data.maxVal  = zeros(this.vars.measLen, this.vars.reps, ch);
            this.data.avgBkg  = zeros(this.vars.measLen, this.vars.reps, ch);
            this.data.stdBkg  = zeros(this.vars.measLen, this.vars.reps, ch);
            
            this.data.normMax = ones(this.vars.measLen, this.vars.reps, ch);
            this.data.normStd = ones(this.vars.measLen, this.vars.reps, ch);
            
            this.data.SNR    = zeros(this.vars.measLen, this.vars.reps, ch);
            this.data.avgSNR = zeros(this.vars.measLen, this.vars.reps, ch);
            this.data.stdSNR = zeros(this.vars.measLen, this.vars.reps, ch);

            this.data.timeNormSNR    = zeros(this.vars.measLen, this.vars.reps, ch);
            this.data.avgTimeNormSNR = zeros(this.vars.measLen, this.vars.reps, ch);
            this.data.stdTimeNormSNR = zeros(this.vars.measLen, this.vars.reps, ch);

            this.data.normSNR    = ones(this.vars.measLen, this.vars.reps, ch);
            this.data.SNRsubset  = zeros(this.vars.measLen, this.vars.subSetLen, ch);
            
            this.data.SNRAvgMeas  = zeros(1, this.vars.measLen, ch);
            this.data.SNRstdMeas  = zeros(1, this.vars.measLen, ch);
            this.data.SNRAvgNorm  = zeros(1, this.vars.measLen, ch);
            this.data.SNRStdNorm  = zeros(1, this.vars.measLen, ch);
            this.data.SNRAvgGain  = zeros(1, this.vars.measLen);
            this.data.SNRStdGain  = zeros(1, this.vars.measLen);
            
            this.data.SNRGain     = zeros(this.vars.measLen, this.vars.reps);
            
            this.data.rawStd = zeros(this.vars.measLen, this.vars.reps, ch);
            
            this.data.dc = [];

            if isfield(this.data, 'res')
                this.data = rmfield(this.data, 'res');
            end
            if isfield(this.data, 'resTotal')
            this.data = rmfield(this.data, 'resTotal');
            end
            if isfield(this.data, 'resMeas')
                this.data = rmfield(this.data, 'resMeas');
            end
        end
    
        function exp = measure(this)
            switch this.vars.measType
                case {'PD'; 'PMT'; 'Point'}
                    exp = this.measureSingle();
                case {'Dual'; 'Dual-Point'}
                    exp = this.measureDual();
                otherwise
                    exp = [];
            end
        end
            
        function exp = measureDual(this) 
            this.initDataStructs();
            this.initFigsDual();
            this.initFS();
            
            for j = 1:this.vars.measLen
                this.vars.idxs(1) = j;
                this.vars.measTime(j) = tic;
                switch this.vars.measType
                    case 'Dual'
                        this.stagesT.moveAbsAx('X', this.vars.thetaVec(j))
                    case 'Dual-Point'
                        if this.vars.rotate
                            this.stagesT.moveAbsAx('X', this.vars.thetaVec(j))
                        end
                end
                
                fprintf("Meas: %d\n", j);
                
                this.data.curVars = this.ao.getVars();
                for i=1:this.vars.reps   
                    this.vars.idxs(2) = i;
                    pause(this.vars.stabTime)
                    
                    this.data.res(i,:)       = this.ao.runAcoustoOptics();
                    for k=1:2
                        if  this.vars.uVarsAO.ao.exportData.meas
                            if i<= this.vars.saveN
                                this.data.rawData(i,:,k) = this.data.res(i,k).export.meas;
                            end
                            this.data.res(i,k).export.rawData = [];
                        end

                        if this.vars.usePhiPP
                            this.data.rawPhi(j,i,:,k)  = this.data.res(i,k).phi;
                        else
                            this.data.rawPhi(j,i,:,k)  = this.data.res(i,k).rawPhi;
                        end                   

                        this.data.fftResSig(j,i,:,k)    = this.data.res(i,k).fittedFFTNorm(:,this.vars.peakIdx);
                        this.data.fftResBkg(j,i,:,k)    = this.data.res(i,k).fittedFFTNorm(:,this.vars.bkgIdx);
                        this.data.fftResBkgAvg(j,i,:,k) = mean(this.data.res(i,k).fittedFFTNorm(:,this.vars.bkgIdxPreCut), 2);
                        this.data.fftResBkgStd(j,i,:,k) = std(this.data.res(i,k).fittedFFTNorm(:,this.vars.bkgIdxPreCut), 0, 2);
                        this.data.rawStd(j,i,k)         = this.data.res(i,k).std;

                        %Extract Net Signal
                        this.data.phiCut(j,i,:,k) = this.data.rawPhi(j,i,this.vars.cutIdxs,k);
                        this.data.phiNet(j,i,:,k) = this.data.rawPhi(j,i,this.vars.keepIdxs,k);
                        this.data.bkg(j,i,:,k)    = this.data.phiNet(j,i,this.vars.bkgIdxPostCut,k);
                        this.data.sig(j,i,:,k)    = this.data.phiNet(j,i,this.vars.sigIdxPostCut,k);

                        % Extract Statistics
                        this.data.maxVal(j,i,k) = max(this.data.phiNet(j,i,:,k));
                        this.data.avgBkg(j,i,k) = mean(this.data.bkg(j,i,:,k));
                        this.data.stdBkg(j,i,k) = std(this.data.bkg(j,i,:,k));

                        this.data.normMax(j,1:i,k) = this.data.maxVal(j,1:i,k)/mean(this.data.maxVal(j,1:i,k));
                        this.data.normStd(j,1:i,k) = this.data.stdBkg(j,1:i,k)/mean(this.data.stdBkg(j,1:i,k));

                        %Calculate SNR
                        this.data.SNR(j,i,k)    = (this.data.maxVal(j,i,k)-this.data.avgBkg(j,i,k))/this.data.stdBkg(j,i,k);
                        this.data.avgSNR(j,i,k) = mean(this.data.SNR(j,1:i,k));
                        this.data.stdSNR(j,i,k) = std(this.data.SNR(j,1:i,k));

                        this.data.timeNormSNR(j,i,k)    = this.data.SNR(j,i,k)/sqrt(this.vars.uVarsAO.ao.timeToSample);
                        this.data.avgTimeNormSNR(j,i,k) = mean(this.data.timeNormSNR(j,1:i,k));
                        this.data.stdTimeNormSNR(j,i,k) = std(this.data.timeNormSNR(j,1:i,k));

                        this.data.normSNR(j, 1:i, k) = this.data.SNR(j,1:i,k)/this.data.avgSNR(j,i,k);

                        this.data.SNRAvgMeas(1,j,k) = this.data.avgSNR(j,i,k);
                        this.data.SNRstdMeas(1,j,k) = this.data.stdSNR(j,i,k);
                        this.data.SNRAvgNorm(1,j,k) = this.data.avgTimeNormSNR(j,i,k);
                        this.data.SNRStdNorm(1,j,k) = this.data.stdTimeNormSNR(j,i,k);
                        
                        if i >= this.vars.subSetSize
                            this.data.SNRsubset(j, i-this.vars.subSetSize+1,k) = mean( this.data.SNR(j,i-this.vars.subSetSize+1:i,k));
                        end
                        
                    end
                    this.data.SNRGain(j,i)    = this.data.SNR(j,i,1)./this.data.SNR(j,i,2);
                    this.data.SNRAvgGain(1,j) = mean(this.data.SNRGain(j,1:i));
                    this.data.SNRStdGain(1,j) = std(this.data.SNRGain(j,1:i));
                    
                    this.plotResultsDual();
                    
                    fprintf("Meas: %d/%d: %d/%d, SNR-PD: %.2f, SNR-PMT: %.2f ; SNRGain %.2f  AvgSNRGain %.2f\n",...
                                j, this.vars.measLen, i, this.vars.reps ,...
                                this.data.SNR(j,i,1), this.data.SNR(j,i,2),...
                                this.data.SNRGain(j,i), this.data.SNRAvgGain(1,j));
                end
                
                %Collect current results:
                if  this.vars.uVarsAO.ao.exportData.meas
                    resMeas.rawData = this.data.rawData;
                end

                resMeas.aoVars  = this.data.curVars;
                resMeas.res     = this.data.res;
                resMeas.rawPhi  = squeeze(this.data.rawPhi(j,:,:,:));

                resMeas.phiCut = squeeze(this.data.phiCut(j,:,:,:));
                resMeas.phiNet = squeeze(this.data.phiNet(j,:,:,:));
                resMeas.bkg    = squeeze(this.data.bkg(j,:,:,:));
                resMeas.sig    = squeeze(this.data.sig(j,:,:,:));

                resMeas.maxVal  = this.data.maxVal(j,:,:);
                resMeas.avgBkg  = this.data.avgBkg(j,:,:);
                resMeas.stdBkg  = this.data.stdBkg(j,:,:);

                resMeas.normMax = this.data.normMax(j,:,:);
                resMeas.normStd = this.data.normStd(j,:,:);

                resMeas.SNR    = this.data.SNR(j,:,:);
                resMeas.avgSNR = this.data.avgSNR(j,:,:);
                resMeas.stdSNR = this.data.stdSNR(j,:,:);

                resMeas.timeNormSNR = this.data.timeNormSNR(j,:,:);
                resMeas.avgTimeNormSNR = this.data.avgTimeNormSNR(j,:,:);
                resMeas.stdTimeNormSNR = this.data.stdTimeNormSNR(j,:,:);

                resMeas.normSNR = this.data.normSNR(j, :, :);
                
                resMeas.SNRGain = this.data.SNRGain(j,:);
                
                % Finalizing Iteration
                str = sprintf("Results: Avg. PD SNR: %.2f, Avg. PMT SNR: %.2f, Avg. SNR Gain: %.2f\n",...
                    this.data.SNRAvgMeas(1,j,1), this.data.SNRAvgMeas(1,j,2), this.data.SNRAvgGain(1,j));
                fprintf(str)

                if this.vars.tgUpdate
                    tgprintf(str);
                    for q=1:3
                        tgprint(this.vars.figs.dual.hFig(q));
                    end
                end
                
                this.saveLoopData(resMeas);

                resMeas.rawData = [];
                resTotal(j) = resMeas;
                this.vars.measTime(j) = toc(this.vars.measTime(j));
                fprintf("Done Iteration : %d, Duration: %.2f \n", j,  this.vars.measTime(j)) 

            end 
            this.data.rawData = [];
            
            exp.vars = this.getVars();
            exp.data = this.data;

            this.saveMeasData(exp, resTotal);
            
            if this.vars.tgUpdate
                str = sprintf("Done Measurement: %s\n", this.vars.measName);
                tgprintf(str);
%                 tgprint(this.vars.figs.hFig);
            end

            fprintf("Done Measurement.\n")
        end
        
        function exp = measureSingle(this)
            this.initDataStructs();
            this.initFigsSingle();
            this.initFS();
            
            for j = 1:this.vars.measLen
                this.vars.idxs(1) = j;
                this.vars.measTime(j) = tic;
                switch this.vars.measType
                    case 'PD'
                        this.stagesT.moveAbsAx('X', this.vars.thetaVec(j))
                    case 'PMT'
                        this.data.dc(j) = this.dialogBox(this.vars.loopVar(j));
                    case 'Point'
                        if this.vars.rotate
                            this.stagesT.moveAbsAx('X', this.vars.thetaVec(j))
                        end
                end
                
                fprintf("Meas: %d\n", j);
                
                this.data.curVars = this.ao.getVars();
                for i=1:this.vars.reps   
                    this.vars.idxs(2) = i;
                    pause(this.vars.stabTime)
                    
%                     this.data.res(i).rawData = [];
                    this.data.res(i)       = this.ao.runAcoustoOptics();

                    if  this.vars.uVarsAO.ao.exportData.rawData
                        if i<= this.vars.saveN
                            this.data.rawData(i, :) = this.data.res(i).rawData;
                        end
                        this.data.res(i).export.rawData = [];
                    end
                    
                    if this.vars.usePhiPP
                        this.data.rawPhi(j,i,:)  = this.data.res(i).phi;
                    else
                        this.data.rawPhi(j,i,:)  = this.data.res(i).rawPhi;
                    end                   
                    
                    this.data.fftResSig(j,i,:)    = this.data.res(i).unFittedFFTShift(:,this.vars.peakIdx);
                    this.data.fftResBkg(j,i,:)    = this.data.res(i).unFittedFFTShift(:,this.vars.bkgIdx);
                    this.data.fftResBkgAvg(j,i,:) = mean(this.data.res(i).unFittedFFTShift(:,this.vars.bkgIdxPreCut), 2);
                    this.data.fftResBkgStd(j,i,:) = std(this.data.res(i).unFittedFFTShift(:,this.vars.bkgIdxPreCut), 0, 2);
                    this.data.rawStd(j,i)         = this.data.res(i).std;
                    
                    %Extract Net Signal
                    this.data.phiCut(j,i,:) = this.data.rawPhi(j,i,this.vars.cutIdxs);
                    this.data.phiNet(j,i,:) = this.data.rawPhi(j,i,this.vars.keepIdxs);
                    this.data.bkg(j,i,:)    = this.data.phiNet(j,i,this.vars.bkgIdxPostCut);
                    this.data.sig(j,i,:)    = this.data.phiNet(j,i,this.vars.sigIdxPostCut);

                    % Extract Statistics
                    this.data.maxVal(j,i) = max(this.data.phiNet(j,i,:));
                    this.data.avgBkg(j,i) = mean(this.data.bkg(j,i,:));
                    this.data.stdBkg(j,i) = std(this.data.bkg(j,i,:));

                    this.data.normMax(j,1:i) = this.data.maxVal(j,1:i)/mean(this.data.maxVal(j,1:i));
                    this.data.normStd(j,1:i) = this.data.stdBkg(j,1:i)/mean(this.data.stdBkg(j,1:i));

                    %Calculate SNR
                    this.data.SNR(j,i)    = (this.data.maxVal(j,i)-this.data.avgBkg(j,i))/this.data.stdBkg(j,i);
                    this.data.avgSNR(j,i) = mean(this.data.SNR(j,1:i));
                    this.data.stdSNR(j,i) = std(this.data.SNR(j,1:i));

                    this.data.timeNormSNR(j,i)    = this.data.SNR(j,i)/sqrt(this.vars.uVarsAO.ao.timeToSample);
                    this.data.avgTimeNormSNR(j,i) = mean(this.data.timeNormSNR(j,1:i));
                    this.data.stdTimeNormSNR(j,i) = std(this.data.timeNormSNR(j,1:i));

                    this.data.normSNR(j, 1:i) = this.data.SNR(j,1:i)/this.data.avgSNR(j,i);

                    this.data.SNRAvgMeas(j) = this.data.avgSNR(j,i);
                    this.data.SNRstdMeas(j) = this.data.stdSNR(j,i);
                    this.data.SNRAvgNorm(j) = this.data.avgTimeNormSNR(j,i);
                    this.data.SNRStdNorm(j) = this.data.stdTimeNormSNR(j,i);
                    this.data.SNRAvgGain(j) = this.data.avgSNR(j,i)./this.vars.refSNR;
                    this.data.SNRStdGain(j) = this.data.stdSNR(j,i)./this.vars.refSNR;

                    if i >= this.vars.subSetSize
                        this.data.SNRsubset(j, i-this.vars.subSetSize+1) = mean( this.data.SNR(j,i-this.vars.subSetSize+1:i));
                    end

                    this.plotResultsSingle()

                    fprintf("Meas: %d/%d: %d/%d, SNR: %.2f, AvgSNR: %.2f ; SNRGain %.2f \n",...
                            j, this.vars.measLen, i, this.vars.reps ,...
                            this.data.SNR(j,i), this.data.avgSNR(j,i),...
                            this.data.SNRAvgGain(j));
                end
                

                %Collect current results:
                if  this.vars.uVarsAO.ao.exportData.rawData
                    resMeas.rawData = this.data.rawData;
                end
                
                resMeas.aoVars  = this.data.curVars;
                resMeas.res     = this.data.res;
                resMeas.rawPhi  = squeeze(this.data.rawPhi(j,:,:));

                resMeas.phiCut = squeeze(this.data.phiCut(j,:,:));
                resMeas.phiNet = squeeze(this.data.phiNet(j,:,:));
                resMeas.bkg    = squeeze(this.data.bkg(j,:,:));
                resMeas.sig    = squeeze(this.data.sig(j,:,:));

                resMeas.maxVa   = this.data.maxVal(j,:);
                resMeas.avgBkg  = this.data.avgBkg(j,:);
                resMeas.stdBkg  = this.data.stdBkg(j,:);

                resMeas.normMax = this.data.normMax(j,:);
                resMeas.normStd = this.data.normStd(j,:);

                resMeas.SNR    = this.data.SNR(j,:);
                resMeas.avgSNR = this.data.avgSNR(j,:);
                resMeas.stdSNR = this.data.stdSNR(j,:);

                resMeas.timeNormSNR = this.data.timeNormSNR(j,:);
                resMeas.avgTimeNormSNR = this.data.avgTimeNormSNR(j,:);
                resMeas.stdTimeNormSNR = this.data.stdTimeNormSNR(j,:);

                resMeas.normSNR = this.data.normSNR(j, :);

                str = sprintf("Results: Average SNR: %.2f, STD: %.2f ; Average Norm SNR: %.2f, STD: %.2f\n",...
                    this.data.SNRAvgMeas(j), this.data.SNRstdMeas(j), this.data.SNRAvgNorm(j), this.data.SNRStdNorm(j));
                fprintf(str)

                if this.vars.tgUpdate
                    tgprintf(str);
                    tgprint(this.vars.figs.hFig);
                end
                this.saveLoopData(resMeas);

                resMeas.rawData = [];
                resTotal(j) = resMeas;
                this.vars.measTime(j) = toc(this.vars.measTime(j));
                fprintf("Done Iteration : %d, Duration: %.2f \n", j,  this.vars.measTime(j)) 

            end 
            this.data.rawData = [];
            
            exp.vars = this.getVars();
            exp.data = this.data;

            this.saveMeasData(exp, resTotal);
            
            if this.vars.tgUpdate
                str = sprintf("Done Measurement: %s\n", this.vars.measName);
                tgprintf(str);
%                 tgprint(this.vars.figs.hFig);
            end

            fprintf("Done Measurement.\n")
        end
        
        function saveLoopData(this, data)
            j = this.vars.idxs(1);
            timeStamp = strrep(datestr(datetime('now')),':','-');
            filename = sprintf('%s/%s-%s-%.2f.mat', this.vars.measDir, timeStamp, this.vars.measName, this.vars.loopVar(j));
            a = tic;
            if this.vars.saveFlag
                fprintf("Saving raw data...\n");
                save(filename, 'data', '-v7.3');
            end
            fprintf("Done Saving: %.2f[s] \n", toc(a))
        end
        
        function saveMeasData(this, expData, loopsData)
            timeStamp = strrep(datestr(datetime('now')),':','-');
            filename = sprintf('%s/%s-%s-ExpVars.mat', this.vars.measDir, timeStamp, this.vars.measName);
            if this.vars.saveFlag
                save(filename, 'expData', 'loopsData', '-v7.3');
            end
        end
        
        function data = getData(this)
           data = this.data;            
        end
        
    end
end