classdef AOGraphics < Graphics
    %AOGRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data
        dataAllCh
        usSignal
        extClkSignal
        croppedSig
        
        tVecUS
        tVecSampleClk
        tVecMeas
        tVecSig
        tVecPos
        fBar
        
        depthNorm
        depthCntr
        depthZero
        depthIdx
        depthNormRaw
    end
    
    methods (Static)
        function gNames = getGraphicsNames()
            gNames = {'extClk';'usSignal';...
                'cropped';'meas';'signal'; 'deMul'; 'reshaped';...
                'rawFFT'; 'calibration'; 'fittingModel'; 'fittedPowerFFT';'finalFFT';...
                'phi'; 'rawPhi'; 'phiLog'; };
        end
        
        function figs  = createGraphicsVars()
            figs.depthIdx      = 1;
            figs.ch            = 1;
            figs.frame         = 1;
            figs.numOfChannels = 1;
            
            figs.depthAxType    = 'Normal';
            figs.fitToShift     = true;
            figs.fIdx           = [];
            figs.displayFullFFT = true;
            figs.FFTenv         = [];
            figs.FFTIdxs        = [];
            
            figs.intExt  = 'int';
            figsNames = AOGraphics.getGraphicsNames();

            figs.fonts.type       = [];
            figs.fonts.titleSize  = 14;
            figs.fonts.labelsSize = 14;
            figs.fonts.axisSize   = 14;
            
            figs.analyzeSingleCh = false;
            figs.singlChIdx = 1;
            figs.dispMuEff  = false;
            
            for i=1:length(figsNames)
                figs.validStruct.(figsNames{i}) = false;
            end
            
            for i = 1:length(figsNames)  
                figs.(figsNames{i}).type      = []; 
                
                figs.(figsNames{i}).handles.int = Graphics.createHandlesStruct();
                figs.(figsNames{i}).handles.ext = Graphics.createHandlesStruct();
                figs.(figsNames{i}).handles.cur = Graphics.createHandlesStruct();
                
                figs.(figsNames{i}).strings.title       = [];
                figs.(figsNames{i}).strings.titleModel  = [];
                figs.(figsNames{i}).strings.xlabel      = [];
                figs.(figsNames{i}).strings.ylabel      = [];
                figs.(figsNames{i}).strings.legend      = [];
                figs.(figsNames{i}).strings.legendModel = [];
                figs.(figsNames{i}).strings.updateLegend = true;

                figs.(figsNames{i}).dims.zDim  = [];
                figs.(figsNames{i}).dims.chDim   = [];
                figs.(figsNames{i}).dims.dataDim = [];

                figs.(figsNames{i}).markers.plotMark = false;

                figs.(figsNames{i}).lims.xlims  = [];
                figs.(figsNames{i}).lims.ylims  = [];
                figs.(figsNames{i}).lims.clims = [];
            end
        end 
        
        function uVars = createUserVars()
            %Variables that comes from the user
            figsNames = AOGraphics.getGraphicsNames();
            
            uVars.ch    = 1;
            uVars.depthIdx  = 1;
            uVars.frame = 1;
            
            uVars.singleChIdx = 1;
            
            uVars.reopenFigures = false;
            uVars.depthAxType = 'Normal';
            uVars.displayFullFFT = false;
            uVars.FFTenv = [];
            
            uVars.intExt = [];
            for i=1:length(figsNames)
                uVars.validStruct.(figsNames{i}) = false;
            end
            
            for i=1:length(figsNames)
                uVars.extH.(figsNames{i}) =  Graphics.createHandlesStruct();
            end
        end
        
        function vars  = createOwnerVars()
            %Variables that comes from the AO Object
            vars = AOGraphics.createUserVars();
            
            vars.numOfChannels   = [];
            vars.df              = [];
            vars.fIdx            = [];
            vars.fitToShift      = [];
            
            vars.tUS       = [];
            vars.tExtClk   = [];
            vars.tMeas     = [];
            vars.tMeasSig  = [];
            vars.tSig      = [];
            vars.tPos      = [];
            vars.depthVec  = [];
            vars.fBar      = [];
            
            vars.depthNorm   = [];
            vars.depthCntr   = [];
            vars.depthZero   = [];
            vars.depthIdx    = [];
            
            vars.depthNormRaw = [];
        end 
    end
    
    methods
        function this = AOGraphics()
            this@Graphics()
            
            this.figsNames = this.getGraphicsNames();
            this.numOfFigs = length(this.figsNames);
            
            this.figs     = this.createGraphicsVars();
            
            this.uVars    = this.createOwnerVars();
            
            this.setGraphicsStaticVars();
        end
        
        % Set Functions
        function setUserVars(this, uVars)
            this.uVars.depthIdx  = uVars.depthIdx;
            this.uVars.frame = uVars.frame;
            this.uVars.ch    = uVars.ch;
            
            this.uVars.displayFullFFT = uVars.displayFullFFT;
            this.uVars.FFTenv         = uVars.FFTenv;
            this.uVars.depthAxType    = uVars.depthAxType;
            this.uVars.reopenFigures  = uVars.reopenFigures;
            
            this.uVars.intExt      = uVars.intExt;
            this.uVars.validStruct = uVars.validStruct;
            this.uVars.extH        = uVars.extH;
            
            this.uVars.numOfChannels  = uVars.numOfChannels;
            this.uVars.df             = uVars.df;
            this.uVars.fIdx           = uVars.fIdx;
            
            this.tVecUS         = uVars.tVecUS*1e6;
            this.tVecSampleClk  = uVars.tVecSampleClk*1e6;
            
            this.tVecMeas   = uVars.tVecMeas*1e6;
            this.tVecSig    = uVars.tVecSig*1e6;
            this.tVecPos    = uVars.tVecPos*1e6;
            
            this.fBar        = uVars.fBar*1e-6;
            
            this.depthNorm   = uVars.depthVecNorm*1e3;
            this.depthCntr   = uVars.depthVecCntr*1e3;
            this.depthZero   = uVars.depthVecZero*1e3;
            this.depthIdx    = uVars.depthVecIdx;
            
            this.depthNormRaw = uVars.depthVecNormRaw;
%             this.depthCntrRaw = uVars.depthVecNormRaw;
%             this.depthZeroRaw = uVars.depthVecNormRaw;
%             this.depthIdxRaw  = uVars.depthVecNormRaw;
            
            this.uVars.analyzeSingleCh = uVars.analyzeSingleCh;
            this.uVars.singleChIdx     = uVars.singleChIdx;
            this.uVars.dispMuEff       = uVars.dispMuEff;
        end
        
        function setGraphicsStaticVars(this)
%             this.graphics.setType(type);
%             this.graphics.setStrings(title, xlabel, ylabel, legend);
%             this.graphics.setDims(gName, zIdxDim, chDim, dataDim)
%             this.graphics.setMarkersEnable(enable);          

            % extClk
            this.setType('extClk', 'stem');
            this.setStrings('extClk', "External Clk (Sampling Clock)", "t[\mus]", "Amp", []);
            
            % usSignal
            this.setType('usSignal', 'stem');
            this.setStrings('usSignal', "Ultrasound Pulse", "t[\mus]", "Amp", []);
           
            %-------------------
            % Time Domain Signal
            %-------------------
            % cropped [ch x samplesPerMeas]
            this.setType('cropped', 'plot');
            this.setStrings('cropped', "Cropped Samples Channel: %d", "t[\mus]", "Amp[V]", []);
            
            % meas [ch x samplesPerMeas]
            this.setType('meas', 'plot');
            this.setStrings('meas', "Measured Signal Channel: %d", "t[\mus]", "Amp[V]", []);
            
            % signal [ch x samplesPerSignal]
            this.setType('signal', 'plot');
            this.setStrings('signal', "Net Signal Channel: %d", "t[\mus]", "Amp[V]", []);
            
            % deMultiplexed
            this.setType('deMul', 'plot');
            this.setStrings('deMul', "Demultiplexed Signal Channel: %d", "t[\mus]", "Amp[V]", []);
            
            % reshaped [ch x samplesPerPos x numOfPos]
            this.setType('reshaped', 'stem');
            this.setStrings('reshaped', "Reshaped Signal Channel: %d, Depth: %.2f[mm] (idx:(%d))", "t[\mu s]", "Amp[V]", []);
            
            %-------------------
            % Frequency Domain Signal
            %-------------------
            % rawFFT 
            this.setType('rawFFT', 'plot');
            this.setStrings('rawFFT', "Raw FFT (Frame-Avg.): Depth: %.2f[mm] (idx: %d)", "f[MHz]", "Power Spectrum", "Ch %d");
            this.setMarkersEnable('rawFFT', true);
            
            % Fitting Model 
            this.setType('fittingModel', 'plot');
            this.setStrings('fittingModel', "Fitting Model Per-Channel: Ch: %d, Depth: %.2f[mm] (idx: %d)", "f[MHz]", "Power Spectrum", []);
            this.figs.fittingModel.strings.legend = {'Raw FFT'; 'Fit'};
%             this.figs.fittingModel.strings.updateLegend = false;
            this.setMarkersEnable('fittingModel', true);
            
            % fittedPowerFFT
            this.setType('fittedPowerFFT', 'plot');
            this.setStrings('fittedPowerFFT', "Fitted FFT (Frame-Avg.): Depth: %.2f[mm] (idx: %d)", "f[MHz]", "Power Spectrum", "Ch %d");
            this.setMarkersEnable('fittedPowerFFT', true);
            
            % finalFFT
            this.setType('finalFFT', 'plot');
            this.setStrings('finalFFT', "Final FFT (Ch.-Avg.): Depth: %.2f[mm] (idx: %d)", "f[MHz]", "Power Spectrum", []);
            this.setMarkersEnable('finalFFT', true);
            
            %-------------------
            % Spatial Recontstruction
            %-------------------
            
            % phi [numOfPos]
            this.setType('phi', 'stem');
            this.setStrings('phi', [], "Depth (US)[Idx]", "\phi", []);
            
            % phiRaw [numOfPos]
            this.setType('rawPhi', 'stem');
            this.setStrings('rawPhi', [], "Depth (US)[Idx]", "\phi", []);
            
            % phi [numOfPos]
            this.setType('phiLog', 'plot');
            this.setStrings('phiLog', [], "Depth (US)[Idx]", "log(\phi)", ["AlgoPhi"; "RawPhi"]);
            
        end
        
        function setGraphicsDynamicVars(this)
            %this function assums that uVars are updated
            this.setChAndPos(this.uVars.ch, this.uVars.depthIdx, this.uVars.frame)
            
            this.figs.fIdx           = this.uVars.fIdx;
            this.figs.displayFullFFT = this.uVars.displayFullFFT;
            this.figs.FFTenv         = this.uVars.FFTenv;
            envIdx                   = floor(this.uVars.FFTenv / this.uVars.df);
            this.figs.FFTIdxs        = this.uVars.fIdx + [-envIdx, envIdx];
            this.figs.depthAxType     = this.uVars.depthAxType;
        
            numOfCh   = this.uVars.numOfChannels;
            if this.figs.numOfChannels ~= numOfCh
                this.figs.rawFFT.update = true;
                this.figs.fittedPowerFFT.update = true;
            end
            this.figs.numOfChannels = numOfCh;
            this.setLegendVariables('rawFFT',  1:numOfCh);
            this.setLegendVariables('fittedPowerFFT', 1:numOfCh);
            
            this.updateGraphicsConstruction()
            this.figs.validStruct = this.uVars.validStruct; % has to be after updating the construction.
%             if this.uVars.frame > 1
%                 this.setType('phi', 'stem');
%             else
%                 this.setType('phi', 'stem');
%             end
            
            this.figs.analyzeSingleCh = this.uVars.analyzeSingleCh;
            this.figs.singleChIdx     = this.uVars.singleChIdx;
            this.figs.dispMuEff       = this.uVars.dispMuEff;
            
            if this.figs.analyzeSingleCh
                this.figs.ch = 1;
            end

        end
        
        function setData(this, data, usSig, extClkSig)
            this.dataAllCh = data;
            if ~isempty(data)
                if this.figs.analyzeSingleCh
                    this.data = data(this.figs.singleChIdx);
                else
                    this.data = data;
                end
            end
                
            if nargin > 2
                this.usSignal = usSig;
                this.extClkSignal = extClkSig;
            end
        end
        
        function setSingleChIdx(this, idx)
           this.figs.singleChIdx = idx;
           if this.figs.analyzeSingleCh
               if ~isempty(this.dataAllCh) 
                    this.data = this.dataAllCh(this.figs.singleChIdx);
               end
               this.figs.ch = 1;
           end
        end
        
        function setCroppedData(this, data)
            this.croppedSig = data;
        end
        
        function depthVec = getDepthAxByType(this)
           switch this.figs.depthAxType
               case 'Normal'
                   depthVec = this.depthNorm;
               case 'Center'
                   depthVec = this.depthCntr;
               case 'Zero'
                   depthVec = this.depthZero;
               case 'Index'
                   depthVec = this.depthIdx;
           end
        end
        
        function signal = extractSignal(this, figName)
            ch = this.figs.ch;
            switch figName
                case 'extClk'
                    signal.xData = this.tVecSampleClk;
                    signal.yData = this.extClkSignal;
                case 'usSignal'
                    signal.xData = this.tVecUS;
                    signal.yData = this.usSignal;
                case 'cropped'
                    signal.xData = 1:size(this.croppedSig, 2);
                    signal.yData = this.croppedSig(ch, :);
                case 'meas'
                    signal.ch    = this.figs.ch;
                    signal.xData = this.tVecMeas;
                    signal.yData = this.data.export.meas(ch, :);
                case 'signal'
                    signal.ch    = this.figs.ch;
                    signal.xData = this.tVecSig;
                    signal.yData = this.data.export.signal(ch, :);
                case 'deMul'
                    signal.ch    = this.figs.ch;
                    signal.xData = this.tVecSig;
                    signal.yData = this.data.export.deMul(ch, :); 
                case 'reshaped'
                    signal.ch        = this.figs.ch;
                    signal.depthIdx  = this.figs.depthIdx;
                    signal.frame     = this.figs.frame;
                    signal.xData     = this.tVecPos;
                    signal.yData     = permute( this.data.export.reshaped(signal.frame ,signal.ch, :, signal.depthIdx), [2,3,1,4]);
            end
        end
        
        function signal = extractPhi(this, figName)
            switch figName
                case 'phi'
                    switch this.figs.depthAxType
                       case 'Normal'
                           signal.xData{1} = this.data.analysis.phi.depthVecBkg;
                           signal.xData{2} = this.data.analysis.phi.depthVecSignal;
                       case 'Index'
                           signal.xData{1} = this.data.analysis.phi.bkgIdx;
                           signal.xData{2} = this.data.analysis.phi.signalIdx;
                    end
                    signal.yData{1} = this.data.analysis.phi.bkg;
                    signal.yData{2} = this.data.analysis.phi.signal;
                    signal.SNR      = this.data.analysis.phi.SNR;
                    signal.noiseStd = this.data.analysis.phi.stdNoise;
                    if this.figs.dispMuEff
                        signal.muEffVal    = this.data.analysis.muEff.muEffVal;
                        signal.muEffRawVal = this.data.analysis.muEff.muEffRawVal;
                    end
                    minVal = min(min(signal.yData{1}(signal.yData{1}~=-inf)), min(signal.yData{2}(signal.yData{2}~=-inf)));
                    maxVal = max(max(signal.yData{1}(signal.yData{1}~=inf)), max(signal.yData{2}(signal.yData{2}~=inf)));
                case 'rawPhi'
                    switch this.figs.depthAxType
                       case 'Normal'
                           signal.xData{1} = this.data.analysis.rawPhi.depthVecBkg;
                           signal.xData{2} = this.data.analysis.rawPhi.depthVecSignal;
                       case 'Index'
                           signal.xData{1} = this.data.analysis.rawPhi.bkgIdx;
                           signal.xData{2} = this.data.analysis.rawPhi.signalIdx;
                    end
                    signal.yData{1} = this.data.analysis.rawPhi.bkg;
                    signal.yData{2} = this.data.analysis.rawPhi.signal;
                    signal.SNR      = this.data.analysis.rawPhi.SNR;
                    signal.noiseStd = this.data.analysis.rawPhi.stdNoise;
                    if this.figs.dispMuEff
                        signal.muEffVal    = this.data.analysis.muEff.muEffVal;
                        signal.muEffRawVal = this.data.analysis.muEff.muEffRawVal;
                    end
                    minVal = min(min(signal.yData{1}(signal.yData{1}~=-inf)), min(signal.yData{2}(signal.yData{2}~=-inf)));
                    maxVal = max(max(signal.yData{1}(signal.yData{1}~=inf)), max(signal.yData{2}(signal.yData{2}~=inf)));
                case 'phiLog'
                    signal.xData{1}    = this.getDepthAxByType();
                    signal.xData{2}    = this.getDepthAxByType();
                    signal.yData{1}    = this.data.phiLog;
                    signal.yData{2}    = this.data.analysis.rawPhi.normTypes.phiLog2;
                    if this.figs.dispMuEff
                        switch this.figs.depthAxType
                            case 'Normal'
                                signal.xData{3}    = this.data.analysis.muEff.xVecInt;
                                signal.xData{4}    = this.data.analysis.muEff.xVecInt;
                            case 'Index'
                                signal.xData{3}    = this.data.analysis.muEff.xVecIntIdx;
                                signal.xData{4}    = this.data.analysis.muEff.xVecIntIdx;
                        end
                        signal.yData{3}    = this.data.analysis.muEff.phiCutInt;
                        signal.yData{4}    = this.data.analysis.muEff.phiRawCutInt;
                        signal.muEffVal    = this.data.analysis.muEff.muEffVal;
                        signal.muEffRawVal = this.data.analysis.muEff.muEffRawVal;
                    end
                    signal.SNR      = this.data.analysis.phi.SNR;
                    signal.noiseStd = this.data.analysis.phi.stdNoise;
                    minVal = min(signal.yData{1}(signal.yData{1}~=-inf));
                    maxVal = max(signal.yData{1});
            end
            span = maxVal - minVal;
            signal.minLim = minVal - 0.1*span;
            signal.maxLim = maxVal + 0.1*span;
            signal.SNRYlim = double(signal.maxLim - 0.15*span);
        end

        function signal = extractFFTSignal(this, figName)
            signal.ch = this.figs.ch;
            signal.depthIdx  = this.figs.depthIdx;
            signal.depthVec  = this.getDepthAxByType();
            if ~this.figs.displayFullFFT
                % plot user defined enviroment around the signal
                signal.fidxs = this.figs.FFTIdxs(1):this.figs.FFTIdxs(2);
            else
                signal.fidxs = 1:length(this.fBar);
            end
            this.figs.(figName).lims.xlims = [this.fBar(signal.fidxs(1)), this.fBar(signal.fidxs(end))]; 
            signal.xData = this.fBar(signal.fidxs);
            switch figName               
                case 'rawFFT'
                    signal.yData     = this.data.frameAvgPowerFFTCut(:, signal.fidxs, signal.depthIdx);
%                     signal.yData     = squeeze(signal.yData(:,signal.fidxs));
                    signal.markerY   = squeeze(this.data.frameAvgPowerFFTCut(:,this.figs.fIdx, signal.depthIdx));
                    signal.titleVars = {{signal.depthVec(signal.depthIdx), signal.depthIdx}};
                case 'fittingModel'
                    signal.yData(1,:,1) = permute(this.data.frameAvgPowerFFTCut(signal.ch, signal.fidxs, signal.depthIdx), [3,1,2]);
                    signal.yData(2,:,1) = permute(this.data.fitModelMat(signal.ch, signal.fidxs), [3,1,2]);
                    signal.titleVars    = {{signal.ch, signal.depthVec(signal.depthIdx), signal.depthIdx}};
                    signal.markerY(1)   = squeeze(this.data.frameAvgPowerFFTCut(signal.ch, this.figs.fIdx, signal.depthIdx));
                    signal.markerY(2)   = squeeze(this.data.fitModelMat(signal.ch, this.figs.fIdx));
                case 'fittedPowerFFT'
                    signal.yData     = this.data.fittedPowerFFT(:, signal.fidxs, signal.depthIdx);
                    signal.titleVars = {{signal.depthVec(signal.depthIdx), signal.depthIdx}};
                    signal.markerY   = squeeze(this.data.fittedPowerFFT(:, this.figs.fIdx, signal.depthIdx));
                case 'finalFFT'
                    signal.yData     = this.data.fittedFFT(signal.fidxs, signal.depthIdx)';
%                     signal.yData     = squeeze(signal.yData(:,signal.fidxs));
                    signal.markerY   = squeeze(this.data.fittedFFT(this.figs.fIdx, signal.depthIdx));
                    signal.titleVars = {{signal.depthVec(signal.depthIdx), signal.depthIdx}};
            end
            signal.markerX = this.fBar(this.figs.fIdx);
        end
        
        % Disply Functions
        function quickPlot(this, gName, xData, yData, std, markerX, markerY, hIdx) 
            if ~exist('hIdx', 'var')
                hIdx = 1;
            end
            switch this.figs.(gName).type
                case 'stem'
                    set(this.figs.(gName).handles.cur.plot(hIdx),...
                        'XData', xData,...
                        'YData', gather(yData'));
                case 'plot'
                    set(this.figs.(gName).handles.cur.plot(hIdx),...
                        'XData', xData, ...
                        'YData', gather(yData'));
                case 'errorbar'
                    set(this.figs.(gName).handles.cur.plot(hIdx),...
                        'XData', xData,...
                        'YData', yData,...
                        'YPositiveDelta', std,...
                        'YNegativeDelta', std);
            end
            if this.figs.(gName).markers.plotMark
                set(this.figs.(gName).handles.cur.plotMarker(hIdx), 'XData', markerX, 'YData', gather(markerY));
            end
            set(this.figs.(gName).handles.cur.title, 'String', this.figs.(gName).strings.title);
%             pause(0.01)
%             drawnow();
        end
        
        function plotAFGSignals(this, figName)
            %yData - [1 x Sig]
            if ~isgraphics(this.figs.(figName).handles.cur.ax) || ~this.figs.validStruct.(figName)
               return
            end
            
            signal = this.extractSignal(figName);
            
            if ~isgraphics(this.figs.(figName).handles.cur.plot)
                switch this.figs.(figName).type
                    case 'stem'
                       this.figs.(figName).handles.cur.plot = ...
                           stem(this.figs.(figName).handles.cur.ax,...
                           signal.xData, signal.yData);    
                    case 'plot'
                       this.figs.(figName).handles.cur.plot = ...
                           plot(this.figs.(figName).handles.cur.ax,...
                           signal.xData, signal.yData);  
                end
                setStringsToPlot(this, figName)
            else
                quickPlot(this, figName, signal.xData, signal.yData,[], [], [])
                pause(0.01);
                drawnow();
            end
        end
                
        function plotPhi(this, figName)
            %yData - [samples per pos, 1]
            
            if ~isgraphics(this.figs.(figName).handles.cur.ax)  || ~this.figs.validStruct.(figName)
               return
            end
            
            signal = this.extractPhi(figName);
            if this.figs.dispMuEff
               SNRstr = sprintf("SNR: %.2f\nBkg STD: %.2E\n\\mu_{Eff} = %.3f\n\\mu_{Eff}^{Raw} = %.3f", signal.SNR, signal.noiseStd, signal.muEffVal, signal.muEffRawVal);
            else
               SNRstr = sprintf("SNR: %.2f\nBkg STD: %.2E", signal.SNR, signal.noiseStd);
            end
            
            if ~isgraphics(this.figs.(figName).handles.cur.plot(1)) ||...
                ~strcmp(this.figs.(figName).type, this.figs.(figName).handles.cur.plot(1).Type)  
                this.figs.(figName).handles.cur.plot = gobjects(1,size(signal.yData,1));
                for i = 1:size(signal.yData, 2)
                    if i == 2
                        hold(this.figs.(figName).handles.cur.ax, 'on');
                    end
                    switch this.figs.(figName).type
                        case 'stem'
                               this.figs.(figName).handles.cur.plot(i) = ...
                                   stem(this.figs.(figName).handles.cur.ax,...
                                   signal.xData{i}, signal.yData{i});
                        case 'plot'
                           this.figs.(figName).handles.cur.plot(i) = ...
                               plot(this.figs.(figName).handles.cur.ax,...
                               signal.xData{i}, signal.yData{i});
                       case 'errorbar'
                           this.figs.(figName).handles.cur.plot(i) = ...
                               errorbar(this.figs.(figName).handles.cur.ax,...
                               signal.xData{i}, signal.yData{i}, signal.std);
                    end
                end
                hold(this.figs.(figName).handles.cur.ax, 'off');
                this.figs.(figName).handles.cur.text = ...
                    text(this.figs.(figName).handles.cur.ax, 0, signal.SNRYlim,...
                    SNRstr, 'FontSize', 14);
                this.setStringsToPlot(figName);
                ylim(this.figs.(figName).handles.cur.ax, [signal.minLim, signal.maxLim]);
            else
                for i = 1:size(signal.yData, 2)
                    quickPlot(this, figName, signal.xData{i}, signal.yData{i}, [], [], [], i);
                end
                set(this.figs.(figName).handles.cur.text, 'String', SNRstr, 'Position', [0, signal.SNRYlim]);
                pause(0.01);
                ylim(this.figs.(figName).handles.cur.ax, [signal.minLim, signal.maxLim]);
                drawnow();
            end
        end
     
        function plotSignal(this, figName)
        % yData - [ch x Sig]
            
            if ~isgraphics(this.figs.(figName).handles.cur.ax)  || ~this.figs.validStruct.(figName)
               return
            end
            signal = extractSignal(this, figName);

            if ~isgraphics(this.figs.(figName).handles.cur.plot)
                switch this.figs.(figName).type
                    case 'stem'
                       this.figs.(figName).handles.cur.plot = ...
                           stem(this.figs.(figName).handles.cur.ax,...
                           signal.xData, signal.yData);    
                    case 'plot'
                       this.figs.(figName).handles.cur.plot = ...
                           plot(this.figs.(figName).handles.cur.ax,...
                           signal.xData, signal.yData);  
                end
                this.setStringsToPlot(figName);
            else
                quickPlot(this, figName, signal.xData, signal.yData,[], [], [])
                pause(0.01);
                drawnow();
            end
        end
        
        function plotReshapedSignal(this) 
        % yData - [ch x Sig x pos]
            
            if ~isgraphics(this.figs.reshaped.handles.cur.ax)  || ~this.figs.validStruct.reshaped
               return
            end
            
            signal = this.extractSignal('reshaped');
            
            if ~isgraphics(this.figs.reshaped.handles.cur.plot)
                switch this.figs.reshaped.type
                    case 'stem'
                       this.figs.reshaped.handles.cur.plot = ...
                           stem(this.figs.reshaped.handles.cur.ax,...
                           signal.xData, signal.yData);    
                    case 'plot'
                       this.figs.reshaped.handles.cur.plot = ...
                           plot(this.figs.reshaped.handles.cur.ax,...
                           signal.xData, signal.yData);  
                end
                this.setStringsToPlot('reshaped');
            else
                quickPlot(this, 'reshaped', signal.xData, signal.yData,[], [], [])
                pause(0.01);
                drawnow()
            end
        end
        
        function plotFFT(this, figName)
            if ~isgraphics(this.figs.(figName).handles.cur.ax) || ~this.figs.validStruct.(figName)
               return
            end
%             fprintf("Plotting: %s\n", figName);
            signal = this.extractFFTSignal(figName);
            this.setTitleVariables(figName, signal.titleVars);          % FFT
 
            if ~isgraphics(this.figs.(figName).handles.cur.plot(1))...
                || this.figs.(figName).update
                
                for i = 1:size(signal.yData, 1)
                    if i == 2
                        hold(this.figs.(figName).handles.cur.ax, 'on');
                    end
                    this.figs.(figName).handles.cur.plot(i) = ...
                        plot(this.figs.(figName).handles.cur.ax,...
                        signal.xData, signal.yData( i, :));  
                end
                hold(this.figs.(figName).handles.cur.ax, 'on');
                for i = 1:length(signal.markerY)            
                    this.plotMarkers(figName, signal.markerX, signal.markerY(i), i)
                end
                hold(this.figs.(figName).handles.cur.ax, 'off');
                
                this.setStringsToPlot(figName)
                this.setLimsToPlot(figName)
                this.figs.(figName).update = false;
            else
                for i = 1:size(signal.yData, 1)
                    this.quickPlot(figName, signal.xData, signal.yData(i, :),[], signal.markerX, signal.markerY(i), i)
                end
                this.setLimsToPlot(figName)
                pause(0.01);
                drawnow()
            end
        end
        
        function plotMarkers(this, gName, xMarker, yMarker, hIdx)
            if this.figs.(gName).markers.plotMark 
                if ~exist('hIdx', 'var')
                    hIdx = 1;
                end
                for i =1:size(yMarker,1)
                    this.figs.(gName).handles.cur.plotMarker(hIdx) = ...
                                    plot(this.figs.(gName).handles.cur.ax,...
                                         xMarker, yMarker, 'g+');
                end
            end
        end
    end
end
