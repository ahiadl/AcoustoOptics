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
                'unFittedFFTAllCh'; 'unFittedFFTSingleCh'; 'fittedFFTSingleCh';'fittedFFTAllCh'; 'avgFFT'; 'normFFT';...
                'phi'; 'rawPhi'; 'phiLog'};
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
            
            this.uVars.analyzeSingleCh = uVars.analyzeSingleCh;
            this.uVars.singleChIdx     = uVars.singleChIdx;
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
            % unFittedFFT All Ch 
            this.setType('unFittedFFTAllCh', 'plot');
            this.setStrings('unFittedFFTAllCh', "UnFitted Frame-Avg FFT (All Ch): Depth: %.2f[mm] (idx: %d)", "f[MHz]", "Power Spectrum", "Ch %d");
            this.setMarkersEnable('unFittedFFTAllCh', false);
            
            % unFittedFFT Single Ch 
            this.setType('unFittedFFTSingleCh', 'plot');
            this.setStrings('unFittedFFTSingleCh', "UnFitted Frame-Avg FFT: Ch: %d, Depth: %.2f[mm] (idx: %d)", "f[MHz]", "Power Spectrum", []);
            this.figs.unFittedFFTSingleCh.strings.legend = {'Raw FFT'; 'Fit'};
            this.figs.unFittedFFTSingleCh.strings.updateLegend = false;
            this.setMarkersEnable('unFittedFFTSingleCh', true);
            
            % fittedFFT Single Ch
            this.setType('fittedFFTSingleCh', 'plot');
            this.setStrings('fittedFFTSingleCh', "Fitted FFT Frame-Avg: Ch: %d, Depth: %.2f[mm] (idx: %d)", "f[MHz]", "Power Spectrum", []);
            this.setMarkersEnable('fittedFFTSingleCh', true);
            
            % fittedFFT All Ch
            this.setType('fittedFFTAllCh', 'plot');
            this.setStrings('fittedFFTAllCh', "Fitted Frame-Avg FFT (All Ch): Depth: %.2f[mm] (idx: %d)", "f[MHz]", "Power Spectrum", "Ch %d");
            this.setMarkersEnable('fittedFFTAllCh', false);
            
            % Avg FFT
            this.setType('avgFFT', 'plot');
            this.setStrings('avgFFT', "Average FFT: Depth: %.2f[mm] (idx: %d)", "f[MHz]", "Power Spectrum", []);
            this.setMarkersEnable('avgFFT', true);
            
            % normalized FFT
            this.setType('normFFT', 'plot');
            this.setStrings('normFFT', "Normalized FFT: Depth: %.2f[mm] (idx: %d)", "f[MHz]", "Power Spectrum", []);
            this.setMarkersEnable('normFFT', true);
            
            %-------------------
            % Spatial Recontstruction
            %-------------------
            
            % phi [numOfPos]
            this.setType('phi', 'stem');
            this.setStrings('phi', "\\phi", "Depth (US)[mm]", "\\phi", []);
            
            % phiRaw [numOfPos]
            this.setType('rawPhi', 'stem');
            this.setStrings('rawPhi', "\\phi (Raw)", "Depth (US)[mm]", "\\phi", []);
            
            % phi [numOfPos]
            this.setType('phiLog', 'plot');
            this.setStrings('phiLog', "log(\\phi)", "Depth (US)[mm]", "log(\\phi)", []);
            
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
                this.figs.qAvgChFFT.update = true;
            end
            this.figs.numOfChannels = numOfCh;
            this.setLegendVariables('unFittedFFTAllCh',  1:numOfCh);
            this.setLegendVariables('fittedFFTAllCh', 1:numOfCh);
            
            this.updateGraphicsConstruction()
            if this.uVars.frame > 1
                this.setType('phi', 'stem');
            else
                this.setType('phi', 'stem');
            end
            
            this.figs.analyzeSingleCh = this.uVars.analyzeSingleCh;
            this.figs.singleChIdx     = this.uVars.singleChIdx;
            
            if this.figs.analyzeSingleCh
                this.figs.ch = 1;
            end
            
            this.figs.validStruct = this.uVars.validStruct;
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
                case 'phi'
                    signal.xData    = this.getDepthAxByType();
                    signal.yData    = this.data.phi;
                    signal.SNR      = this.data.SNR.val;
                    signal.noiseStd = this.data.SNR.noiseStd;
                case 'rawPhi'
                    signal.xData    = this.getDepthAxByType();
                    signal.yData    = this.data.rawPhi;
                    signal.SNR      = this.data.SNR.valRaw;
                    signal.noiseStd = this.data.SNR.noiseStdRaw;
                case 'phiLog'
                    signal.xData    = this.getDepthAxByType();
%                     signal.yData    = real(log(this.data.phi-min(this.data.phi)/(max(this.data.phi) - min(this.data.phi)))); %TDOD: calculate log phi in algo.
                    signal.yData    = this.data.phiLog;
                    signal.SNR      = this.data.SNR.val;
                    signal.noiseStd = this.data.SNR.noiseStd;
            end
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
                case 'unFittedFFTAllCh'
                    signal.yData     = this.data.qAvgChFFT(:, signal.fidxs, signal.depthIdx);
                    signal.markerY   = squeeze(this.data.qAvgChFFT(:,this.figs.fIdx, signal.depthIdx));
                    signal.yData     = squeeze(signal.yData(:,signal.fidxs));
                    signal.titleVars = {{signal.depthVec(signal.depthIdx), signal.depthIdx}};
                case 'unFittedFFTSingleCh'
                    signal.yData(1,:,1) = permute(this.data.qAvgChFFT(signal.ch, signal.fidxs, signal.depthIdx), [3,1,2]);
                    signal.yData(2,:,1) = permute(this.data.fitModelMat(signal.ch, signal.fidxs, signal.depthIdx), [3,1,2]);
                    signal.titleVars    = {{signal.ch, signal.depthVec(signal.depthIdx), signal.depthIdx}};
                    signal.markerY(1)   = squeeze(this.data.qAvgChFFT(signal.ch, this.figs.fIdx, signal.depthIdx));
                    signal.markerY(2)   = squeeze(this.data.fitModelMat(signal.ch, this.figs.fIdx, signal.depthIdx));
                case 'fittedFFTSingleCh'
                    signal.yData     = this.data.fittedQAvgFFT(signal.ch, signal.fidxs, signal.depthIdx);
                    signal.titleVars = {{signal.ch, signal.depthVec(signal.depthIdx), signal.depthIdx}};
                    signal.markerY   = squeeze(this.data.fittedQAvgFFT(signal.ch, this.figs.fIdx, signal.depthIdx));
                case 'fittedFFTAllCh'
                    signal.yData     = this.data.fittedQAvgFFT(:, signal.fidxs, signal.depthIdx);
                    signal.markerY   = squeeze(this.data.fittedQAvgFFT(:,this.figs.fIdx, signal.depthIdx));
                    signal.yData     = squeeze(signal.yData(:,signal.fidxs));
                    signal.titleVars = {{signal.depthVec(signal.depthIdx), signal.depthIdx}};
                case 'avgFFT'
                    signal.yData     = this.data.fittedChAvgFFT(signal.fidxs, signal.depthIdx)';
                    signal.titleVars = {{signal.depthVec(signal.depthIdx), signal.depthIdx}};
                    signal.markerY   = squeeze(this.data.fittedChAvgFFT(this.figs.fIdx, signal.depthIdx));
                case 'normFFT'
                    signal.yData     = this.data.fittedFFTNorm(signal.fidxs, signal.depthIdx)';
                    signal.titleVars = {{signal.depthVec(signal.depthIdx), signal.depthIdx}};
                    signal.markerY   = squeeze(this.data.fittedFFTNorm(this.figs.fIdx, signal.depthIdx));
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
                        'YData', gather(yData));
                case 'plot'
                    set(this.figs.(gName).handles.cur.plot(hIdx), 'XData',...
                        xData,...
                        'YData', gather(yData));
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
            
            signal = this.extractSignal(figName);
            SNRstr = sprintf("SNR: %.2f\nBkg STD: %.2E", signal.SNR, signal.noiseStd);
            minVal = min(signal.yData); maxVal = max(signal.yData);
            span = maxVal-minVal;
            minLim = minVal - 0.1*span;
            maxLim = maxVal + 0.1*span;
            SNRYlim = double(maxLim - 0.05*span);
            
            if ~isgraphics(this.figs.(figName).handles.cur.plot) ||...
                ~strcmp(this.figs.(figName).type, this.figs.(figName).handles.cur.plot.Type)  
                switch this.figs.(figName).type
                    case 'stem'
                       this.figs.(figName).handles.cur.plot = ...
                           stem(this.figs.(figName).handles.cur.ax,...
                           signal.xData, signal.yData);    
                    case 'plot'
                       this.figs.(figName).handles.cur.plot = ...
                           plot(this.figs.(figName).handles.cur.ax,...
                           signal.xData, signal.yData);
                   case 'errorbar'
                       this.figs.(figName).handles.cur.plot = ...
                           errorbar(this.figs.(figName).handles.cur.ax,...
                           signal.xData, signal.yData, signal.std);
                end
                this.figs.(figName).handles.cur.text = ...
                    text(this.figs.(figName).handles.cur.ax, 0, SNRYlim,...
                    SNRstr, 'FontSize', 14);
                setStringsToPlot(this, figName) 
                ylim(this.figs.(figName).handles.cur.ax, [minLim, maxLim]);
            else
                quickPlot(this, figName, signal.xData, signal.yData, [], [], []);
                set(this.figs.(figName).handles.cur.text, 'String', SNRstr, 'Position', [0, SNRYlim]);
                pause(0.01);
                ylim(this.figs.(figName).handles.cur.ax, [minLim, maxLim]);
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
