classdef AOGraphics < Graphics
    %AOGRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data
        usSignal
        extClkSignal
        tUS
        tExtClk
        tRawData
        tMeasSig
        tNetSig
        tPos
        zAxis
    end
    
    methods (Static)
        function gNames = getGraphicsNames()
            gNames = {'extClk';'usSignal';'fullSignal';'measSamples';'netSignal'; 'deMul'; 'reshapedSignal';'FFT'; 'qAvgFFT'; 'phiCh';'phi'; 'deMul'};
        end
        
        function figs  = createGraphicsVars()
            figs.zIdx  = 1;
            figs.ch    = 1;
            figs.quant = 1;
            
            figs.fBar      = [];
            figs.fBarShift = [];
            figs.fIdx      = [];
            figs.zVec      = [];
            figs.displayFullFFT = false;
            figs.FFTenv = [];
            figs.FFTIdxs = [];
            
            figs.intExt  = 'int';
            figsNames = AOGraphics.getGraphicsNames();
            
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

                figs.(figsNames{i}).dims.zDim  = [];
                figs.(figsNames{i}).dims.chDim   = [];
                figs.(figsNames{i}).dims.dataDim = [];

                figs.(figsNames{i}).markers.plotMark = false;

                figs.(figsNames{i}).lims.xlims  = [];
                figs.(figsNames{i}).lims.ylims  = [];
                figs.(figsNames{i}).lims.clims = [];

                figs.(figsNames{i}).fonts.type       = [];
                figs.(figsNames{i}).fonts.titleSize  = 18;
                figs.(figsNames{i}).fonts.labelsSize = 18;
                figs.(figsNames{i}).fonts.axisSize   = 18;
            end
        end 
        
        function uVars = createGraphicsUserVars() 
            figsNames = AOGraphics.getGraphicsNames();
            
            uVars.ch    = 1;
            uVars.zIdx  = 1;
            uVars.quant = 1;
             
            uVars.fBar = [];
            uVars.fBarShift = [];
            uVars.fIdx = [];
            uVars.zVec = [];
            uVars.displayFullFFT = false;
            uVars.FFTenv = [];
            
            uVars.intExt = [];
            for i=1:length(figsNames)
                uVars.validStruct.(figsNames{i}) = false;
            end
            
            for i=1:length(figsNames)
                uVars.extH.(figsNames{i}) =  Graphics.createHandlesStruct();
            end
            
            uVars.fonts.type       = [];
            uVars.fonts.titleSize  = 14;
            uVars.fonts.labelsSize = 14;
            uVars.fonts.axisSize   = 14;
        end
        
    end
    
    methods
        function this = AOGraphics()
            this@Graphics()
            
            this.figsNames = this.getGraphicsNames();
            this.numOfFigs = length(this.figsNames);
            
            this.figs     = AOGraphics.createGraphicsVars();
            
            this.uVars    = AOGraphics.createGraphicsUserVars();
            
            this.setGraphicsStaticVars();
        end
        
        % Set Functions
        
        function setFrequencyBar(this, fBar, fBarShift, fIdx)
            this.figs.fBar      = fBar;
            this.figs.fBarShift = fBarShift;
            this.figs.fIdx      = fIdx;            
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
           
            % fullSignal [ch x samplesPerAcq]
            this.setType('fullSignal', 'plot');
            this.setStrings('fullSignal', "Full Signal Channel: %d", "t[\mus]", "Amp[V]", []);

            % measSamples [ch x samplesPerMeas]
            this.setType('measSamples', 'plot');
            this.setStrings('measSamples', "Measured Signal Channel: %d", "t[\mus]", "Amp[V]", []);
            
            % netSignal [ch x samplesPerSignal]
            this.setType('netSignal', 'plot');
            this.setStrings('netSignal', "Net Signal Channel: %d", "t[\mus]", "Amp[V]", []);
            
            this.setType('deMul', 'plot');
            this.setStrings('deMul', "Demultiplexed Signal Channel: %d", "t[\mus]", "Amp[V]", []);
            
            % reshapedSignal [ch x samplesPerPos x numOfPos]
            this.setType('reshapedSignal', 'stem');
            this.setStrings('reshapedSignal', "Reshaped Signal Channel: %d, Z: %.2f[mm] (idx:(%d))", "t[\mu s]", "Amp[V]", []);
                        
            %FFT [ch x samplesPerPos x numOfPos]
            this.setType('FFT', 'plot');
            this.setStrings('FFT', "FFT Z: %.2f[mm] (idx:(%d)), Quant: %d", "f[MHz]", "Power Spectrum", "Ch %d");
            this.setMarkersEnable('FFT', true);
            
            this.setType('qAvgFFT', 'plot');
            this.setStrings('qAvgFFT', "Quant AVG FFT Z: %.2f[mm] (idx: %d)", "f[MHz]", "Power Spectrum", "Ch %d");
            this.setMarkersEnable('qAvgFFT', true);
            
            % phiCh [ch x numOfPos]
            this.setType('phiCh', 'stem');
            this.setStrings('phiCh', "\\phi_{ch}", "z[mm]", "\phi", "Ch %d");
            
            % phi [numOfPos]
            this.setType('phi', 'stem');
            this.setStrings('phi', "\\phi", "z[mm]", "\phi", []);
        end
        
        function setGraphicsDynamicVars(this, algoVars)
            %this function assums that uVars are updated
            this.setChAndPos(this.uVars.ch, this.uVars.zIdx, this.uVars.quant)
            
%             z    = algoVars.len.zVecUSRes(this.uVars.zIdx)*1e3;
            ch   = algoVars.digitizer.channels;
%             fSin = algoVars.usSignal.fSin;
            
            this.figs.displayFullFFT = this.uVars.displayFullFFT;
            this.figs.FFTenv = this.uVars.FFTenv;
            envIdx = floor(this.uVars.FFTenv / algoVars.freq.df);
            this.figs.FFTIdxs = algoVars.freq.fSinIdx + [-envIdx, envIdx];
            
            this.setFrequencyBar(algoVars.freq.frequencyBar*1e-6, algoVars.freq.frequencyBarShifted*1e-6, algoVars.freq.fSinIdx);
            this.setZVec(algoVars.len.zVecUSRes);
            
            this.setLegendVariables('FFT', 1:ch);
            this.setLegendVariables('qAvgFFT', 1:ch);
            this.setLegendVariables('phiCh', 1:ch);
            
            this.figs.fonts = this.uVars.fonts;
            
            this.updateGraphicsConstruction()
            if algoVars.uVars.useQuant
                this.setType('phi', 'stem');
            else
                this.setType('phi', 'stem');
            end
            
            this.tUS         = algoVars.timing.tUS*1e6;
            this.tExtClk     = algoVars.timing.tExtClk*1e6;
            this.tRawData    = algoVars.timing.tAcqVec*1e6;
            this.tMeasSig    = algoVars.timing.tMeasVec*1e6;
            this.tNetSig     = algoVars.timing.tSigVec*1e6;
            this.tPos        = algoVars.timing.tPosVec*1e6;
            this.zAxis       = algoVars.len.zVecUSRes*1e3;
        end
        
        function setZVec(this, zVec)
            this.figs.zVec = zVec;
        end
        
        function setData(this, data, usSig, extClkSig)
            this.data = data;
            if nargin > 2
                this.usSig = usSig;
                this.extClkSig = extClkSig;
            end
        end
         
        function signal = extractSignal(this, figName)
            switch figName
                case 'extClk'
                    signal.xData = this.tExtClk;
                    signal.yData = this.extClkSignal;
                case 'usSignal'
                    signal.xData = this.tUS;
                    signal.yData = this.usSignal;
                case 'fullSignal'
                    signal.ch  = this.figs.ch;
                    signal.xData = this.tRawData;
                    signal.yData = this.data.rawData(signal.ch, :);
                case 'measSamples'
                    signal.ch    = this.figs.ch;
                    signal.xData = this.tMeasSig;
                    signal.yData = this.data.measSignal(signal.ch, :);
                case 'netSignal'
                    signal.ch    = this.figs.ch;
                    signal.xData = this.tNetSig;
                    signal.yData = this.data.netSignal(signal.ch, :);
                case 'deMul'
                    signal.ch    = this.figs.ch;
                    signal.xData = this.tNetSig;
                    signal.yData = this.data.deMultiplexed(signal.ch, :); 
                case 'reshapedSignal'
                    signal.ch    = this.figs.ch;
                    signal.zIdx  = this.figs.zIdx;
                    signal.quant = this.figs.quant;
                    signal.xData = this.tPos;
                    signal.yData = permute( this.data.reshapedSignal(signal.quant ,signal.ch, :, signal.zIdx), [2,3,1,4]);
                case 'phiCh'
                    signal.quant = this.figs.quant;
                    signal.xData = this.zAxis;
                    signal.yData = this.data.phiCh(signal.quant, :, :);
                case 'phi'
                    signal.xData = this.zAxis;
                    signal.yData = this.data.phi;
                    signal.std   = this.data.phiStd;
            end
        end
        
        function signal = extractFFTSignal(this, figName)
            signal.zIdx  = this.figs.zIdx;
            signal.zVec  = this.figs.zVec;
            switch figName
                case 'FFT'
                    signal.yData = abs(this.data.fftRes);
                    signal.quant = this.figs.quant;
                    signal.titleVars = {{signal.zVec(signal.zIdx)*1e3, signal.zIdx, signal.quant}};
                case 'qAvgFFT'
                    signal.yData = abs(this.data.qAvgFFT);
                    signal.quant = 1;
                    signal.titleVars = {{signal.zVec(signal.zIdx)*1e3, signal.zIdx}};
            end
            
            if ~this.figs.displayFullFFT
                % plot user defined enviroment around the signal
                signal.fidxs = this.figs.FFTIdxs(1):this.figs.FFTIdxs(2);
                this.figs.(figName).lims.xlims = ...
                    [this.figs.fBar(this.figs.FFTIdxs(1)), this.figs.fBar(this.figs.FFTIdxs(2))];
                signal.xData = this.figs.fBar(signal.fidxs);
            else
                signal.fidxs = 1:length(this.figs.fBar);
                this.figs.(figName).lims.xlims = ...
                    [this.figs.fBarShift(1), this.figs.fBarShift(end)];
                signal.xData = this.figs.fBarShift;
                signal.yData = fftshift(signal.yData, 3);
            end
                     
            signal.markerX = this.figs.fBar(this.figs.fIdx);
            signal.markerY = squeeze(signal.yData(signal.quant, :, this.figs.fIdx, signal.zIdx));
            
            signal.yData = squeeze(signal.yData(signal.quant, :, signal.fidxs, signal.zIdx));
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
            if ~isgraphics(this.figs.(figName).handles.cur.ax)
               return
            end
            
            signal = this.extractSingal(figName);
            
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
                
        function plotPhi(this)
            %yData - [samples per pos, 1]
            
            if ~isgraphics(this.figs.phi.handles.cur.ax)
               return
            end
            
            signal = this.extractSignal('phi');
            
            if ~isgraphics(this.figs.phi.handles.cur.plot) ||...
                ~strcmp(this.figs.phi.type, this.figs.phi.handles.cur.plot.Type)  
                switch this.figs.phi.type
                    case 'stem'
                       this.figs.phi.handles.cur.plot = ...
                           stem(this.figs.phi.handles.cur.ax,...
                           signal.xData, signal.yData);    
                    case 'plot'
                       this.figs.phi.handles.cur.plot = ...
                           plot(this.figs.phi.handles.cur.ax,...
                           signal.xData, signal.yData);
                   case 'errorbar'
                       this.figs.phi.handles.cur.plot = ...
                           errorbar(this.figs.phi.handles.cur.ax,...
                           signal.xData, signal.yData, signal.std);
                end
                setStringsToPlot(this, 'phi')
            else
                quickPlot(this, 'phi', signal.xData, signal.yData, signal.std, [], []);
                pause(0.01);
                drawnow();
            end
        end
        
        function plotPhiCh(this)
            %yData - [ch x numOfPos]
            if ~isgraphics(this.figs.phiCh.handles.cur.ax)
               return
            end
            
            signal = this.extractSignal('phiCh');
                        
            if ~isgraphics(this.figs.phiCh.handles.cur.plot)
                hold(this.figs.phiCh.handles.cur.ax, 'on');
                for i = 1:size(signal.yData, 2)
                    switch this.figs.phiCh.type
                        case 'stem'
                           this.figs.phiCh.handles.cur.plot(i) = ...
                               stem(this.figs.phiCh.handles.cur.ax,...
                               signal.xData, squeeze(signal.yData(:, i, :)));    
                        case 'plot'
                           this.figs.(gName).handles.cur.plot(i) = ...
                               plot(this.figs.(gName).handles.cur.ax,...
                               signal.xData, squeeze(signal.yData(:, i, :)));  
                    end
                end
%                 hold(this.figs.(gName).handles.cur.ax, 'off');
                this.setStringsToPlot('phiCh');
            else
                for i = 1:size(this.data.phiCh, 2)
                    quickPlot(this, 'phiCh', signal.xData, squeeze(signal.yData(:, i, :)),[], [], [], i)
                end
                pause(0.01);
                drawnow();
            end
        end
        
        function plotSignal(this, figName)
        % yData - [ch x Sig]
            
            if ~isgraphics(this.figs.(figName).handles.cur.ax)
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
            
            if ~isgraphics(this.figs.reshapedSignal.handles.cur.ax)
               return
            end
            
            signal = this.extractSignal('reshapedSignal');
            
            if ~isgraphics(this.figs.reshapedSignal.handles.cur.plot)
                switch this.figs.reshapedSignal.type
                    case 'stem'
                       this.figs.reshapedSignal.handles.cur.plot = ...
                           stem(this.figs.reshapedSignal.handles.cur.ax,...
                           signal.xData, signal.yData);    
                    case 'plot'
                       this.figs.reshapedSignal.handles.cur.plot = ...
                           plot(this.figs.reshapedSignal.handles.cur.ax,...
                           signal.xData, signal.yData);  
                end
                this.setStringsToPlot('reshapedSignal');
            else
                quickPlot(this, 'reshapedSignal', signal.xData, signal.yData,[], [], [])
                pause(0.01);
                drawnow()
            end
        end
        
        function plotFFT(this, figName)
            if ~isgraphics(this.figs.(figName).handles.cur.ax)
               return
            end
            
            signal = this.extractFFTSignal(figName);
            this.setTitleVariables(figName, signal.titleVars);          % FFT
 
            if ~isgraphics(this.figs.(figName).handles.cur.plot(1))
            	hold(this.figs.(figName).handles.cur.ax, 'on');
                for i = 1:size(signal.yData, 1)
                    switch this.figs.(figName).type
                        case 'stem'
                           this.figs.(figName).handles.cur.plot(i) = ...
                               stem(this.figs.(figName).handles.cur.ax,...
                               signal.xData, signal.yData( i, :, :));    
                        case 'plot'
                           this.figs.(figName).handles.cur.plot(i) = ...
                               plot(this.figs.(figName).handles.cur.ax,...
                               signal.xData, signal.yData( i, :, :));  
                    end
                end
                
                for i = 1:length(signal.markerY)            
                    this.plotMarkers(figName, signal.markerX, signal.markerY(i), i)
                end
                hold(this.figs.(figName).handles.cur.ax, 'off');
                
                this.setStringsToPlot(figName)
                this.setLimsToPlot(figName)
            else
                for i = 1:size(signal.yData, 1)
                    this.quickPlot(figName, signal.xData, signal.yData(i, :, :),[], signal.markerX, signal.markerY(i), i)
                end
                this.setLimsToPlot(figName)
                pause(0.01);
                drawnow()
            end
        end
        
        function plotMarkers(this, gName, xMarker, yMarker, hIdx)
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
