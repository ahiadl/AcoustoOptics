classdef AOGraphics < Graphics
    %AOGRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
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
            
            figs.intExt  = [];
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
            
            this.uVarsOld = AOGraphics.createGraphicsUserVars();
            this.uVars    = AOGraphics.createGraphicsUserVars();
            
            this.setGraphicsStaticVars();
        end
        
        % Set Functions
        
        function setDims(this, gName, zDim, chDim, dataDim)
            this.figs.(gName).dims.zDim  = zDim;
            this.figs.(gName).dims.chDim   = chDim;
            this.figs.(gName).dims.dataDim = dataDim;
        end
        
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
            
%             setDims - when dim is irrelevant send 1
            % extClk
            this.setType('extClk', 'stem');
            this.setStrings('extClk', {"External Clk (Sampling Clock)"}, "t[\mus]", "Amp", []);
            this.setDims('extClk', [], [], 1);
            
            % usSignal
            this.setType('usSignal', 'stem');
            this.setStrings('usSignal', {"Ultrasound Pulse"}, "t[\mus]", "Amp", []);
            this.setDims('usSignal', [], [], 1);
           
            % fullSignal [ch x samplesPerAcq]
            this.setType('fullSignal', 'plot');
            this.setStrings('fullSignal', {"Full Signal Channel: %d"}, "t[\mus]", "Amp[V]", []);
            this.setDims('fullSignal', [], 1, 2);

            % measSamples [ch x samplesPerMeas]
            this.setType('measSamples', 'plot');
            this.setStrings('measSamples', {"Measured Signal Channel: %d"}, "t[\mus]", "Amp[V]", []);
            this.setDims('measSamples', [], 1, 2);
            
            % netSignal [ch x samplesPerSignal]
            this.setType('netSignal', 'plot');
            this.setStrings('netSignal', {"Net Signal Channel: %d"}, "t[\mus]", "Amp[V]", []);
            this.setDims('netSignal', [], 1, 2);
            
            this.setType('deMul', 'plot');
            this.setStrings('deMul', {"Demultiplexed Signal Channel: %d"}, "t[\mus]", "Amp[V]", []);
            this.setDims('deMul', [], 1, 2);
            
            % reshapedSignal [ch x samplesPerPos x numOfPos]
            this.setType('reshapedSignal', 'stem');
            this.setStrings('reshapedSignal', {"Reshaped Signal Channel: %d, Z: %.2f[mm] (idx:(%d))"}, "t[\mu s]", "Amp[V]", []);
            this.setDims('reshapedSignal', 3, 1, 2);
                        
            %FFT [ch x samplesPerPos x numOfPos]
            this.setType('FFT', 'plot');
            this.setStrings('FFT', {"FFT Z: %.2f[mm] (idx:(%d)), Quant: %d"}, "f[MHz]", "Power Spectrum", "Ch %d");
            this.setDims('FFT', 3, 1, 2);
            this.setMarkersEnable('FFT', true);
            
            this.setType('qAvgFFT', 'plot');
            this.setStrings('qAvgFFT', {"Quant AVG FFT Z: %.2f[mm] (idx:(%d))"}, "f[MHz]", "Power Spectrum", "Ch %d");
            this.setDims('qAvgFFT', 3, 1, 2);
            this.setMarkersEnable('qAvgFFT', true);
            
            % phiCh [ch x numOfPos]
            this.setType('phiCh', 'stem');
            this.setStrings('phiCh', {"\\phi_{ch}"}, "z[mm]", "\phi", "Ch %d");
            this.setDims('phiCh', 2, 1, 2);
            
            % phi [numOfPos]
            this.setType('phi', 'stem');
            this.setStrings('phi', {"\\phi"}, "z[mm]", "\phi", []);
            this.setDims('phi', 1, [], 1);
        end
        
        function setGraphicsDynamicVars(this, algoVars)
            %this function assums that uVars are updated
            this.setChAndPos(this.uVars.ch, this.uVars.zIdx, this.uVars.quant)
            
            z    = algoVars.len.zVecUSRes(this.uVars.zIdx)*1e3;
            ch   = algoVars.digitizer.channels;
            fSin = algoVars.usSignal.fSin;
            
            this.figs.displayFullFFT = this.uVars.displayFullFFT;
            this.figs.FFTenv = this.uVars.FFTenv;
            envIdx = floor(this.uVars.FFTenv / algoVars.freq.df);
            this.figs.FFTIdxs = algoVars.freq.fSinIdx + [-envIdx, envIdx];
            
            this.setFrequencyBar(algoVars.freq.frequencyBar*1e-6, algoVars.freq.frequencyBarShifted*1e-6, algoVars.freq.fSinIdx);
            this.setZVec(algoVars.len.zVecUSRes);
%             this.setTitleVariables('extClk',      {[]});
%             this.setTitleVariables('usSignal',    {[]});
%             this.setTitleVariables('fullSignal',  {[this.uVars.ch]});              % fullSignal 
%             this.setTitleVariables('measSamples', {[this.uVars.ch]});              % measSamples
%             this.setTitleVariables('netSignal',   {[this.uVars.ch]});
%             this.setTitleVariables('deMul',       {[this.uVars.ch]});% netSignal
%             this.setTitleVariables('reshapedSignal', {[this.uVars.ch, z, this.uVars.zIdx]}); % reshapedSignal
%             this.setTitleVariables('FFT', {[z, this.uVars.zIdx]});          % FFT
%             this.setTitleVariables('qAvgFFT', {[z, this.uVars.zIdx]});          % FFT
%             this.setTitleVariables('phiCh', {[]});
%             this.setTitleVariables('phi', {[]});
            
            this.setLegendVariables('FFT', 1:ch);
            this.setLegendVariables('qAvgFFT', 1:ch);
            this.setLegendVariables('phiCh', 1:ch);
             
%             this.setLimits('FFT', (fSin + fSin*[-0.001, 0.001])*1e-6, []);
%             this.setLimits('FFT', (fSin + fSin*[-0.001, 0.001])*1e-6, []);
            
            this.updateGraphicsConstruction()
            if algoVars.uVars.useQuant
                this.setType('phi', 'stem');
            else
                this.setType('phi', 'stem');
            end
        end
        
        function setZVec(this, zVec)
            this.figs.zVec = zVec;
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
        end
        
        function plotAFGSignals(this, gName, xData, yData)
            %yData - [1 x Sig]
            if ~isgraphics(this.figs.(gName).handles.cur.ax)
               return
            end
            
            if ~isgraphics(this.figs.(gName).handles.cur.plot)
                switch this.figs.(gName).type
                    case 'stem'
                       this.figs.(gName).handles.cur.plot = ...
                           stem(this.figs.(gName).handles.cur.ax,...
                           xData, yData);    
                    case 'plot'
                       this.figs.(gName).handles.cur.plot = ...
                           plot(this.figs.(gName).handles.cur.ax,...
                           xData, yData);  
                end
                setStringsToPlot(this, gName)
            else
                quickPlot(this, gName, xData, yData,[], [], [])
                drawnow();
            end
        end
                
        function plotPhi(this, gName, xData, yData, std)
            %yData - [samples per pos, 1]
            
            if ~isgraphics(this.figs.(gName).handles.cur.ax)
               return
            end
            if ~isgraphics(this.figs.(gName).handles.cur.plot) ||...
                ~strcmp(this.figs.(gName).type, this.figs.(gName).handles.cur.plot.Type)  
                switch this.figs.(gName).type
                    case 'stem'
                       this.figs.(gName).handles.cur.plot = ...
                           stem(this.figs.(gName).handles.cur.ax,...
                           xData, yData);    
                    case 'plot'
                       this.figs.(gName).handles.cur.plot = ...
                           plot(this.figs.(gName).handles.cur.ax,...
                           xData, yData);
                   case 'errorbar'
                       this.figs.(gName).handles.cur.plot = ...
                           errorbar(this.figs.(gName).handles.cur.ax,...
                           xData, yData, std);
                end
                setStringsToPlot(this, gName)
            else
                quickPlot(this, gName, xData, gather(yData), gather(std), [], []);
                drawnow();
            end
        end
        
        function plotPhiCh(this, gName, xData, yData)
            %yData - [ch x numOfPos]
            if ~isgraphics(this.figs.(gName).handles.cur.ax)
               return
            end
            
            quant = this.figs.quant;
                        
            if ~isgraphics(this.figs.(gName).handles.cur.plot)
                hold(this.figs.(gName).handles.cur.ax, 'on');
                for i = 1:size(yData, 2)
                    switch this.figs.(gName).type
                        case 'stem'
                           this.figs.(gName).handles.cur.plot(i) = ...
                               stem(this.figs.(gName).handles.cur.ax,...
                               xData, squeeze(yData(quant, i, :)));    
                        case 'plot'
                           this.figs.(gName).handles.cur.plot(i) = ...
                               plot(this.figs.(gName).handles.cur.ax,...
                               xData, squeeze(yData(quant, i, :)));  
                    end
                end
%                 hold(this.figs.(gName).handles.cur.ax, 'off');
                this.setStringsToPlot(gName);
            else
                for i = 1:size(yData, 2)
                    quickPlot(this, gName, xData, squeeze(yData(quant, i, :)),[], [], [], i)
                end
                drawnow();
            end
        end
        
        function plotSignal(this, gName, xData, yData)
        % yData - [ch x Sig]
            
            if ~isgraphics(this.figs.(gName).handles.cur.ax)
               return
            end
     
            ch  = this.figs.ch;

            if ~isgraphics(this.figs.(gName).handles.cur.plot)
                switch this.figs.(gName).type
                    case 'stem'
                       this.figs.(gName).handles.cur.plot = ...
                           stem(this.figs.(gName).handles.cur.ax,...
                           xData, yData(ch, :));    
                    case 'plot'
                       this.figs.(gName).handles.cur.plot = ...
                           plot(this.figs.(gName).handles.cur.ax,...
                           xData, yData(ch, :));  
                end
                this.setStringsToPlot(gName);
            else
                quickPlot(this, gName, xData, gather(yData(ch, :)),[], [], [])
                drawnow();
            end
        end
        
        function plotReshapedSignal(this, gName, xData, yData)
        % yData - [ch x Sig x pos]
            
            if ~isgraphics(this.figs.(gName).handles.cur.ax)
               return
            end
     
            ch  = this.figs.ch;
            zIdx = this.figs.zIdx;
            quant = this.figs.quant;
            
            yData = permute( yData(quant ,ch, :, zIdx), [2,3,1,4]);
            
            if ~isgraphics(this.figs.(gName).handles.cur.plot)
                switch this.figs.(gName).type
                    case 'stem'
                       this.figs.(gName).handles.cur.plot = ...
                           stem(this.figs.(gName).handles.cur.ax,...
                           xData, yData);    
                    case 'plot'
                       this.figs.(gName).handles.cur.plot = ...
                           plot(this.figs.(gName).handles.cur.ax,...
                           xData, yData);  
                end
                this.setStringsToPlot(gName);
            else
                quickPlot(this, gName, xData, gather(yData),[], [], [])
                drawnow()
            end
        end
        
        function plotFFT(this, gName, xData, yData, markerX, markerYIdx)
            if ~isgraphics(this.figs.(gName).handles.cur.ax)
               return
            end
            
            zIdx  = this.figs.zIdx;
            if ~this.figs.displayFullFFT
                % plot user defined enviroment around the signal
                fidxs = this.figs.FFTIdxs(1):this.figs.FFTIdxs(2);
                this.figs.(gName).lims.xlims = ...
                    [this.figs.fBar(this.figs.FFTIdxs(1)), this.figs.fBar(this.figs.FFTIdxs(2))];
                xDataReduced = this.figs.fBar(fidxs);
            else
                fidxs = 1:length(xData);
                this.figs.(gName).lims.xlims = ...
                    [this.figs.fBarShift(1), this.figs.fBarShift(end)];
                xDataReduced = this.figs.fBarShift;
                yData = fftshift(yData, 3);
            end
            
            zVec  = this.figs.zVec;
            quant = this.figs.quant;
            
%             xDataReduced = xData(fidxs);
            yDataReduced = squeeze(yData(quant, :, fidxs, zIdx));
            
            this.setTitleVariables(gName, {[zVec(zIdx)*1e3, zIdx, quant]});          % FFT
            
            markerY = squeeze(yData(quant, :, this.figs.fIdx, zIdx));
            
            if ~isgraphics(this.figs.(gName).handles.cur.plot(1))
            	hold(this.figs.(gName).handles.cur.ax, 'on');
                for i = 1:size(yDataReduced, 1)
                    switch this.figs.(gName).type
                        case 'stem'
                           this.figs.(gName).handles.cur.plot(i) = ...
                               stem(this.figs.(gName).handles.cur.ax,...
                               xDataReduced, yDataReduced( i, :, :));    
                        case 'plot'
                           this.figs.(gName).handles.cur.plot(i) = ...
                               plot(this.figs.(gName).handles.cur.ax,...
                               xDataReduced, yDataReduced( i, :, :));  
                    end
                end
                
                for i = 1:length(markerY)            
                    this.plotMarkers(gName, markerX, markerY(i), i)
                end
                hold(this.figs.(gName).handles.cur.ax, 'off');
                
                this.setStringsToPlot(gName)
                this.setLimsToPlot(gName)
            else

                for i = 1:size(yData, 2)
                    this.quickPlot(gName, xDataReduced, yDataReduced(i, :, :),[], markerX, markerY(i), i)
%                     this.quickPlot(gName, xData(fidxs), yData(i, fidxs, zIdx), markerX, markerY(i, zIdx), i)
                end
                this.setLimsToPlot(gName)
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
