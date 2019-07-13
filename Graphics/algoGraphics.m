classdef algoGraphics < Graphics
    %AOGRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function graphReq = getGraphicsRequest(plotNames)
            graphReq.zIdx  = 1;
            graphReq.ch    = 1;
            graphReq.quant = 1;
            
            graphReq.fBar = [];
            graphReq.fIdx = [];
            graphReq.zVec = [];
            
            graphReq.intExt  = [];
            
            for i = 1:length(plotNames)  
                graphReq.(plotNames{i}).type      = []; 
                
                graphReq.(plotNames{i}).handles.int = Graphics.createHandlesStruct();
                graphReq.(plotNames{i}).handles.ext = Graphics.createHandlesStruct();
                graphReq.(plotNames{i}).handles.cur = Graphics.createHandlesStruct();
                
                graphReq.(plotNames{i}).strings.title       = [];
                graphReq.(plotNames{i}).strings.titleModel  = [];
                graphReq.(plotNames{i}).strings.xlabel      = [];
                graphReq.(plotNames{i}).strings.ylabel      = [];
                graphReq.(plotNames{i}).strings.legend      = [];
                graphReq.(plotNames{i}).strings.legendModel = [];

                graphReq.(plotNames{i}).dims.zDim  = [];
                graphReq.(plotNames{i}).dims.chDim   = [];
                graphReq.(plotNames{i}).dims.dataDim = [];

                graphReq.(plotNames{i}).markers.plotMark = false;

                graphReq.(plotNames{i}).lims.xlims  = [];
                graphReq.(plotNames{i}).lims.ylims  = [];
                graphReq.(plotNames{i}).lims.clims = [];

                graphReq.(plotNames{i}).fonts.type       = [];
                graphReq.(plotNames{i}).fonts.titleSize  = 18;
                graphReq.(plotNames{i}).fonts.labelsSize = 18;
                graphReq.(plotNames{i}).fonts.axisSize   = 18;
            end
        end 
        
        function gReq = createGraphicsRunVars() %for Reduced Operation
            gNames = AlgoNew.getGraphicsNames();
            
            gReq.ch    = 1;
            gReq.zIdx  = 1;
            gReq.quant = 1;
             
            gReq.fBar = [];
            gReq.fIdx = [];
            gReq.zVec = [];
            
            gReq.intExt = [];
            for i=1:length(gNames)
                gReq.validStruct.(gNames{i}) = false;
            end
            
            for i=1:length(gNames)
                gReq.extH.(gNames{i}) =  Graphics.createHandlesStruct();
            end
            
            gReq.fonts.type       = [];
            gReq.fonts.titleSize  = 14;
            gReq.fonts.labelsSize = 14;
            gReq.fonts.axisSize   = 14;
        end
        
    end
    
    methods
        function this = algoGraphics(gNames)
            this@Graphics(gNames)
            this.requests = algoGraphics.getGraphicsRequest(gNames);
            this.globalReqOld = algoGraphics.createGraphicsRunVars();
            this.globalReq = algoGraphics.createGraphicsRunVars();
        end
        
        % Set Functions
        function setDims(this, gName, zDim, chDim, dataDim)
            this.requests.(gName).dims.zDim  = zDim;
            this.requests.(gName).dims.chDim   = chDim;
            this.requests.(gName).dims.dataDim = dataDim;
        end
        
        function setFrequencyBar(this, fBar, fIdx)
            this.requests.fBar = fBar;
            this.requests.fIdx = fIdx;            
        end
        
        function setGraphicsStaticVars(this)
%             this.graphics.setType(type);
%             this.graphics.setStrings(title, xlabel, ylabel, legend);
%             this.graphics.setDims(gName, zIdxDim, chDim, dataDim)
%             this.graphics.setMarkersEnable(enable);          
            
%             setDims - when dim is irrelevant send 1

            % extClk
            this.setType(this.graphNames{1}, 'stem');
            this.setStrings(this.graphNames{1}, {"External Clk (Sampling Clock)"}, "t[\mus]", "Amp", []);
            this.setDims(this.graphNames{1}, [], [], 1);
            
            % usSignal
            this.setType(this.graphNames{2}, 'stem');
            this.setStrings(this.graphNames{2}, {"Ultrasound Pulse"}, "t[\mus]", "Amp", []);
            this.setDims(this.graphNames{2}, [], [], 1);
           
            % fullSignal [ch x samplesPerAcq]
            this.setType(this.graphNames{3}, 'plot');
            this.setStrings(this.graphNames{3}, {"Full Signal Channel: %d"}, "t[\mus]", "Amp[V]", []);
            this.setDims(this.graphNames{3}, [], 1, 2);

            % measSamples [ch x samplesPerMeas]
            this.setType(this.graphNames{4}, 'plot');
            this.setStrings(this.graphNames{4}, {"Measured Signal Channel: %d"}, "t[\mus]", "Amp[V]", []);
            this.setDims(this.graphNames{4}, [], 1, 2);
            
            % netSignal [ch x samplesPerSignal]
            this.setType(this.graphNames{5}, 'plot');
            this.setStrings(this.graphNames{5}, {"Net Signal Channel: %d"}, "t[\mus]", "Amp[V]", []);
            this.setDims(this.graphNames{5}, [], 1, 2);

            % reshapedSignal [ch x samplesPerPos x numOfPos]
            this.setType(this.graphNames{6}, 'stem');
            this.setStrings(this.graphNames{6}, {"Reshaped Signal Channel: %d, Z: %.2f[mm] (idx:(%d))"}, "t[\mu s]", "Amp[V]", []);
            this.setDims(this.graphNames{6}, 3, 1, 2);
                        
            %FFT [ch x samplesPerPos x numOfPos]
            this.setType(this.graphNames{7}, 'plot');
            this.setStrings(this.graphNames{7}, {"FFT Z: %.2f[mm] (idx:(%d))"}, "f[MHz]", "Power Spectrum", "Ch %d");
            this.setDims(this.graphNames{7}, 3, 1, 2);
            this.setMarkersEnable(this.graphNames{7}, true);
            
            % phiCh [ch x numOfPos]
            this.setType(this.graphNames{8}, 'stem');
            this.setStrings(this.graphNames{8}, {"\\phi_{ch}"}, "z[mm]", "\phi", "Ch %d");
            this.setDims(this.graphNames{8}, 2, 1, 2);
            
            % phi [numOfPos]
            this.setType(this.graphNames{9}, 'stem');
            this.setStrings(this.graphNames{9}, {"\\phi"}, "z[mm]", "\phi", []);
            this.setDims(this.graphNames{1}, 1, [], 1);
        end
        
        function setZVec(this, zVec)
            this.requests.zVec = zVec;
        end
        
        % Disply Functions
        
        function quickPlot(this, gName, xData, yData, markerX, markerY, hIdx) 
            if ~exist('hIdx', 'var')
                hIdx = 1;
            end
            set(this.requests.(gName).handles.cur.plot(hIdx), 'XData', xData, 'YData', gather(yData));
            if this.requests.(gName).markers.plotMark
                set(this.requests.(gName).handles.cur.plotMarker(hIdx), 'XData', markerX, 'YData', gather(markerY));
            end
            set(this.requests.(gName).handles.cur.title, 'String', this.requests.(gName).strings.title);
%             tic
%             drawnow();
%             toc
        end
        
        function plotAFGSignals(this, gName, xData, yData)
            %yData - [1 x Sig]
            if ~isgraphics(this.requests.(gName).handles.cur.ax)
               return
            end
            
            if ~isgraphics(this.requests.(gName).handles.cur.plot)
                switch this.requests.(gName).type
                    case 'stem'
                       this.requests.(gName).handles.cur.plot = ...
                           stem(this.requests.(gName).handles.cur.ax,...
                           xData, yData);    
                    case 'plot'
                       this.requests.(gName).handles.cur.plot = ...
                           plot(this.requests.(gName).handles.cur.ax,...
                           xData, yData);  
                end
                setStringsToPlot(this, gName)
            else
                quickPlot(this, gName, xData, yData, [], [])
                drawnow();
            end
        end
                
        function plotPhi(this, gName, xData, yData)
            %yData - [samples per pos, 1]
            
            if ~isgraphics(this.requests.(gName).handles.cur.ax)
               return
            end
            if ~isgraphics(this.requests.(gName).handles.cur.plot)
                switch this.requests.(gName).type
                    case 'stem'
                       this.requests.(gName).handles.cur.plot = ...
                           stem(this.requests.(gName).handles.cur.ax,...
                           xData, yData);    
                    case 'plot'
                       this.requests.(gName).handles.cur.plot = ...
                           plot(this.requests.(gName).handles.cur.ax,...
                           xData, yData);  
                end
                setStringsToPlot(this, gName)
            else
                quickPlot(this, gName, xData, yData, [], []);
                drawnow();
            end
        end
        
        function plotPhiCh(this, gName, xData, yData)
            %yData - [ch x numOfPos]
            if ~isgraphics(this.requests.(gName).handles.cur.ax)
               return
            end
            
            if ~isgraphics(this.requests.(gName).handles.cur.plot)
                hold(this.requests.(gName).handles.cur.ax, 'on');
                for i = 1:size(yData, 1)
                    switch this.requests.(gName).type
                        case 'stem'
                           this.requests.(gName).handles.cur.plot(i) = ...
                               stem(this.requests.(gName).handles.cur.ax,...
                               xData, yData(i, :));    
                        case 'plot'
                           this.requests.(gName).handles.cur.plot(i) = ...
                               plot(this.requests.(gName).handles.cur.ax,...
                               xData, yData(i, :));  
                    end
                end
                hold(this.requests.(gName).handles.cur.ax, 'off');
                this.setStringsToPlot(gName);
            else
                for i = 1:size(yData, 1)
                    quickPlot(this, gName, xData, yData(i, :), [], [])
                end
                drawnow();
            end
        end
        
        function plotSignal(this, gName, xData, yData)
        % yData - [ch x Sig]
            
            if ~isgraphics(this.requests.(gName).handles.cur.ax)
               return
            end
     
            ch  = this.requests.ch;

            if ~isgraphics(this.requests.(gName).handles.cur.plot)
                switch this.requests.(gName).type
                    case 'stem'
                       this.requests.(gName).handles.cur.plot = ...
                           stem(this.requests.(gName).handles.cur.ax,...
                           xData, yData(ch, :));    
                    case 'plot'
                       this.requests.(gName).handles.cur.plot = ...
                           plot(this.requests.(gName).handles.cur.ax,...
                           xData, yData(ch, :));  
                end
                this.setStringsToPlot(gName);
            else
                quickPlot(this, gName, xData, gather(yData(ch, :)), [], [])
                drawnow();
            end
        end
        
        function plotReshapedSignal(this, gName, xData, yData)
        % yData - [ch x Sig x pos]
            
            if ~isgraphics(this.requests.(gName).handles.cur.ax)
               return
            end
     
            ch  = this.requests.ch;
            zIdx = this.requests.zIdx;
            
            if ~isgraphics(this.requests.(gName).handles.cur.plot)
                switch this.requests.(gName).type
                    case 'stem'
                       this.requests.(gName).handles.cur.plot = ...
                           stem(this.requests.(gName).handles.cur.ax,...
                           xData, yData(ch, :, zIdx));    
                    case 'plot'
                       this.requests.(gName).handles.cur.plot = ...
                           plot(this.requests.(gName).handles.cur.ax,...
                           xData, yData(ch, :, zIdx));  
                end
                this.setStringsToPlot(gName);
            else
                quickPlot(this, gName, xData, gather(yData(ch, :, zIdx)), [], [])
                drawnow()
            end
        end
        
        function plotFFT(this, gName, xData, yData, markerX, markerY)
            if ~isgraphics(this.requests.(gName).handles.cur.ax)
               return
            end
            
            zIdx = this.requests.zIdx;
            fidxs =  this.requests.fIdx-200:this.requests.fIdx+200;
            zVec = this.requests.zVec;
            
            xDataReduced = xData(fidxs);
            yDataReduced = yData(:, fidxs, zIdx);
            
            this.requests.(gName).lims.xlims = ...
                [this.requests.fBar(fidxs(1)), this.requests.fBar(fidxs(end))];
            this.setTitleVariables(gName, {[zVec(zIdx)*1e3, zIdx]});          % FFT
            
            if ~isgraphics(this.requests.(gName).handles.cur.plot(1))
            	hold(this.requests.(gName).handles.cur.ax, 'on');
                for i = 1:size(yData, 1)
                    switch this.requests.(gName).type
                        case 'stem'
                           this.requests.(gName).handles.cur.plot(i) = ...
                               stem(this.requests.(gName).handles.cur.ax,...
                               xDataReduced, yDataReduced(i, :, :));    
                        case 'plot'
                           this.requests.(gName).handles.cur.plot(i) = ...
                               plot(this.requests.(gName).handles.cur.ax,...
                               xDataReduced, yDataReduced(i, :, :));  
                    end
                end
                
                for i = 1:size(yData, 1)            
                    this.plotMarkers(gName, markerX, markerY(i, zIdx), i)
                end
                hold(this.requests.(gName).handles.cur.ax, 'off');
                
                this.setStringsToPlot(gName)
                this.setLimsToPlot(gName)
            else

                for i = 1:size(yData, 1)
                    this.quickPlot(gName, xDataReduced, yDataReduced(i, :, :), markerX, markerY(i, zIdx), i)
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
                this.requests.(gName).handles.cur.plotMarker(hIdx) = ...
                                plot(this.requests.(gName).handles.cur.ax,...
                                     xMarker, yMarker, 'g+');
            end
        end
        
        
    end
end
