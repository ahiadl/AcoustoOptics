classdef consGraphicsNew < Graphics
    %CONSGRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end
    
    methods (Static)
        function gNames = getGraphicsNames()
            gNames = {'allFrames'; 'curFrame'; 'curSet'; 'navFrame'; 'navSet'};
        end
        
        function graphReq = getGraphicsRequest()
            
            plotNames = consGraphics.getGraphicsNames();
            
            graphReq.quant = 1;
            graphReq.zIdx = 1;
            graphReq.ch   = 1;
            
            graphReq.curFrameIdx = 1;
            graphReq.curSetIdx   = 1;
            graphReq.navFrameIdx = 1;
            graphReq.navSetIdx   = 1;
            graphReq.numOfQuant  = 1;
            graphReq.timeFrames  = 1;
            
            graphReq.intExt  = 'int';
            
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

                graphReq.(plotNames{i}).dims.posDim   = [];
                graphReq.(plotNames{i}).dims.frameDim = [];
                graphReq.(plotNames{i}).dims.setDim   = [];
                graphReq.(plotNames{i}).dims.quantDim = [];
                graphReq.(plotNames{i}).dims.chDim    = [];
                graphReq.(plotNames{i}).dims.dataDim  = [];

                graphReq.(plotNames{i}).markers.plotMark = false;

                graphReq.(plotNames{i}).lims.xlims  = [];
                graphReq.(plotNames{i}).lims.ylims  = [];
                graphReq.(plotNames{i}).lims.colors = [];

                graphReq.(plotNames{i}).fonts.type       = [];
                graphReq.(plotNames{i}).fonts.titleSize  = 18;
                graphReq.(plotNames{i}).fonts.labelsSize = 18;
                graphReq.(plotNames{i}).fonts.axisSize   = 18;
            end
        end 
        
        function gReq = createGraphicsRunVars() %for Reduced Operation
            gNames = consGraphicsNew.getGraphicsNames();
            
            gReq.ch  = 1;
            gReq.zIdx = 1;
            gReq.quant = 1;
            
            gReq.navFrame = 0;
            gReq.navSet   = 0;
            
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
        function this = consGraphicsNew()
            this@Graphics(consGraphics.getGraphicsNames())
            
            this.requests     = consGraphicsNew.getGraphicsRequest();
            this.globalReqOld = consGraphicsNew.createGraphicsRunVars();
            this.globalReq    = consGraphicsNew.createGraphicsRunVars();
            
            this.setGraphicsStaticVars();
        end
        
        %Display Functions
        
        function dispTotalFrame(this, xData, yData, stdVec)
            if ~isgraphics(this.requests.allFrames.handles.cur.ax)
               return
            end
       
            zIdx = this.requests.zIdx;
            
            if ~isgraphics(this.requests.allFrames.handles.cur.plot)
                this.requests.allFrames.handles.cur.plot = ...
                           errorbar(this.requests.allFrames.handles.cur.ax,...
                           xData, yData(zIdx, :), stdVec(zIdx, :));  
                setStringsToPlot(this, 'allFrames')
                drawnow();
            else
                quickPlot(this, 'allFrames', xData, yData(zIdx, :), stdVec(zIdx, :))
            end
        end
        
        function dispCurrentFrame(this, xData, yData, stdVec)
            if ~isgraphics(this.requests.curFrame.handles.cur.ax)
               return
            end
            
            zIdx   = this.requests.zIdx;
            frame = this.requests.curFrameIdx;
            
            this.setTitleVariables('curFrame', {[this.requests.timeFrames(frame)]});
            
            if ~isgraphics(this.requests.curFrame.handles.cur.plot) ||...
                ~strcmp(this.requests.curFrame.type, this.requests.curFrame.handles.cur.plot.Type)
                switch this.requests.curFrame.type
                    case 'stem'
                       this.requests.curFrame.handles.cur.plot = ...
                           stem(this.requests.curFrame.handles.cur.ax,...
                           xData, squeeze(yData(zIdx, frame, :)));    
                    case 'errorbar'
                       this.requests.curFrame.handles.cur.plot = ...
                           errorbar(this.requests.curFrame.handles.cur.ax,...
                           xData, squeeze(yData(zIdx, frame, :)), squeeze(stdVec(zIdx, frame, :)));  
                end
                setStringsToPlot(this, 'curFrame')
                drawnow();
            else
                quickPlot(this, 'curFrame', xData, squeeze(yData(zIdx, frame, :)), squeeze(stdVec(zIdx, frame, :)));
            end
        end
        
        function dispNavFrame(this, xData, yData, stdVec)
            if ~isgraphics(this.requests.navFrame.handles.cur.ax)
               return
            end
            
            zIdx   = this.requests.zIdx;
            frame  = this.requests.navFrameIdx;

            this.setTitleVariables('navFrame', {[this.requests.timeFrames(frame)]});
            
            if ~isgraphics(this.requests.navFrame.handles.cur.plot) ||...
               ~strcmp(this.requests.navFrame.type, this.requests.navFrame.handles.cur.plot.Type)
                switch this.requests.navFrame.type
                    case 'stem'
                       this.requests.navFrame.handles.cur.plot = ...
                           stem(this.requests.navFrame.handles.cur.ax,...
                           xData, squeeze(yData(zIdx, frame, :)));    
                    case 'errorbar'
                       this.requests.navFrame.handles.cur.plot = ...
                           errorbar(this.requests.navFrame.handles.cur.ax,...
                           xData, squeeze(yData(zIdx, frame, :)), squeeze(stdVec(zIdx, frame, :)));  
                end
                setStringsToPlot(this, 'navFrame');
                drawnow();
            else
                quickPlot(this, 'navFrame', xData, yData(zIdx, frame, :), stdVec(zIdx, frame, :))
            end
        end
        
        function dispCurrentSet(this, xData, yData)
            if ~isgraphics(this.requests.curSet.handles.cur.ax)
               return
            end
            
            zIdx  = this.requests.zIdx;
            frame = this.requests.curFrameIdx;
            set   = this.requests.curSetIdx;   
            
            this.setTitleVariables('curSet', {[set]});
            
            
            if ~isgraphics(this.requests.curSet.handles.cur.plot)
                this.requests.curSet.handles.cur.plot = ...
                    stem(this.requests.curSet.handles.cur.ax,...
                    xData, squeeze(yData(zIdx, frame, set, 1:this.requests.numOfQuant(frame)))); 
                setStringsToPlot(this, 'curSet');
                drawnow();
            else
                quickPlot(this, 'curSet', xData,  squeeze(yData(zIdx, frame, set, 1:this.requests.numOfQuant(frame))),  []);
            end
        end
        
        function dispNavSet(this, xData, yData)
            if ~isgraphics(this.requests.navSet.handles.cur.ax)
               return
            end
            
            zIdx  = this.requests.zIdx;
            frame = this.requests.navFrameIdx;
            set   = this.requests.navSetIdx;
            
            this.setTitleVariables('navSet', {[this.requests.timeFrames(frame), set]});
            
            if ~isgraphics(this.requests.navSet.handles.cur.plot)
%                ~strcmp(this.requests.navSet.type, this.requests.navSet.handles.cur.plot.Type)
                this.requests.navSet.handles.cur.plot = ...
                    stem(this.requests.navSet.handles.cur.ax,...
                    xData, squeeze(yData(zIdx, frame, set, 1:this.requests.numOfQuant(frame))));
                setStringsToPlot(this, 'navSet');
                drawnow();
            else
                quickPlot(this, 'navSet', xData,  yData(zIdx, frame, set, 1:this.requests.numOfQuant(frame)), [])
            end
        end
        
        function quickPlot(this, gName, xData, yData, stdDev) 
            switch this.requests.(gName).type
                case 'stem'
                    set(this.requests.(gName).handles.cur.plot,...
                        'XData', xData,...
                        'YData', yData);
                case 'plot'    
                    set(this.requests.(gName).handles.cur.plot,...
                        'XData', xData,...
                        'YData', yData);
                case 'errorbar'
                    set(this.requests.(gName).handles.cur.plot,...
                        'XData', xData,...
                        'YData', yData,...
                        'YPositiveDelta', stdDev,...
                        'YNegativeDelta', stdDev);
            end 
            set(this.requests.(gName).handles.cur.title, 'String', this.requests.(gName).strings.title);
            drawnow();
         end
        
        % Set Functiond
        
        function setDims(this, gName, posDim, frameDim, setDim, quantDim, chDim, dataDim)
            this.requests.(gName).dims.posDim   = posDim;
            this.requests.(gName).dims.chDim    = chDim;
            this.requests.(gName).dims.frameDim = frameDim;
            this.requests.(gName).dims.setDim   = setDim;
            this.requests.(gName).dims.quantDim = quantDim;
            this.requests.(gName).dims.dataDim  = dataDim;
        end
        
        function setCurrentFrameAndSet(this, curFrame, curSet)
            this.requests.curFrameIdx = curFrame;
            this.requests.curSetIdx   = curSet;
        end
        
        function setNavFrame(this, frame)
            this.requests.navFrameIdx = frame;
        end
        
        function setNavSet(this, set)
            this.requests.navSetIdx = set;
        end
        
        function setTimeFramesAndQuants(this, timeFrames, numOfQuant)
            this.requests.timeFrames = timeFrames;
            this.requests.numOfQuant = numOfQuant;
        end
        
        function setGraphicsStaticVars(this)
            % All frames
            this.setType(this.graphNames{1}, 'errorbar');
            this.setStrings(this.graphNames{1}, {"All Frames Statistics"}, "Frame Duration[s]", "Amp", []);
            this.setDims(this.graphNames{1}, 1, 2, 3, 4, 5, 2)
            
            % Cur Frame
            this.setType(this.graphNames{2}, 'errorbar');
            this.setStrings(this.graphNames{2}, {"Current Frame [%.2f s] Sets"}, "Set[#]", "Amp", []);
            this.setDims(this.graphNames{2}, 1, 2, 3, 4, 5, 3)
           
            % Cur Set
            this.setType(this.graphNames{3}, 'stem');
            this.setStrings(this.graphNames{3}, {"Current Set [%d] Quants"}, "Quant", "Amp[V]", []);
            this.setDims(this.graphNames{3}, 1, 2, 3, 4, 5, 4)

            % navigated Frame
            this.setType(this.graphNames{4}, 'errorbar');
            this.setStrings(this.graphNames{4}, {"Frames Navigator: Frame %.2f[s]"}, "Frame Duration[s]", "Amp[V]", []);
            this.setDims(this.graphNames{4}, 1, 2, 3, 4, 5, 2)
            
            % Navigated Set
            this.setType(this.graphNames{5}, 'stem');
            this.setStrings(this.graphNames{5}, {"Set Navigator: Frame %.2f[s], Set: %d"}, "Set [#]", "Amp[V]", []);
            this.setDims(this.graphNames{5}, 1, 2, 3, 4, 5, 3)
        end
        
    end
end

