classdef scan2dGraphics < Graphics
    %SACN2DGRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods (Static)
  
       function gNames = getGraphicsNames()
            gNames = {'curMainAxis'; 'curMainAxisRep'; 'curMainPlain'; 'curMainPlainRep'; 'navPlain'; 'navPlainRep'};
       end 
        
       function graphReq = getGraphicsRequest()
            graphReq.zIdx    = 1;
            graphReq.ch      = 1;
            
            graphReq.mainAx = 'Y';
            
            graphReq.xAxis = [];
            graphReq.yAxis = [];
            graphReq.zAxis = [];
            
            graphReq.curXidx  = 0;
            graphReq.curYidx  = 0;
            graphReq.curQuant = 0;
            graphReq.curRep   = 0;
            
            graphReq.navXidx  = 0;
            graphReq.navYidx  = 0;
            graphReq.navZidx  = 0;
            graphReq.navRep   = 0;
            
            graphReq.navPlainCoor    = 'Z';
            graphReq.curNavPlainCoor = 'Z';
            
            graphReq.intExt  = 'int';
            plotNames = scan2dGraphics.getGraphicsNames();
            for i = 1:length(plotNames)  
                graphReq.(plotNames{i}).type      = []; 
                graphReq.(plotNames{i}).update    = true;
 
                graphReq.(plotNames{i}).handles.int = Graphics.createHandlesStruct();
                graphReq.(plotNames{i}).handles.ext = Graphics.createHandlesStruct();
                graphReq.(plotNames{i}).handles.cur = Graphics.createHandlesStruct();
                
                graphReq.(plotNames{i}).strings.title       = [];
                graphReq.(plotNames{i}).strings.titleModel  = [];
                graphReq.(plotNames{i}).strings.xlabel      = [];
                graphReq.(plotNames{i}).strings.xlabelModel = [];
                graphReq.(plotNames{i}).strings.ylabel      = [];
                graphReq.(plotNames{i}).strings.ylabelModel = [];
                graphReq.(plotNames{i}).strings.legend      = [];
                graphReq.(plotNames{i}).strings.legendModel = [];

                graphReq.(plotNames{i}).dims.xDim     = [];
                graphReq.(plotNames{i}).dims.yDim     = [];
                graphReq.(plotNames{i}).dims.zDim     = [];
                graphReq.(plotNames{i}).dims.chDim    = [];
                graphReq.(plotNames{i}).dims.quantDim = [];
                graphReq.(plotNames{i}).dims.repDim   = [];
                graphReq.(plotNames{i}).dims.mainDim  = [];

                graphReq.(plotNames{i}).markers.plotMark = false;

                graphReq.(plotNames{i}).lims.xlims  = [];
                graphReq.(plotNames{i}).lims.ylims  = [];
                graphReq.(plotNames{i}).lims.zlims  = [];
                graphReq.(plotNames{i}).lims.clims = [];

                graphReq.(plotNames{i}).fonts.type       = [];
                graphReq.(plotNames{i}).fonts.titleSize  = 18;
                graphReq.(plotNames{i}).fonts.labelsSize = 18;
                graphReq.(plotNames{i}).fonts.axisSize   = 18;
            end
       end  
        
       function gReq = createGraphicsRunVars() %for Reduced Operation
            gReq.ch  = 1;
            gReq.zIdx = 1;
     
            gReq.navXidx = 0;
            gReq.navYidx = 0;
            gReq.navZidx = 0;
            gReq.navRep  = 0;
            
            gReq.mainAx = 'Y';
            
            gNames = scan2dGraphics.getGraphicsNames();
            
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
        function this = scan2dGraphics()
            this@Graphics(scan2dGraphics.getGraphicsNames())
            this.requests     = scan2dGraphics.getGraphicsRequest();
            this.globalReqOld = scan2dGraphics.createGraphicsRunVars();
            this.globalReq    = scan2dGraphics.createGraphicsRunVars();
            this.setGraphicsStaticVars();
        end
 
        function setGraphicsStaticVars(this)
            % curMainAxis
            this.setType(this.graphNames{1}, 'stem');
            this.setStrings(this.graphNames{1}, {"Current Main Axis"; "(R,X,Y,Z) = (%d, %.2f, %.2f, %.2f)"}, "%s[mm]", "Amp[v]", []);
            this.setDims(this.graphNames{1}, 1, 2, 3, 4, 5)
            
            % curMainAxisRep
            this.setType(this.graphNames{2}, 'errorbar');
            this.setStrings(this.graphNames{2}, {"Current Main Axis (Averaged)"; "(R,X,Y,Z) = (%d, %.2f, %.2f, %.2f)"}, "%s[mm]",  "Amp[v]", []);
            this.setDims(this.graphNames{2}, 1, 2, 3, 4, 5)
            
            % curMainPlain
            this.setType(this.graphNames{3}, 'imagesc');
            this.setStrings(this.graphNames{3}, {"Current Main Plain"; "(R,X,Y,Z) = (%d, %.2f, %.2f, %.2f)"}, "%s [mm]", "%s[mm]", []);
            this.setDims(this.graphNames{3}, 1, 2, 3, 4, [])
            
            % curMainPlainRep
            this.setType(this.graphNames{4}, 'imagesc');
            this.setStrings(this.graphNames{4}, {"Current Main Plain (Averaged)"; "(R,X,Y,Z) = (%d, %.2f, %.2f, %.2f)"}, "%s [mm]", "%s[mm]", []);
            this.setDims(this.graphNames{4}, 1, 2, 3, 4, [])
           
            % navPlane
            this.setType(this.graphNames{5}, 'imagesc');
            this.setStrings(this.graphNames{5}, {"Navigator Plain"; "(R,%s) = (%d, %.2f)"}, "%s[mm]", "%s[mm]", []);
            this.setDims(this.graphNames{5}, 1, 2, 3, [], [])
            
            % navPlaneRep
            this.setType(this.graphNames{6}, 'imagesc');
            this.setStrings(this.graphNames{6}, {"Navigator Plain (Averaged)"; "(%s) = (%.2f)"}, "%s[mm]", "%s[mm]", []);
            this.setDims(this.graphNames{6}, 1, 2, 3, [], [])
        end
        
        function setDims(this, gName, xDim, yDim, zDim, quantDim, chDim)
            this.requests.(gName).dims.xDim     = xDim;
            this.requests.(gName).dims.yDim     = yDim;
            this.requests.(gName).dims.zDim     = zDim;
            this.requests.(gName).dims.quantDim = quantDim;
            this.requests.(gName).dims.chDim    = chDim;
        end
        
        % Display Function
        function dispCurMainAxis(this, yData, stdMat)
            % yData - (xLen, yLen, zLen, r)
            % stdMat - (xLen, yLen, zLen, r)
            if ~isgraphics(this.requests.curMainAxis.handles.cur.ax)
               return
            end
            
            zIdx = this.requests.zIdx;
            xIdx = this.requests.curXidx;
            yIdx = this.requests.curYidx;
            rIdx = this.requests.curRep;
            
            xCoor = this.requests.xAxis(xIdx);
            yCoor = this.requests.yAxis(yIdx);
            zCoor = this.requests.zAxis(zIdx)*1e3;
            
            this.setTitleVariables('curMainAxis', {[]; [rIdx, xCoor, yCoor, zCoor]})
            
            if strcmp(this.requests.mainAx, 'Y')
                yData = yData(xIdx, :, zIdx, rIdx);
                stdMat = stdMat(xIdx, :, zIdx, rIdx);
                xData = this.requests.yAxis;
                xName = 'Y'; yName = 'X';
                dx = abs(this.requests.yAxis(1) - this.requests.yAxis(2));
                xlim = [this.requests.yAxis(1)-dx, this.requests.yAxis(end)+dx];
            elseif strcmp(this.requests.mainAx, 'X')
                yData = yData(:, yIdx, zIdx, rIdx);
                stdMat = stdMat(:, yIdx, zIdx, rIdx);
                xData = this.requests.xAxis;
                xName = 'X'; yName = 'Y'; 
                dx = abs(this.requests.xAxis(1) - this.requests.xAxis(2));
                xlim = [this.requests.xAxis(end)-dx, this.requests.xAxis(1)+dx];
            end

            if ~isgraphics(this.requests.curMainAxis.handles.cur.plot) ||...
               this.requests.curMainAxis.update
                
                this.setAxesVar('curMainAxis', xName, yName);
                this.setLimits('curMainAxis', xlim, []);
                
                switch this.requests.curMainAxis.type
                    case 'stem'
                        this.requests.curMainAxis.handles.cur.plot = ...
                            stem(this.requests.curMainAxis.handles.cur.ax,...
                            xData, yData); 
                    case 'errorbar'
                        this.requests.curMainAxis.handles.cur.plot = ...
                            errorbar(this.requests.curMainAxis.handles.cur.ax,...
                            xData, yData, stdMat);
                end
                
                this.setLimsToPlot('curMainAxis')
                this.setStringsToPlot('curMainAxis');
                drawnow();
                this.requests.curMainAxis.update = false;
            else
                quickPlot(this, 'curMainAxis', xData, yData, [], stdMat)
            end
        end
        
        function dispCurMainAxisRep(this, yData, stdMat)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.requests.curMainAxisRep.handles.cur.ax)
               return
            end
            
            zIdx = this.requests.zIdx;
            xIdx = this.requests.curXidx;
            yIdx = this.requests.curYidx;
            rIdx = this.requests.curRep;
            
            xCoor = this.requests.xAxis(xIdx);
            yCoor = this.requests.yAxis(yIdx);
            zCoor = this.requests.zAxis(zIdx)*1e3;
            
            this.setTitleVariables('curMainAxisRep', {[]; [rIdx, xCoor, yCoor, zCoor]})
            
            if strcmp(this.requests.mainAx, 'Y')
                yData  = yData(xIdx, :, zIdx);
                stdMat = stdMat(xIdx, :, zIdx);
                xData  = this.requests.yAxis;
                xName = 'Y'; yName = 'X';
                dx = abs(this.requests.yAxis(1) - this.requests.yAxis(2));
                xlim = [this.requests.yAxis(1)-dx, this.requests.yAxis(end)+dx];
            elseif strcmp(this.requests.mainAx, 'X')
                yData  = yData(:, yIdx, zIdx);
                stdMat = stdMat(:, yIdx, zIdx);
                xData  = this.requests.xAxis;
                xName = 'X'; yName = 'Y'; 
                dx = abs(this.requests.xAxis(1) - this.requests.xAxis(2));
                xlim = [this.requests.xAxis(end)-dx, this.requests.xAxis(1)+dx];
            end

            if ~isgraphics(this.requests.curMainAxisRep.handles.cur.plot) ||...
               this.requests.curMainAxisRep.update
                
                this.setAxesVar('curMainAxisRep', xName, yName)
                this.setLimits('curMainAxisRep', xlim, [])
                switch this.requests.curMainAxisRep.type
                    case 'stem'
                        this.requests.curMainAxisRep.handles.cur.plot = ...
                            stem(this.requests.curMainAxisRep.handles.cur.ax,...
                            xData, yData); 
                    case 'errorbar'
                        this.requests.curMainAxisRep.handles.cur.plot = ...
                            errorbar(this.requests.curMainAxisRep.handles.cur.ax,...
                            xData, yData, stdMat);
                end
                this.setLimsToPlot('curMainAxisRep')
                this.setStringsToPlot('curMainAxisRep');
                drawnow();
                this.requests.curMainAxisRep.update = false;
            else
                quickPlot(this, 'curMainAxisRep', xData, yData, [], stdMat)
            end
        end
        
        function dispCurMainPlain(this, cData)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.requests.curMainPlain.handles.cur.ax)
               return
            end
            
            zIdx = this.requests.zIdx;
            xIdx = this.requests.curXidx;
            yIdx = this.requests.curYidx;
            rIdx = this.requests.curRep;
            
            xCoor = this.requests.xAxis(xIdx);
            yCoor = this.requests.yAxis(yIdx);
            zCoor = this.requests.zAxis(zIdx)*1e3;
            
            this.setTitleVariables('curMainPlain', {[]; [rIdx, xCoor, yCoor, zCoor]})
            
            clims = [min(min(min(min(cData)))), max(max(max(max(cData))))];
            if clims(1) == clims(2)
                clims(2) = clims(1)+ 0.1;
            end
            this.requests.curMainPlain.lims.clims = clims;

            if strcmp(this.requests.mainAx, 'Y')
                cData  = permute(cData(xIdx, :, :, rIdx), [3,2,1,4]);
                xData  = this.requests.yAxis;
                xName = 'Y'; yName = 'Z';
            elseif strcmp(this.requests.mainAx, 'X')
                cData  = permute(cData(:, yIdx, :, rIdx), [3,1,2,4]);                
                xData  = this.requests.xAxis;
                xName = 'X'; yName = 'Z'; 
            end
            yData  = this.requests.zAxis;
            
            if ~isgraphics(this.requests.curMainPlain.handles.cur.plot) ||...
                this.requests.curMainPlain.update
                
                cla(this.requests.curMainPlain.handles.cur.ax)
                this.setAxesVar('curMainPlain', xName, yName);
                
                this.requests.curMainPlain.handles.cur.plot = ...
                    imagesc(this.requests.curMainPlain.handles.cur.ax,...
                    'XData', xData, 'Ydata', yData, 'CData', cData,...
                    this.requests.curMainPlain.lims.clims);
                axis(this.requests.curMainPlain.handles.cur.ax, 'tight')
                colorbar(this.requests.curMainPlain.handles.cur.ax);
                this.setLimsToPlot('curMainPlain')
                this.setStringsToPlot('curMainPlain');
                drawnow();
                this.requests.curMainPlain.update = false;
            else
                quickPlot(this, 'curMainPlain', xData, yData, cData)
            end
        end
        
        function dispCurMainPlainRep(this, cData)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.requests.curMainPlainRep.handles.cur.ax)
               return
            end
            
            zIdx = this.requests.zIdx;
            xIdx = this.requests.curXidx;
            yIdx = this.requests.curYidx;
            rIdx = this.requests.curRep;
            
            xCoor = this.requests.xAxis(xIdx);
            yCoor = this.requests.yAxis(yIdx);
            zCoor = this.requests.zAxis(zIdx)*1e3;
            
            this.setTitleVariables('curMainPlain', {[]; [rIdx, xCoor, yCoor, zCoor]})
            
            clims = [min(min(min(cData))), max(max(max(cData)))];
            if (clims(1) == clims(2))
                clims(2) = clims(1)+0.1;
            end
            this.requests.curMainPlainRep.lims.clims = clims;
            
            if strcmp(this.requests.mainAx, 'Y')
                cData  = permute(cData(xIdx, :, :), [3,2,1]);
                xData  = this.requests.yAxis;
                xName = 'Y'; yName = 'Z';
            elseif strcmp(this.requests.mainAx, 'X')
                cData  = permute(cData(:, yIdx, :), [3,1,2]);                
                xData  = this.requests.xAxis;
                xName = 'X'; yName = 'Z'; 
            end
            yData  = this.requests.zAxis;
            
            if ~isgraphics(this.requests.curMainPlainRep.handles.cur.plot) ||...
                this.requests.curMainPlain.update
                
                cla(this.requests.curMainPlainRep.handles.cur.ax)
                this.setAxesVar('curMainPlainRep', xName, yName);
                
                this.requests.curMainPlainRep.handles.cur.plot = ...
                    imagesc(this.requests.curMainPlainRep.handles.cur.ax,...
                    'XData', xData, 'Ydata', yData, 'CData', cData,...
                    this.requests.curMainPlainRep.lims.clims);
                axis(this.requests.curMainPlainRep.handles.cur.ax, 'tight')
                colorbar(this.requests.curMainPlainRep.handles.cur.ax);
                this.setLimsToPlot('curMainPlainRep')
                this.setStringsToPlot('curMainPlainRep');
                drawnow();
                this.requests.curMainPlainRep.update = false;
            else
                quickPlot(this, 'curMainPlainRep', xData, yData, cData)
            end
        end
        
        function dispNavPlain(this, cData)
            % yData  - (xLen, yLen, zLen, rep)
            if ~isgraphics(this.requests.navPlain.handles.cur.ax)
               return
            end
            
            clims = [min(min(min(min(cData)))), max(max(max(max(cData))))];
            if (clims(1) == clims(2)); clims(2) = clims(1)+0.1; end
            this.requests.navPlain.lims.clims = clims;
            
            rIdx = this.requests.navRep;
%             this.requests.navPlain.update = strcmp(this.requests.curNavPlainCoor, this.requests.navPlainCoor);
            switch this.requests.navPlainCoor
                case 'X'
                    xIdx = this.requests.navXidx;
                    xCoor = this.requests.xAxis(xIdx);
                    
                    str{1} = this.requests.navPlain.strings.titleModel{1};
                    str{2} = sprintf(this.requests.navPlain.strings.titleModel{2}, 'X', rIdx, xCoor);
                    this.requests.navPlain.strings.title = str;
                    
                    xData  = this.requests.yAxis;
                    yData  = this.requests.zAxis;
                    cData  = permute(cData(xIdx, :, :, rIdx), [3,2,1,4]);
                    
                    xName = 'Y';
                    yName = 'Z';
                case 'Y'
                    yIdx = this.requests.navYidx;
                    yCoor = this.requests.yAxis(yIdx);
                    
                    str{1} = this.requests.navPlain.strings.titleModel{1};
                    str{2} = sprintf(this.requests.navPlain.strings.titleModel{2}, 'Y', rIdx, yCoor);
                    this.requests.navPlain.strings.title = str;
                    
                    xData  = this.requests.xAxis;
                    yData  = this.requests.zAxis;
                    cData  = permute(cData(:, yIdx, :, rIdx), [3,1,2,4]);
                    xName = 'X';
                    yName = 'Z';
                case 'Z'
                    zIdx  = this.requests.navZidx;
                    zCoor = this.requests.zAxis(zIdx)*1e3;
                    
                    str{1} = this.requests.navPlain.strings.titleModel{1};
                    str{2} = sprintf(this.requests.navPlain.strings.titleModel{2}, 'Z', rIdx, zCoor);
                    this.requests.navPlain.strings.title = str;
                    
                    xData  = this.requests.xAxis;
                    yData  = this.requests.yAxis;
                    cData  = permute(cData(:, :, zIdx, rIdx), [2,1,3,4]);
                    
                    xName = 'X';
                    yName = 'Y';
            end
            
            if ~isgraphics(this.requests.navPlain.handles.cur.plot) ||...
                this.requests.navPlain.update
                
                cla(this.requests.navPlain.handles.cur.ax)
                this.setAxesVar('navPlain', xName, yName);
                
                this.requests.navPlain.handles.cur.plot = ...
                    imagesc(this.requests.navPlain.handles.cur.ax,...
                    'XData', xData, 'Ydata', yData, 'CData', cData,...
                    this.requests.navPlain.lims.clims);
                axis(this.requests.navPlain.handles.cur.ax, 'tight')
                colorbar(this.requests.navPlain.handles.cur.ax);
                this.setLimsToPlot('navPlain')
                this.setStringsToPlot('navPlain');
                drawnow();
                this.requests.navPlain.update = false;
            else
                quickPlot(this, 'navPlain', xData, yData, cData)
            end    
        end
        
        function dispNavPlainRep(this, cData)
            % yData  - (xLen, yLen, zLen, rep)
            if ~isgraphics(this.requests.navPlainRep.handles.cur.ax)
               return
            end
            
            clims = [min(min(min(cData))), max(max(max(cData)))];
            if (clims(1) == clims(2)); clims(2) = clims(1)+0.1; end
            this.requests.navPlainRep.lims.clims = clims;
            
%             this.requests.navPlain.update = strcmp(this.requests.curNavPlainCoor, this.requests.navPlainCoor);
            switch this.requests.navPlainCoor
                case 'X'
                    xIdx = this.requests.navXidx;
                    xCoor = this.requests.xAxis(xIdx);
                    
                    str{1} = this.requests.navPlainRep.strings.titleModel{1};
                    str{2} = sprintf(this.requests.navPlainRep.strings.titleModel{2}, 'X', xCoor);
                    this.requests.navPlainRep.strings.title = str;
                    
                    xData  = this.requests.yAxis;
                    yData  = this.requests.zAxis;
                    cData  = permute(cData(xIdx, :, :), [3,2,1]);
                    
                    xName = 'Y';
                    yName = 'Z';                    
                case 'Y'
                    yIdx = this.requests.navYidx;
                    yCoor = this.requests.yAxis(yIdx);
                    
                    str{1} = this.requests.navPlainRep.strings.titleModel{1};
                    str{2} = sprintf(this.requests.navPlainRep.strings.titleModel{2}, 'Y', yCoor);
                    this.requests.navPlainRep.strings.title = str;
                    
                    xData  = this.requests.xAxis;
                    yData  = this.requests.zAxis;
                    cData  = permute(cData(:, yIdx, :), [3,1,2]);
                    
                    xName = 'X';
                    yName = 'Z';
                case 'Z'
                    zIdx  = this.requests.navZidx;
                    zCoor = this.requests.zAxis(zIdx)*1e3;
                    
                    str{1} = this.requests.navPlainRep.strings.titleModel{1};
                    str{2} = sprintf(this.requests.navPlainRep.strings.titleModel{2}, 'Z', zCoor);
                    this.requests.navPlainRep.strings.title = str;
                    
                    xData  = this.requests.xAxis;
                    yData  = this.requests.yAxis;
                    cData  = permute(cData(:, :, zIdx), [2,1,3]);
                    
                    xName = 'X';
                    yName = 'Y';
            end

            if ~isgraphics(this.requests.navPlainRep.handles.cur.plot) ||...
                this.requests.navPlainRep.update
                
                cla(this.requests.navPlainRep.handles.cur.ax)
                this.setAxesVar('navPlainRep', xName, yName);
                
                this.requests.navPlainRep.handles.cur.plot = ...
                    imagesc(this.requests.navPlainRep.handles.cur.ax,...
                    'XData', xData, 'Ydata', yData, 'CData', cData,...
                    this.requests.navPlainRep.lims.clims);
                axis(this.requests.navPlainRep.handles.cur.ax, 'tight')
                colorbar(this.requests.navPlainRep.handles.cur.ax);
                this.setLimsToPlot('navPlainRep')
                this.setStringsToPlot('navPlainRep');
                drawnow();
                this.requests.navPlainRep.update = false;
            else
                quickPlot(this, 'navPlainRep', xData, yData, cData)
            end    
        end
        
        %Controlling Functions
        function updateCurScan(this, curScan)
            this.requests.curRep   = curScan(1);
            this.requests.curXidx  = curScan(2);
            this.requests.curYidx  = curScan(3);
%             this.requests.curQuant = curScan(4);
        end
        
        function setAxesVec(this, xAxis, yAxis, zAxis)
            this.requests.xAxis = xAxis;
            this.requests.yAxis = yAxis;
            this.requests.zAxis = zAxis;
        end
        
        function setMainAx(this, ax)
           this.requests.mainAx = ax; 
        end
        
        function quickPlot(this, gName, xData, yData, cData, stdDev) 
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
                case 'imagesc'
                    set(this.requests.(gName).handles.cur.plot,...
                        'XData', xData,...
                        'YData', yData,...
                        'CData', cData);
                    this.setLimsToPlot(gName)
            end 
            set(this.requests.(gName).handles.cur.title, 'String', this.requests.(gName).strings.title);
            drawnow();
        end
        
        function setAxesVar(this, gName, xName, yName)
           this.requests.(gName).strings.xlabel = ...
               sprintf(this.requests.(gName).strings.xlabelModel, xName);
           this.requests.(gName).strings.ylabel = ...
               sprintf(this.requests.(gName).strings.ylabelModel, yName);
        end
        
        function setUpdate(this, name, update)
           this.requests.(name).update = update; 
        end
        
        function setNavParams(this, plain, idx, rep)
            if ~strcmp(this.requests.navPlainCoor, plain)
                this.requests.navPlainRep.update = true;
                this.requests.navPlain.update    = true;
            end
            this.requests.navPlainCoor = plain;
            this.requests.navRep = rep;
            
            switch plain
                case 'X'
                    this.requests.navXidx = idx;
                case 'Y'
                    this.requests.navYidx = idx;
                case 'Z'
                    this.requests.navZidx = idx;
            end
        end
    end
    
end

