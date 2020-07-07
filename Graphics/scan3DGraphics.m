classdef scan3DGraphics < Graphics
    %SACN2DGRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data
    end
    
    methods (Static)
  
       function figsNames = getGraphicsNames()
            figsNames = {'nav'; 'navAvg'};
       end 
        
       function figs  = createGraphicsVars()    
            % General Scan Parameters
            figs.firstAxType     = 'normal';
            figs.firstAxis       = [];
            figs.firstAxisCntr   = [];
            figs.firstAxisZero   = [];
            figs.firstAxisNorm   = [];
            
            figs.secondAxType     = 'normal';
            figs.secondAxis       = [];
            figs.secondAxisCntr   = [];
            figs.secondAxisZero   = [];
            figs.secondAxisNorm   = [];
            
            figs.depthAxType     = 'normal';
            figs.depthAxis       = [];
            figs.depthAxisCntr   = [];
            figs.depthAxisZero   = [];
            figs.depthAxisNorm   = [];
            
            figs.repeats = [];
            figs.normColorsToPlane = true;
            
            figs.firstAxLabel   = 'Y';
            figs.secondAxLabel  = 'X';
            figs.depthAxLabel   = 'Z';
            figs.mainPlaneLabel = 'YZ';
            
            figs.dxFirstAx         = [];
            
            % Current Position
            figs.curPosFirst  = 0;
            figs.curPosSecond = 0;
            figs.curIdxFirst  = 0;
            figs.curIdxSecond = 0;
            figs.curRep       = 0;            
            
            % Navigator Position
            figs.navAx    = 'Z';  % navigating on this axes
            figs.navPlane = 'XY'; % displaying these planes

            figs.navPos   = 0;
            figs.navIdx   = 0;
            figs.navRep   = 0;
            
            figs.navXaxis = [];
            figs.navYaxis = [];
            
            figs.navXname = 'X';
            figs.navYname = 'Y';
            
            % Color Limits
            figs.scan3DClims     = [0, 1];
            figs.scanAvg3DClims  = [0, 1];
            figs.scan2DClims     = [0, 1];
            figs.scanAvg2DClims  = [0, 1];
            
            %Specific Figures parameters
            figs.intExt  = 'int';
            figsNames    = scan3DGraphics.getGraphicsNames();
            
            figs.fonts.type       = [];
            figs.fonts.titleSize  = 14;
            figs.fonts.labelsSize = 14;
            figs.fonts.axisSize   = 14;
            
            for i = 1:length(figsNames)  
                figs.(figsNames{i}).type      = []; 
                figs.(figsNames{i}).update    = true;
 
                figs.(figsNames{i}).handles.int = Graphics.createHandlesStruct();
                figs.(figsNames{i}).handles.ext = Graphics.createHandlesStruct();
                figs.(figsNames{i}).handles.cur = Graphics.createHandlesStruct();
                
                figs.(figsNames{i}).strings.title       = [];
                figs.(figsNames{i}).strings.titleModel  = [];
                figs.(figsNames{i}).strings.xlabel      = [];
                figs.(figsNames{i}).strings.xlabelModel = [];
                figs.(figsNames{i}).strings.ylabel      = [];
                figs.(figsNames{i}).strings.ylabelModel = [];
                figs.(figsNames{i}).strings.legend      = [];
                figs.(figsNames{i}).strings.legendModel = [];

                figs.(figsNames{i}).markers.plotMark = false;

                figs.(figsNames{i}).lims.xlims  = [];
                figs.(figsNames{i}).lims.ylims  = [];
                figs.(figsNames{i}).lims.zlims  = [];
                figs.(figsNames{i}).lims.clims  = [];
            end
       end  
       
       function figs  = createUserVars()
           figs.normColorsToPlane = false;
           figs.reopenFigures     = false;
           
           figs.firstAxType  = 'normal';
           figs.secondAxType = 'normal';
           figs.depthAxType  = 'normal';
           
           figsNames = scan3DGraphics.getGraphicsNames();
           figs.intExt       = 'int';
           for i=1:length(figsNames)
               figs.validStruct.(figsNames{i}) = false;
           end

           for i=1:length(figsNames)
               figs.extH.(figsNames{i}) =  Graphics.createHandlesStruct();
           end
           
           figs.fonts.type       = [];
           figs.fonts.titleSize  = 14;
           figs.fonts.labelsSize = 14;
           figs.fonts.axisSize   = 14;
       end
       
       function uVars = createOwnerVars() %for Reduced Operation
            uVars = scan3DGraphics.createUserVars();
           
            uVars.repeats    = 1;
            
            uVars.firstAxisNorm   = [];
            uVars.secondAxisNorm  = [];
            uVars.depthAxisNorm   = [];

            uVars.firstAxLabel   = 'Y';
            uVars.secondAxLabel  = 'X';
            uVars.depthAxLabel   = 'Z';
            uVars.mainPlaneLabel = 'YZ';
       end
       
       function navVars = createNavigatorVars()
            navVars.navAx   = 'Z';
            navVars.navIdx  = 0;
            navVars.navRep  = 0;
       end
    end
    
    methods
        function this = scan3DGraphics()
            this@Graphics()
            this.figsNames = scan3DGraphics.getGraphicsNames();
            this.figs      = scan3DGraphics.createGraphicsVars();
            this.uVars     = scan3DGraphics.createOwnerVars();
            this.setGraphicsStaticVars();
            this.numOfFigs = length(this.figsNames);
        end
 
        % Set vars Functions
        function setGraphicsStaticVars(this)
            % navPlane
            this.setType('nav', 'imagesc');
            this.setStrings('nav', "Navigator Plane %s: (R, %s) = (%d, %.2f)", "%s[mm]", "%s[mm]", []);
            
            % navPlaneRep
            this.setType('navAvg', 'imagesc');
            this.setStrings('navAvg', "Navigator Plane %s (Averaged): (R, %s) = (%d, %.2f)", "%s[mm]", "%s[mm]", []);

        end
                 
        function setGraphicsScanVars(this)
            this.figs.repeats    = this.uVars.repeats;
            
            this.figs.firstAxLabel   = this.uVars.firstAxLabel;
            this.figs.secondAxLAbel  = this.uVars.secondAxLabel;
            this.figs.depthAxLabel   = this.uVars.depthAxLabel;
            this.figs.mainPlaneLabel = this.uVars.mainPlane;
            
            this.figs.firstAxType   = this.uVars.firstAxType;
            this.figs.firstAxisNorm = this.uVars.firstAxisNorm;
            this.figs.firstAxisCntr = this.figs.firstAxisNorm - mean(this.figs.firstAxisNorm);
            this.figs.firstAxisZero = abs(this.figs.firstAxisNorm - this.figs.firstAxisNorm(1));
            this.setAxisType(this.figs.firstAxLabel, this.figs.firstAxType);
            
            this.figs.secondAxType   = this.uVars.secondAxType;
            this.figs.secondAxisNorm = this.uVars.secondAxisNorm;
            this.figs.secondAxisCntr = this.figs.secondAxisNorm - mean(this.figs.secondAxisNorm);
            this.figs.secondAxisZero = abs(this.figs.secondAxisNorm - this.figs.secondAxisNorm(1));
            this.setAxisType(this.figs.secondAxLabel, this.figs.secondAxType);
            
            this.figs.depthAxType     = this.uVars.depthAxType;
            this.figs.depthAxisNorm   = this.uVars.depthAxisNorm*1e3;
            this.figs.depthAxisCntr   = this.figs.depthAxisNorm - mean(this.figs.depthAxisNorm);
            this.figs.depthAxisZero   = abs(this.figs.depthAxisNorm - this.figs.depthAxisNorm(1));
            this.setAxisType(this.figs.depthAxLabel, this.figs.depthAxType);
            
            this.figs.normColorsToPlane = this.uVars.normColorsToPlane;
            
            % reset data arrays:
            this.initDataArray();

            this.figs.dxFirstAx         = abs(this.figs.firstAxis(1) - this.figs.firstAxis(2));
            
            this.figs.validStruct = this.uVars.validStruct;

            %extH and intExt flag should not be transferred at this point.
            %it is updated by the updateGraphicsConstruction function
            this.resetNavigator();
        end
        
        function updateS3DCurPosAndIdx(this, curIdx, curPos)
            this.figs.curIdxSecond = curIdx;
            this.figs.curPosSecond = curPos;
        end
        
        function updateS2DCurPosAndIdx(this, pns)
            this.figs.curRep       = pns.curPosIdx(1);
            this.figs.curIdxFirst  = pns.curPosIdx(2);
            this.figs.curPosFirst  = pns.curPos(1);
        end

        function setNavVars(this, navVars)
            this.figs.navIdx =  navVars.navIdx;
            this.figs.navRep =  navVars.navRep;
            if ~strcmp( this.figs.navAx, navVars.navAx) || this.figs.navUpdate
                this.figs.nav.update = true;
                this.figs.navAvg.update = true;
                this.figs.navAx = navVars.navAx;
                if strcmp(this.figs.navAx, this.figs.firstAxLabel)
                    this.figs.navPlane = sprintf("%s%s", this.figs.secondAxLabel, this.figs.depthAxLabel);
                    this.figs.navAxis  = this.figs.firstAxis;
                    this.figs.navXaxis = this.figs.secondAxis;
                    this.figs.navYaxis = this.figs.depthAxis;
                    this.figs.navXname = this.figs.secondAxLabel;
                    this.figs.navYname = this.figs.depthAxLabel;
                elseif strcmp(this.figs.navAx, this.figs.secondAxLabel)
                    this.figs.navPlane = sprintf("%s%s", this.figs.firstAxLabel, this.figs.depthAxLabel);
                    this.figs.navAxis  = this.figs.secondAxis;
                    this.figs.navXaxis = this.figs.firstAxis;
                    this.figs.navYaxis = this.figs.depthAxis;
                    this.figs.navXname = this.figs.firstAxLabel;
                    this.figs.navYname = this.figs.depthAxLabel;
                elseif strcmp(this.figs.navAx, this.figs.depthAxLabel)
                    this.figs.navPlane = sprintf("%s%s", this.figs.firstAxLabel, this.figs.secondAxLabel);
                    this.figs.navAxis  = this.figs.depthAxis;
                    this.figs.navXaxis = this.figs.firstAxis;
                    this.figs.navYaxis = this.figs.secondAxis;
                    this.figs.navXname = this.figs.firstAxLabel;
                    this.figs.navYname = this.figs.secondAxLabel;
                else
                    fprintf("S2D Graphics: ERROR: no compatible navigator Axis");
                end
            end
            this.figs.navPos   = this.figs.navAxis(this.figs.navIdx);
        end
        
        function resetNavigator(this)
            navVars.navAx = this.figs.depthAxLabel;
            navVars.navIdx = 1;
            navVars.navRep = 1;
            this.figs.navUpdate = true;
            this.setNavVars(navVars);
            this.dispNavPlane();
            this.dispNavPlaneAvg();
        end
        
        function setLoadedDataClims(this, data)
            this.figs.scan3DClims(1) = min(data.phi(:));
            this.figs.scan3DClims(2) = max(data.phi(:));
            this.figs.scanAvg3DClims(1) = min(data.phiAvg(:));
            this.figs.scanAvg3DClims(2) = max(data.phiAvg(:));             
        end
        
        function setColorsScale(this, normColorsToPlane)
           this.figs.normColorsToPlane = normColorsToPlane; 
        end
        
        function setAxisType(this, ax, type)
            if strcmp(ax, this.figs.firstAxLabel)
                this.figs.firstAxType = type;
                switch type
                    case 'Normal'
                        this.figs.firstAxis = this.figs.firstAxisNorm;
                    case 'Center'
                        this.figs.firstAxis = this.figs.firstAxisCntr;
                    case 'Zero'
                        this.figs.firstAxis = this.figs.firstAxisZero;
                    case 'Index'
                        this.figs.firstAxis = 1:length(this.figs.firstAxisNorm);
                end
                this.figs.dxFirstAx         = abs(this.figs.firstAxis(1) - this.figs.firstAxis(2));
                this.figs.firstAxLimsToPlot = [min(this.figs.firstAxis) - this.figs.dxFirstAx,...
                               max(this.figs.firstAxis) + this.figs.dxFirstAx];
            elseif strcmp(ax, this.figs.secondAxLabel)
                this.figs.depthAxType = type;
                switch type
                    case 'Normal'
                        this.figs.secondAxis = this.figs.secondAxisNorm;
                    case 'Center'
                        this.figs.secondAxis = this.figs.secondAxisCntr;
                    case 'Zero'
                        this.figs.secondAxis = this.figs.secondAxisZero;
                    case 'Index' 
                        this.figs.secondAxis = 1:length(this.figs.secondAxisNorm);
                end
            elseif strcmp(ax, this.figs.depthAxLabel)
                this.figs.depthAxType = type;
                switch type
                    case 'Normal'
                        this.figs.depthAxis = this.figs.depthAxisNorm;
                    case 'Center'
                        this.figs.depthAxis = this.figs.depthAxisCntr;
                    case 'Zero'
                        this.figs.depthAxis = this.figs.depthAxisZero;
                    case 'Index' 
                        this.figs.firstAxis = 1:length(this.figs.depthAxisNorm);
                end
            end  
        end
        
        %Set Data functions
        function initDataArray(this)
            rstPhi    = zeros(length(this.figs.firstAxis), length(this.figs.secondAxis), length(this.figs.depthAxis), this.uVars.repeats);
            rstPhiAvg = zeros(length(this.figs.firstAxis), length(this.figs.secondAxis), length(this.figs.depthAxis));
            
            this.data.phi       = rstPhi;
            this.data.phiAvg    = rstPhiAvg;
            this.data.phiAVgStd = rstPhiAvg;
        end
        
        function set1DData(this, data)
            this.data.phi(this.figs.curIdxFirst, this.figs.curIdxSecond, :, this.figs.curRep) = data.phi;

            if (this.figs.curIdxFirst == 1)
               this.figs.scan2DClims  = [inf, 0];
            end
            
            if (this.figs.curIdxFirst == 1) && (this.figs.curIdxSecond == 1) && (this.figs.curRep == 1)
                this.figs.scan3DClims  = [inf, 0];
            end    

            curMin = min(this.data.phi(this.figs.curIdxFirst, this.figs.curIdxSecond, :, this.figs.curRep));
            curMax = max(this.data.phi(this.figs.curIdxFirst, this.figs.curIdxSecond, :, this.figs.curRep));
            this.figs.scan2DClims(1) = min(this.figs.scan2DClims(1), curMin);
            this.figs.scan2DClims(2) = max(this.figs.scan2DClims(2), curMax);
            this.figs.scan3DClims(1) = min(this.figs.scan3DClims(1), curMin);
            this.figs.scan3DClims(2) = max(this.figs.scan3DClims(2), curMax);
        end

        function setAvg1DData(this, data)
            this.data.phiAvg(this.figs.curIdxFirst, this.figs.curIdxSecond, :)     = data.phiAvg;
            this.data.phiAvgStd(this.figs.curIdxFirst, this.figs.curIdxSecond, :) = data.phiAvgStd;
            
            if (this.figs.curRep == 1) && (this.figs.curIdxFirst == 1)
               this.figs.scanAvg2DClims  = [inf, 0];
            end
            
            if (this.figs.curIdxFirst == 1) && (this.figs.curIdxSecond == 1) && (this.figs.curRep == 1)
                this.figs.scanAvg3DClims  = [inf, 0];
            end 
            
            curMin = min(this.data.phiAvg(this.figs.curIdxFirst, this.figs.curIdxSecond, :));
            curMax = max(this.data.phiAvg(this.figs.curIdxFirst, this.figs.curIdxSecond, :));
            this.figs.scanAvg2DClims(1) = min(this.figs.scanAvg2DClims(1), curMin);
            this.figs.scanAvg2DClims(2) = max(this.figs.scanAvg2DClims(2), curMax);
            this.figs.scanAvg3DClims(1) = min(this.figs.scanAvg3DClims(1), curMin);
            this.figs.scanAvg3DClims(2) = max(this.figs.scanAvg3DClims(2), curMax);
        end

        function setLoadedData(this, data)
            this.data.phi = data.phi;
            this.data.phiAvg = data.phiAvg;
            this.data.phiAvgStd = data.phiAvgStd;
            
            this.setLoadedDataClims(data);
        end
        
        % Extract data functions
        function data = extractNavPlaneFromData(this, isAvg)
            if ~isAvg
%                 data.clims = [min(min(min(min(this.data.phi)))), max(max(max(max(this.data.phi))))];
                if strcmp(this.figs.navAx, this.figs.firstAxLabel)
                    data.cData = permute(this.data.phi(this.figs.navIdx, :, :, this.figs.navRep), [3,2,1,4]);
                elseif strcmp(this.figs.navAx, this.figs.secondAxLabel)
                    data.cData  = permute(this.data.phi(:, this.figs.navIdx, :, this.figs.navRep), [3,1,2,4]);
                elseif strcmp(this.figs.navAx, this.figs.depthAxLabel)
                    data.cData  = permute(this.data.phi(:, :, this.figs.navIdx, this.figs.navRep), [2,1,3,4]); 
                else
                    fprintf("S2D ERROR: no compatible navigator axis");
                end
            else
%                 data.clims = [min(min(min(min(this.data.phiAvg)))), max(max(max(max(this.data.phiAvg))))];
                if strcmp(this.figs.navAx, this.figs.firstAxLabel)
                    data.cData  = permute(this.data.phiAvg(this.figs.navIdx, :, :), [3,2,1]);
                elseif strcmp(this.figs.navAx, this.figs.secondAxLabel)
                    data.cData  = permute(this.data.phiAvg(:, this.figs.navIdx, :), [3,1,2]);
                elseif strcmp(this.figs.navAx, this.figs.depthAxLabel)
                    data.cData  = permute(this.data.phiAvg(:, :, this.figs.navIdx), [2,1,3]); 
                else
                    fprintf("S2D ERROR: no compatible navigator axis");
                end
            end
            
            if this.figs.normColorsToPlane  
                uniData = unique(data.cData(:));
                minVars = mink(uniData, 2);
                if minVars(1) ~= 0
                    minVal = minVars(1);
                else
                    minVal = minVars(2);
                end
                data.clims(1) = minVal;
                data.clims(2) = max(uniData);
            else
                if ~isAvg
                    data.clims = this.figs.scan3DClims;
                else
                    data.clims = this.figs.scanAvg3DClims;
                end
            end
                        
            if (data.clims(1) == data.clims(2)) 
                data.clims(2) = data.clims(1)+0.1; 
            end
        end
        
        % Display Functions
        function updatePlots(this)
            this.dispNavPlane();
            this.dispNavPlaneAvg();
        end
        
        function dispNavPlane(this)
            % yData  - (xLen, yLen, zLen, rep)
            if ~isgraphics(this.figs.nav.handles.cur.ax)
               return
            end
            
            this.setTitleVariables('nav',...
                                    {{ this.figs.navPlane, this.figs.navAx, ...
                                       this.figs.navRep,   this.figs.navPos}});

            plotData = this.extractNavPlaneFromData(false);
            this.figs.nav.lims.clims = plotData.clims;
            
            if ~isgraphics(this.figs.nav.handles.cur.plot) ||...
                this.figs.nav.update
                
                cla(this.figs.nav.handles.cur.ax)
                this.setAxesVar('nav', this.figs.navXname, this.figs.navYname);
                
                this.figs.nav.handles.cur.plot = ...
                    imagesc(this.figs.nav.handles.cur.ax, 'XData', this.figs.navXaxis,...
                    'Ydata', this.figs.navYaxis, 'CData', plotData.cData,...
                    plotData.clims);
                
                axis(this.figs.nav.handles.cur.ax, 'tight')
                colorbar(this.figs.nav.handles.cur.ax);
                this.setLimsToPlot('nav')
                this.setStringsToPlot('nav');
                drawnow();
                this.figs.nav.update = false;
            else
                quickPlot(this, 'nav', this.figs.navXaxis, this.figs.navYaxis, plotData.cData)
            end
        end
        
        function dispNavPlaneAvg(this)
            % yData  - (xLen, yLen, zLen, rep)
            if ~isgraphics(this.figs.navAvg.handles.cur.ax)
               return
            end
            this.setTitleVariables('navAvg',...
                                   {{ this.figs.navPlane, this.figs.navAx,...
                                      this.figs.curRep,   this.figs.navPos}});

            plotData = this.extractNavPlaneFromData(true);
            this.figs.nav.lims.clims = plotData.clims;
                                        
            if ~isgraphics(this.figs.navAvg.handles.cur.plot) ||...
                this.figs.navAvg.update
                
                cla(this.figs.navAvg.handles.cur.ax)
                this.setAxesVar('navAvg', this.figs.navXname, this.figs.navYname);
                
                this.figs.navAvg.handles.cur.plot = ...
                    imagesc(this.figs.navAvg.handles.cur.ax, 'XData', this.figs.navXaxis,...
                    'Ydata', this.figs.navYaxis, 'CData', plotData.cData, plotData.clims);
                
                axis(this.figs.navAvg.handles.cur.ax, 'tight')
                colorbar(this.figs.navAvg.handles.cur.ax);
                this.setLimsToPlot('navAvg')
                this.setStringsToPlot('navAvg');
                drawnow();
                this.figs.navAvg.update = false;
            else
                quickPlot(this, 'navAvg', this.figs.navXaxis, this.figs.navYaxis, plotData.cData)
            end    
        end
        
        function quickPlot(this, gName, xData, yData, cData, stdDev) 
            switch this.figs.(gName).type
                case 'stem'
                    set(this.figs.(gName).handles.cur.plot,...
                        'XData', xData,...
                        'YData', yData);
                case 'plot'    
                    set(this.figs.(gName).handles.cur.plot,...
                        'XData', xData,...
                        'YData', yData);
                case 'errorbar'
                    set(this.figs.(gName).handles.cur.plot,...
                        'XData', xData,...
                        'YData', yData,...
                        'YPositiveDelta', stdDev,...
                        'YNegativeDelta', stdDev);
                case 'imagesc'
                    set(this.figs.(gName).handles.cur.plot,...
                        'XData', xData,...
                        'YData', yData,...
                        'CData', cData);
                    this.setLimsToPlot(gName)
            end 
            set(this.figs.(gName).handles.cur.title, 'String', this.figs.(gName).strings.title);
            pause(0.01);
            drawnow();
        end

    end
    
end