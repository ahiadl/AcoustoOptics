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
            figs.scan1AxType     = 'normal';
            figs.scan1Axis       = [];
            figs.scan1AxisCntr   = [];
            figs.scan1AxisZero   = [];
            figs.scan1AxisNorm   = [];
            
            figs.scan2AxType     = 'normal';
            figs.scan2Axis       = [];
            figs.scan2AxisCntr   = [];
            figs.scan2AxisZero   = [];
            figs.scan2AxisNorm   = [];
            
            figs.depthAxType     = 'normal';
            figs.depthAxis       = [];
            figs.depthAxisCntr   = [];
            figs.depthAxisZero   = [];
            figs.depthAxisNorm   = [];
            
            figs.repeats = [];
            figs.normColorsToPlane = true;
            
            figs.scan1AxLabel   = 'Y';
            figs.scan2AxLabel   = 'X';
            figs.depthAxLabel   = 'Z';
            figs.mainPlaneLabel = 'YZ';
            
            figs.dxscanAx         = [];
            
            % Current Position
            figs.curPosScan1  = 0;
            figs.curPosScan2 = 0;
            figs.curIdxScan1  = 0;
            figs.curIdxScan2 = 0;
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
           
           figs.scan1AxType  = 'normal';
           figs.scan2AxType = 'normal';
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
            
            uVars.scan1AxisNorm   = [];
            uVars.scan2AxisNorm   = [];
            uVars.depthAxisNorm   = [];

            uVars.scan1AxLabel   = 'Y';
            uVars.scan2AxLabel   = 'X';
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
            
            this.figs.scan1AxLabel    = this.uVars.scan1AxLabel;
            this.figs.scan2AxLabel   = this.uVars.scan2AxLabel;
            this.figs.depthAxLabel   = this.uVars.depthAxLabel;
            this.figs.mainPlaneLabel = this.uVars.mainPlane;
            
            this.figs.scan1AxType   = this.uVars.scan1AxType;
            this.figs.scan1AxisNorm = this.uVars.scan1AxisNorm;
            this.figs.scan1AxisCntr = this.figs.scan1AxisNorm - mean(this.figs.scan1AxisNorm);
            this.figs.scan1AxisZero = abs(this.figs.scan1AxisNorm - this.figs.scan1AxisNorm(1));
            this.setAxisType(this.figs.scan1AxLabel, this.figs.scan1AxType);
            
            this.figs.scan2AxType   = this.uVars.scan2AxType;
            this.figs.scan2AxisNorm = this.uVars.scan2AxisNorm;
            this.figs.scan2AxisCntr = this.figs.scan2AxisNorm - mean(this.figs.scan2AxisNorm);
            this.figs.scan2AxisZero = abs(this.figs.scan2AxisNorm - this.figs.scan2AxisNorm(1));
            this.setAxisType(this.figs.scan2AxLabel, this.figs.scan2AxType);
            
            this.figs.depthAxType     = this.uVars.depthAxType;
            this.figs.depthAxisNorm   = this.uVars.depthAxisNorm*1e3;
            this.figs.depthAxisCntr   = this.figs.depthAxisNorm - mean(this.figs.depthAxisNorm);
            this.figs.depthAxisZero   = abs(this.figs.depthAxisNorm - this.figs.depthAxisNorm(1));
            this.setAxisType(this.figs.depthAxLabel, this.figs.depthAxType);
            
            this.figs.normColorsToPlane = this.uVars.normColorsToPlane;
            
            % reset data arrays:
            this.initDataArray();

            this.figs.dxscanAx         = abs(this.figs.scan1Axis(1) - this.figs.scan1Axis(2));
            
            this.figs.validStruct = this.uVars.validStruct;

            %extH and intExt flag should not be transferred at this point.
            %it is updated by the updateGraphicsConstruction function
            this.resetNavigator();
        end
        
        function updateS3DCurPosAndIdx(this, curIdx, curPos)
            this.figs.curIdxScan2 = curIdx;
            this.figs.curPosScan2 = curPos;
        end
        
        function updateS2DCurPosAndIdx(this, pns)
            this.figs.curRep       = pns.curPosIdx(1);
            this.figs.curIdxScan1  = pns.curPosIdx(2);
            this.figs.curPosScan1  = pns.curPos(1);
        end

        function setNavVars(this, navVars)
            this.figs.navIdx =  navVars.navIdx;
            this.figs.navRep =  navVars.navRep;
            if ~strcmp( this.figs.navAx, navVars.navAx) || this.figs.navUpdate
                this.figs.nav.update = true;
                this.figs.navAvg.update = true;
                this.figs.navAx = navVars.navAx;
                if strcmp(this.figs.navAx, this.figs.scan1AxLabel)
                    this.figs.navPlane = sprintf("%s%s", this.figs.scan2AxLabel, this.figs.depthAxLabel);
                    this.figs.navAxis  = this.figs.scan1Axis;
                    this.figs.navXaxis = this.figs.scan2Axis;
                    this.figs.navYaxis = this.figs.depthAxis;
                    this.figs.navXname = this.figs.scan2AxLabel;
                    this.figs.navYname = this.figs.depthAxLabel;
                elseif strcmp(this.figs.navAx, this.figs.scan2AxLabel)
                    this.figs.navPlane = sprintf("%s%s", this.figs.scan1AxLabel, this.figs.depthAxLabel);
                    this.figs.navAxis  = this.figs.scan2Axis;
                    this.figs.navXaxis = this.figs.scan1Axis;
                    this.figs.navYaxis = this.figs.depthAxis;
                    this.figs.navXname = this.figs.scan1AxLabel;
                    this.figs.navYname = this.figs.depthAxLabel;
                elseif strcmp(this.figs.navAx, this.figs.depthAxLabel)
                    this.figs.navPlane = sprintf("%s%s", this.figs.scan1AxLabel, this.figs.scan2AxLabel);
                    this.figs.navAxis  = this.figs.depthAxis;
                    this.figs.navXaxis = this.figs.scan1Axis;
                    this.figs.navYaxis = this.figs.scan2Axis;
                    this.figs.navXname = this.figs.scan1AxLabel;
                    this.figs.navYname = this.figs.scan2AxLabel;
                else
                    fprintf("S3D Graphics: ERROR: no compatible navigator Axis");
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
            if strcmp(ax, this.figs.scan1AxLabel)
                this.figs.scan1AxType = type;
                switch type
                    case 'Normal'
                        this.figs.scan1Axis = this.figs.scan1AxisNorm;
                    case 'Center'
                        this.figs.scan1Axis = this.figs.scan1AxisCntr;
                    case 'Zero'
                        this.figs.scan1Axis = this.figs.scan1AxisZero;
                    case 'Index'
                        this.figs.scan1Axis = 1:length(this.figs.scan1AxisNorm);
                end
                this.figs.dxscanAx         = abs(this.figs.scan1Axis(1) - this.figs.scan1Axis(2));
                this.figs.scanAxLimsToPlot = [min(this.figs.scan1Axis) - this.figs.dxscanAx,...
                               max(this.figs.scan1Axis) + this.figs.dxscanAx];
            elseif strcmp(ax, this.figs.scan2AxLabel)
                this.figs.depthAxType = type;
                switch type
                    case 'Normal'
                        this.figs.scan2Axis = this.figs.scan2AxisNorm;
                    case 'Center'
                        this.figs.scan2Axis = this.figs.scan2AxisCntr;
                    case 'Zero'
                        this.figs.scan2Axis = this.figs.scan2AxisZero;
                    case 'Index' 
                        this.figs.scan2Axis = 1:length(this.figs.scan2AxisNorm);
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
                        this.figs.depthAxis = 1:length(this.figs.depthAxisNorm);
                end
            end  
        end
        
        %Set Data functions
        function initDataArray(this)
            rstPhi    = zeros(length(this.figs.scan1Axis), length(this.figs.scan2Axis), length(this.figs.depthAxis), this.uVars.repeats);
            rstPhiAvg = zeros(length(this.figs.scan1Axis), length(this.figs.scan2Axis), length(this.figs.depthAxis));
            
            this.data.phi       = rstPhi;
            this.data.phiAvg    = rstPhiAvg;
            this.data.phiAVgStd = rstPhiAvg;
        end
        
        function set1DData(this, data)
            this.data.phi(this.figs.curIdxScan1, this.figs.curIdxScan2, :, this.figs.curRep) = data.phi;

            if (this.figs.curIdxScan1 == 1)
               this.figs.scan2DClims  = [inf, 0];
            end
            
            if (this.figs.curIdxScan1 == 1) && (this.figs.curIdxScan2 == 1) && (this.figs.curRep == 1)
                this.figs.scan3DClims  = [inf, 0];
            end    

            curMin = min(this.data.phi(this.figs.curIdxScan1, this.figs.curIdxScan2, :, this.figs.curRep));
            curMax = max(this.data.phi(this.figs.curIdxScan1, this.figs.curIdxScan2, :, this.figs.curRep));
            this.figs.scan2DClims(1) = min(this.figs.scan2DClims(1), curMin);
            this.figs.scan2DClims(2) = max(this.figs.scan2DClims(2), curMax);
            this.figs.scan3DClims(1) = min(this.figs.scan3DClims(1), curMin);
            this.figs.scan3DClims(2) = max(this.figs.scan3DClims(2), curMax);
        end

        function setAvg1DData(this, data)
            this.data.phiAvg(this.figs.curIdxScan1, this.figs.curIdxScan2, :)     = data.phiAvg;
            this.data.phiAvgStd(this.figs.curIdxScan1, this.figs.curIdxScan2, :) = data.phiAvgStd;
            
            if (this.figs.curRep == 1) && (this.figs.curIdxScan1 == 1)
               this.figs.scanAvg2DClims  = [inf, 0];
            end
            
            if (this.figs.curIdxScan1 == 1) && (this.figs.curIdxScan2 == 1) && (this.figs.curRep == 1)
                this.figs.scanAvg3DClims  = [inf, 0];
            end 
            
            curMin = min(this.data.phiAvg(this.figs.curIdxScan1, this.figs.curIdxScan2, :));
            curMax = max(this.data.phiAvg(this.figs.curIdxScan1, this.figs.curIdxScan2, :));
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
                if strcmp(this.figs.navAx, this.figs.scan1AxLabel)
                    data.cData = permute(this.data.phi(this.figs.navIdx, :, :, this.figs.navRep), [3,2,1,4]);
                elseif strcmp(this.figs.navAx, this.figs.scan2AxLabel)
                    data.cData  = permute(this.data.phi(:, this.figs.navIdx, :, this.figs.navRep), [3,1,2,4]);
                elseif strcmp(this.figs.navAx, this.figs.depthAxLabel)
                    data.cData  = permute(this.data.phi(:, :, this.figs.navIdx, this.figs.navRep), [2,1,3,4]); 
                else
                    fprintf("S2D ERROR: no compatible navigator axis");
                end
            else
%                 data.clims = [min(min(min(min(this.data.phiAvg)))), max(max(max(max(this.data.phiAvg))))];
                if strcmp(this.figs.navAx, this.figs.scan1AxLabel)
                    data.cData  = permute(this.data.phiAvg(this.figs.navIdx, :, :), [3,2,1]);
                elseif strcmp(this.figs.navAx, this.figs.scan2AxLabel)
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