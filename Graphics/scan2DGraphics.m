classdef scan2DGraphics < Graphics
    %SACN2DGRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data
        hOwnerGraphics
        validOwnerGraphics
    end
    
    methods (Static)
  
       function figsNames = getGraphicsNames()
            figsNames = {'curMainAxis'; 'curMainAxisAvg'; 'curMainPlane'; 'curMainPlaneAvg'};
       end 
        
       function figs  = createGraphicsVars()    
            % General Scan Parameters 
            figs.firstAxType       = 'normal';
            figs.firstAxis       = [];
            figs.firstAxisCntr   = [];
            figs.firstAxisZero   = [];
            figs.firstAxisNorm   = [];
            
            figs.depthAxType       = 'normal';
            figs.depthAxis       = [];
            figs.depthAxisCntr   = [];
            figs.depthAxisZero   = [];
            figs.depthAxisNorm   = [];
            
            figs.secondPos   = [];
            
            figs.repeats     = [];
            figs.useExtClims = false;
            
            figs.firstAxLabel   = 'Y';
            figs.secondAxLabel  = 'X';
            figs.depthAxLabel   = 'Z';
            figs.mainPlaneLabel = 'YZ';
            
            figs.firstAxType = 'normal';
            figs.depthAxType = 'normal';
            
            figs.firstAxLimsToPlot = [];
            figs.dxFirstAx         = [];
            
            figs.depthIdx    = 1; 
            figs.depthPos    = 0; % the z value related to zIdx;
            
            % Current Position
            figs.curPosFirst  = 0;
            figs.curPosSecond = 0;
            figs.curIdxFirst  = 0;
            figs.curRep       = 0;            
            
            % Color Limits
            figs.scan2DClims     = [0, 1];
            figs.scanAvg2DClims  = [0, 1];
            figs.extClims        = [0, 1];
            figs.extAvgClims     = [0, 1];
            
            %Specific Figures parameters
            figs.intExt  = 'int';
            figsNames    = scan2DGraphics.getGraphicsNames();
            
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
            figs.depthIdx           = 0;
            figs.firstAxType        = 'normal';
            figs.depthAxType        = 'normal';
            figs.validOwnerGraphics = false; % not realy user vars...
            figs.hOwnerGraphics     = [];
            figs.reopenFigures      = false;
            
            figsNames = scan2DGraphics.getGraphicsNames();
            figs.intExt      = 'int';
            for i=1:length(figsNames)
                figs.validStruct.(figsNames{i}) = true;
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
            uVars = scan2DGraphics.createUserVars();
 
            uVars.useQuant   = true;
            uVars.repeats    = 1;
            
            uVars.firstAxisNorm  = [];
            uVars.depthAxisNorm  = [];
            
            uVars.secondPos  = [];

            uVars.firstAxLabel   = 'Y';
            uVars.secondAxLabel  = 'X';
            uVars.depthAxLabel   = 'Z';
            uVars.mainPlaneLabel = 'YZ';
       end
       
    end
    
    methods
        function this = scan2DGraphics()
            this@Graphics()
            this.validOwnerGraphics = false;
            this.hOwnerGraphics     = [];
            
            this.figsNames = scan2DGraphics.getGraphicsNames();
            this.figs      = scan2DGraphics.createGraphicsVars();
            this.uVars     = scan2DGraphics.createOwnerVars();
            
            this.setGraphicsStaticVars();
            this.numOfFigs = length(this.figsNames);
        end

        function setGraphicsStaticVars(this)
            % curMainAxis
            this.setType('curMainAxis', 'stem');
            this.setStrings('curMainAxis', "Current Main Axis %s: (R, %s, %s) = (%d, %.2f, %.2f)", "%s[mm]", "\\phi[v]", []);
            
            % curMainAxisAvg
            this.setType('curMainAxisAvg', 'errorbar');
            this.setStrings('curMainAxisAvg', "Current Main Axis %s (Averaged): (R, %s, %s) = (%d, %.2f, %.2f)", "%s[mm]",  "\\phi[v]", []);
            
            % curMainPlane
            this.setType('curMainPlane', 'imagesc');
            this.setStrings('curMainPlane', "Current Main Plane %s: (R, %s) = (%d, %.2f)", "%s [mm]", "%s[mm]", []);
            
            % curMainPlaneAvg
            this.setType('curMainPlaneAvg', 'imagesc');
            this.setStrings('curMainPlaneAvg', "Current Main Plane %s (Averaged): (R, %s) = (%d, %.2f)", "%s [mm]", "%s[mm]", []);
        end
       
        function setGraphicsScanVars(this)
            this.figs.useQuant   = this.uVars.useQuant;
            this.figs.repeats    = this.uVars.repeats;
            this.figs.useExtClims = this.uVars.useExtClims;
            
            this.validOwnerGraphics = this.uVars.validOwnerGraphics;
            this.hOwnerGraphics     = this.uVars.hOwnerGraphics;  
            
            this.figs.firstAxLabel   = this.uVars.firstAxLabel;
            this.figs.secondAxLAbel  = this.uVars.secondAxLabel;
            this.figs.depthAxLabel   = this.uVars.depthAxLabel;
            this.figs.mainPlaneLabel = this.uVars.mainPlane;
            
            this.figs.firstAxType   = this.uVars.firstAxType;
            this.figs.firstAxisNorm = this.uVars.firstAxisNorm;
            this.figs.firstAxisCntr = this.figs.firstAxisNorm - mean(this.figs.firstAxisNorm);
            this.figs.firstAxisZero = abs(this.figs.firstAxisNorm - this.figs.firstAxisNorm(1));
            this.setAxisType(this.figs.firstAxLabel, this.figs.firstAxType);
            
            this.figs.depthAxType     = this.uVars.depthAxType;
            this.figs.depthAxisNorm = this.uVars.depthAxisNorm*1e3;
            this.figs.depthAxisCntr = this.figs.depthAxisNorm - mean(this.figs.depthAxisNorm);
            this.figs.depthAxisZero = abs(this.figs.depthAxisNorm - this.figs.depthAxisNorm(1));
            this.setAxisType(this.figs.depthAxLabel, this.figs.depthAxType);
            
            this.figs.secondPos  = this.uVars.secondPos;            
            if this.uVars.depthIdx > length(this.figs.depthAxis)
                this.figs.depthIdx = 1;
            else
                this.figs.depthIdx = this.uVars.depthIdx;
            end
            
            this.figs.depthPos = this.figs.depthAxis(this.figs.depthIdx);
            
            this.figs.validStruct = this.uVars.validStruct;
            if this.figs.repeats > 1
                this.figs.validStruct.curMainAxisAvg  = true;
                this.figs.validStruct.curMainPlaneAvg = true;
            else
                this.figs.validStruct.curMainAxisAvg  = false;
                this.figs.validStruct.curMainPlaneAvg = false;
            end
            
            this.figs.curRep       = 1;
            this.figs.curIdxFirst  = 1;
            this.figs.curIdxSecond = 1;
            
            this.updateGraphicsConstruction();
            this.resetPlots();         

            this.figs.fonts = this.uVars.fonts;
        end
        
        function updateCurPosAndIdx(this, pns)
            this.figs.curRep       = pns.curPosIdx(1);
            this.figs.curIdxFirst  = pns.curPosIdx(2);
            this.figs.curIdxSecond = 1;
            this.figs.curPosFirst  = pns.curPos(1);
            this.figs.curPosSecond = pns.curPos(2);
            if this.validOwnerGraphics
                this.hOwnerGraphics.updateS2DCurPosAndIdx(pns)
            end
        end
        
        function setDepthToMainAxis(this, idx)
            if idx > length(this.figs.depthAxis)
                this.figs.depthIdx = 1;
            else
                this.figs.depthIdx = idx;
            end
            this.figs.depthPos = this.figs.depthAxis(this.figs.depthIdx); 
        end
        
        function setLoadedDataClims(this, data)
            this.figs.scan2DClims(1) = min(data.phi(:));
            this.figs.scan2DClims(2) = max(data.phi(:));
            this.figs.scanAvg2DClims(1) = min(data.phiAvg(:));
            this.figs.scanAvg2DClims(2) = max(data.phiAvg(:));
        end
        
        function setColorsScale(this, useExtClims)
           this.figs.useExtClims = useExtClims; 
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
        function resetPlots(this)
            this.initDataArray();
            this.dispCurMainAxis();
            this.dispCurMainAxisAvg();
            this.dispCurMainPlane();
            this.dispCurMainPlaneAvg();
        end
        
        function resetGraphics(this)
            this.initDataArray()
        end
        
        function initDataArray(this)
            rstPhi    = zeros(length(this.figs.firstAxis), 1, length(this.figs.depthAxis), this.uVars.repeats);
            rstPhiAvg = zeros(length(this.figs.firstAxis), 1, length(this.figs.depthAxis));
            
            this.data.phi       = rstPhi;
            this.data.phiAvg    = rstPhiAvg;
            this.data.phiAvgStd = rstPhiAvg;
        end
        
        function set1DData(this, data)
            this.data.phi(this.figs.curIdxFirst, 1, :, this.figs.curRep)    = data.phi;

            if (this.figs.curIdxFirst == 1)
               this.figs.scan2DClims  = [inf, 0];
            end

            curMin = min(this.data.phi(this.figs.curIdxFirst, this.figs.curIdxSecond, :, this.figs.curRep));
            curMax = max(this.data.phi(this.figs.curIdxFirst, this.figs.curIdxSecond, :, this.figs.curRep));
            this.figs.scan2DClims(1) = min(this.figs.scan2DClims(1), curMin);
            this.figs.scan2DClims(2) = max(this.figs.scan2DClims(2), curMax);
            
            if this.figs.scan2DClims(1) < this.figs.extClims(1)
                this.figs.extClims(1) = this.figs.scan2DClims(1);
            end
            if this.figs.scan2DClims(2) > this.figs.extClims(2)
                this.figs.extClims(2) = this.figs.scan2DClims(2);
            end
            
            if this.validOwnerGraphics
                this.hOwnerGraphics.set1DData(data);
                this.hOwnerGraphics.updatePlots();
            end
        end

        function setAvg1DData(this, data)
            this.data.phiAvg(this.figs.curIdxFirst, this.figs.curIdxSecond, :)    = data.phiAvg;
            this.data.phiAvgStd(this.figs.curIdxFirst, this.figs.curIdxSecond, :) = data.phiAvgStd;
            
            if (this.figs.curRep == 1) && (this.figs.curIdxFirst == 1)
               this.figs.scanAvg2DClims  = [inf, 0];
            end
            
            curMin = min(this.data.phiAvg(this.figs.curIdxFirst, this.figs.curIdxSecond, :));
            curMax = max(this.data.phiAvg(this.figs.curIdxFirst, this.figs.curIdxSecond, :));
            this.figs.scanAvg2DClims(1) = min(this.figs.scanAvg2DClims(1), curMin);
            this.figs.scanAvg2DClims(2) = max(this.figs.scanAvg2DClims(2), curMax);
            
            if this.figs.scanAvg2DClims(1) < this.figs.extAvgClims(1)
                this.figs.extAvgClims(1) = this.figs.scanAvg2DClims(1);
            end
            if this.figs.scanAvg2DClims(2) > this.figs.extAvgClims(2)
                this.figs.extAvgClims(2) = this.figs.scanAvg2DClims(2);
            end
            
            if this.validOwnerGraphics
                this.hOwnerGraphics.setAvg1DData(data);
                this.hOwnerGraphics.updatePlots();
            end
        end

        function setLoadedData(this, data)
            this.data.phi       = data.phi;
            this.data.phiAvg    = data.phiAvg;
            this.data.phiAvgStd = data.phiAvgStd;
            
            this.setLoadedDataClims(data);
        end
        
        function setExternalClims(this, clims, climsAvg)
            this.figs.extClims    = clims;
            this.figs.extAvgClims = climsAvg;
        end
        
        % Extract data functions
        function data = extractMainAxisFromData(this, isAvg)
            if ~isAvg
                data.yData  = this.data.phi(:, this.figs.curIdxSecond, this.figs.depthIdx, this.figs.curRep);
            else
                data.yData  = this.data.phiAvg   (:, this.figs.curIdxSecond, this.figs.depthIdx);
                data.stdVec = this.data.phiAvgStd(:, this.figs.curIdxSecond, this.figs.depthIdx);
            end
        end
        
        function data = extractMainPlaneFromData(this, isAvg)
            if ~isAvg
                if ~this.figs.useExtClims 
                    data.clims = this.figs.scan2DClims;
                else
                    data.clims = this.figs.extClims;
                end
                data.cData = permute(this.data.phi(:, this.figs.curIdxSecond, :, this.figs.curRep), [3,1,2,4]);
            else
                if ~this.figs.useExtClims 
                    data.clims = this.figs.scanAvg2DClims;
                else
                    data.clims = this.figs.extAvgClims;
                end
                data.cData = permute(this.data.phiAvg(:, this.figs.curIdxSecond, :), [3,1,2]);
            end
            if data.clims(1) == data.clims(2)
                data.clims(2) = data.clims(1)+ 0.1;
            end
        end
        
        % Display Functions
        function dispCurMainAxis(this)
            if ~isgraphics(this.figs.curMainAxis.handles.cur.ax)
               return
            end
            
            this.setTitleVariables('curMainAxis',...
                                   {{this.figs.firstAxLabel, this.figs.secondAxLabel, this.figs.depthAxLabel...
                                     this.figs.curRep, this.figs.curPosSecond, this.figs.depthPos}})
            
            %extract relevent data
            plotData = this.extractMainAxisFromData(false);

            if ~isgraphics(this.figs.curMainAxis.handles.cur.plot) ||...
               this.figs.curMainAxis.update
                
                this.setAxesVar('curMainAxis', this.figs.firstAxLabel, 'Fluence');
                this.setLimits('curMainAxis', this.figs.firstAxLimsToPlot, []);
                
                switch this.figs.curMainAxis.type
                    case 'stem'
                        this.figs.curMainAxis.handles.cur.plot = ...
                            stem(this.figs.curMainAxis.handles.cur.ax,...
                            this.figs.firstAxis, plotData.yData); 
                    case 'errorbar'
                        this.figs.curMainAxis.handles.cur.plot = ...
                            errorbar(this.figs.curMainAxis.handles.cur.ax,...
                            this.figs.firstAxis, plotData.yData, plotData.stdVec);
                end
                
                this.setLimsToPlot('curMainAxis'); %apply limits to the plot
                this.setStringsToPlot('curMainAxis'); %allocate title and labels handles
                drawnow();
                this.figs.curMainAxis.update = false;
            else
                quickPlot(this, 'curMainAxis', this.figs.firstAxis, plotData.yData, [], [])
            end
        end
        
        function dispCurMainAxisAvg(this)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.curMainAxisAvg.handles.cur.ax)
               return
            end
            %"Current Main Axis %s (Averaged): (R, %s, %s) = (%d, %.2f, %.2f)"
            this.setTitleVariables('curMainAxisAvg', ...
                                   {{this.figs.firstAxLabel, this.figs.secondAxLabel, this.figs.depthAxLabel,...
                                     this.figs.curRep,       this.figs.curPosSecond,  this.figs.depthPos}});            
            %extract relevent data
            plotData = this.extractMainAxisFromData(true);
            
            if ~isgraphics(this.figs.curMainAxisAvg.handles.cur.plot) ||...
                this.figs.curMainAxisAvg.update
                
                this.setAxesVar('curMainAxisAvg', this.figs.firstAxLabel, 'Fluence')
                this.setLimits('curMainAxisAvg', this.figs.firstAxLimsToPlot, [])
                switch this.figs.curMainAxisAvg.type
                    case 'stem'
                        this.figs.curMainAxisAvg.handles.cur.plot = ...
                            stem(this.figs.curMainAxisAvg.handles.cur.ax,...
                            this.figs.firstAxis, plotData.yData); 
                    case 'errorbar'
                        this.figs.curMainAxisAvg.handles.cur.plot = ...
                            errorbar(this.figs.curMainAxisAvg.handles.cur.ax,...
                            this.figs.firstAxis, plotData.yData, plotData.stdVec);
                end
                this.setLimsToPlot('curMainAxisAvg')
                this.setStringsToPlot('curMainAxisAvg');
                drawnow();
                this.figs.curMainAxisAvg.update = false;
            else
                quickPlot(this, 'curMainAxisAvg', this.figs.firstAxis, plotData.yData, [], plotData.stdVec)
            end
        end
        
        function dispCurMainPlane(this)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.curMainPlane.handles.cur.ax)
               return
            end
            % "Current Main Plane %s: (R, %s) = (%d, %.2f)"
            this.setTitleVariables('curMainPlane',...
                                   {{ this.figs.mainPlaneLabel, this.figs.secondAxLabel,...
                                      this.figs.curRep,         this.figs.curPosSecond}});
            
            plotData = this.extractMainPlaneFromData(false);
                                                  
            this.figs.curMainPlane.lims.clims = plotData.clims;

            if ~isgraphics(this.figs.curMainPlane.handles.cur.plot) ||...
                this.figs.curMainPlane.update
                
                cla(this.figs.curMainPlane.handles.cur.ax)
                this.setAxesVar('curMainPlane', this.figs.firstAxLabel, 'Z');
                
                this.figs.curMainPlane.handles.cur.plot = ...
                    imagesc(this.figs.curMainPlane.handles.cur.ax,...
                    'XData', this.figs.firstAxis, 'Ydata', this.figs.depthAxis, 'CData', plotData.cData,...
                    plotData.clims);
                
                %title
                axis(this.figs.curMainPlane.handles.cur.ax, 'tight')
                colorbar(this.figs.curMainPlane.handles.cur.ax);
                this.setLimsToPlot('curMainPlane');
                this.setStringsToPlot('curMainPlane');
                drawnow();
                this.figs.curMainPlane.update = false;
            else
                quickPlot(this, 'curMainPlane', this.figs.firstAxis, this.figs.depthAxis, plotData.cData)
            end
        end
        
        function dispCurMainPlaneAvg(this)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.curMainPlaneAvg.handles.cur.ax)
               return
            end
            %"Current Main Plane %s (Averaged): (R, %s) = (%d, %.2f)"
            this.setTitleVariables('curMainPlane', ...
                                   {{ this.figs.mainPlaneLabel, this.figs.secondAxLabel,...
                                      this.figs.curRep,        this.figs.curPosSecond}});
                                                                 
            plotData = this.extractMainPlaneFromData(true);
            this.figs.curMainPlaneAvg.lims.clims = plotData.clims;

            if ~isgraphics(this.figs.curMainPlaneAvg.handles.cur.plot) ||...
                this.figs.curMainPlane.update
                
                cla(this.figs.curMainPlaneAvg.handles.cur.ax)
                this.setAxesVar('curMainPlaneAvg', this.figs.firstAxLabel, 'Z');
                
                this.figs.curMainPlaneAvg.handles.cur.plot = ...
                    imagesc(this.figs.curMainPlaneAvg.handles.cur.ax,...
                    'XData', this.figs.firstAxis, 'Ydata', this.figs.depthAxis, 'CData', plotData.cData,...
                    plotData.clims);
                
                axis(this.figs.curMainPlaneAvg.handles.cur.ax, 'tight')
                colorbar(this.figs.curMainPlaneAvg.handles.cur.ax);
                this.setLimsToPlot('curMainPlaneAvg')
                this.setStringsToPlot('curMainPlaneAvg');
                
                drawnow();
                this.figs.curMainPlaneAvg.update = false;
            else
                quickPlot(this, 'curMainPlaneAvg', this.figs.firstAxis, this.figs.depthAxis, plotData.cData)
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