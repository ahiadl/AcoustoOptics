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
            figsNames = {'curMainAxis'; 'curMainAxisLog'; 'curMainAxisAvg'; 'curMainPlane'; 'curMainPlaneLog'; 'curMainPlaneAvg'};
       end 
        
       function figs  = createGraphicsVars()    
            % General Scan Parameters 
            figs.scanAxType       = 'normal';
            figs.scanAxis       = [];
            figs.scanAxisCntr   = [];
            figs.scanAxisZero   = [];
            figs.scanAxisNorm   = [];
            
            figs.depthAxType       = 'normal';
            figs.depthAxis       = [];
            figs.depthAxisCntr   = [];
            figs.depthAxisZero   = [];
            figs.depthAxisNorm   = [];
            
            figs.secondPos   = [];
            figs.thirdPos   = [];
            
            figs.repeats     = [];
            figs.useExtClims = false;
            
            figs.scanAxLabel   = 'Y';
            figs.secondAxLabel = 'X';
            figs.thirdAxLabel  = 'Z';
            figs.depthAxLabel   = 'Z';
            figs.mainPlaneLabel = 'YZ';
            
            figs.scanAxType  = 'Normal';
            figs.depthAxType = 'Normal';
            
            figs.scanAxLimsToPlot = [];
            figs.dxFirstAx         = [];
            
            figs.depthIdx    = 1; 
            figs.depthPos    = 0; % the z value related to zIdx;
            
            % Current Position
            figs.curScanPos    = 0;
            figs.curSecondPos  = 0;
            figs.curThirdPos  = 0;
            figs.curScanFirst  = 0;
            figs.curRep        = 0;            
            
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
            figs.scanAxType         = 'Normal';
            figs.depthAxType        = 'Normal';
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
 
            uVars.useFrame   = true;
            uVars.repeats    = 1;
            
            uVars.scanAxisNorm  = [];
            uVars.depthAxisNorm  = [];
            
            uVars.thirdPos  = [];

            uVars.scanAxLabel    = 'Y';
            uVars.secondAxLabel  = 'X';
            uVars.thirdAxLabel   = 'Z';
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
            this.setStrings('curMainAxis', "Current Main Axis %s: (R, %s, %s) = (%d, %.2f, %.2f) - Ch: %d", "Scan %s[mm]", "\\phi[v]", []);
            
            % curMainAxisAvg
            this.setType('curMainAxisAvg', 'errorbar');
            this.setStrings('curMainAxisAvg', "Current Main Axis %s (Averaged): (R, %s, %s) = (%d, %.2f, %.2f) - Ch: %d", "Scan %s[mm]",  "\\phi[v]", []);
            
            % curMainAxis
            this.setType('curMainAxisLog', 'stem');
            this.setStrings('curMainAxisLog', "Current Main Axis (Log) %s: (R, %s, %s) = (%d, %.2f, %.2f) - Ch: %d", "Scan %s[mm]", "log(\\phi)[a.u]", []);
            
            % curMainPlane
            this.setType('curMainPlane', 'imagesc');
            this.setStrings('curMainPlane', "Current Main Plane %s: (R, %s) = (%d, %.2f) - Ch: %d", "Depth (US) %s [mm]", "Scan %s[mm]", []);
            
            % curMainPlaneAvg
            this.setType('curMainPlaneAvg', 'imagesc');
            this.setStrings('curMainPlaneAvg', "Current Main Plane %s (Averaged): (R, %s) = (%d, %.2f) - Ch: %d", "Depth (US) %s [mm]", "Scan %s[mm]", []);
            
            % curMainPlaneLog
            this.setType('curMainPlaneLog', 'imagesc');
            this.setStrings('curMainPlaneLog', "Current Main Plane (Log) %s: (R, %s) = (%d, %.2f)- Ch: %d", "Depth (US) %s [mm] ", "Scan %s[mm]", []);
        end
       
        function setGraphicsScanVars(this)
            this.figs.useFrame    = this.uVars.useFrame;
            this.figs.repeats     = this.uVars.repeats;
            this.figs.useExtClims = this.uVars.useExtClims;
            
            this.figs.channelsInRes         = this.uVars.channelsInRes;
            this.figs.singleChannelAnalysis = this.uVars.singleChannelAnalysis;
            this.setSepChIdx(this.uVars.sepChIdx);

            this.validOwnerGraphics = this.uVars.validOwnerGraphics;
            this.hOwnerGraphics     = this.uVars.hOwnerGraphics;  
            
            this.figs.scanAxLabel    = this.uVars.scanAxLabel;
            this.figs.secondAxLabel  = this.uVars.secondAxLabel;
            this.figs.thirdAxLabel   = this.uVars.thirdAxLabel;
            this.figs.depthAxLabel   = this.uVars.depthAxLabel;
            this.figs.mainPlaneLabel = this.uVars.mainPlane;
            
            this.figs.scanAxType   = this.uVars.scanAxType;
            this.figs.scanAxisNorm = this.uVars.scanAxisNorm;
            this.figs.scanAxisCntr = this.figs.scanAxisNorm - mean(this.figs.scanAxisNorm);
            this.figs.scanAxisZero = abs(this.figs.scanAxisNorm - this.figs.scanAxisNorm(1));
            this.setAxisType(this.figs.scanAxLabel, this.figs.scanAxType);
            
            this.figs.depthAxType   = this.uVars.depthAxType;
            this.figs.depthAxisNorm = this.uVars.depthAxisNorm*1e3;
            this.figs.depthAxisCntr = this.figs.depthAxisNorm - mean(this.figs.depthAxisNorm);
            this.figs.depthAxisZero = abs(this.figs.depthAxisNorm - this.figs.depthAxisNorm(1));
            this.setAxisType(this.figs.depthAxLabel, this.figs.depthAxType);
            
            this.figs.secondPos = this.uVars.secondPos;
            this.figs.thirdPos  = this.uVars.thirdPos;            
            
            if this.uVars.depthIdx > length(this.figs.depthAxis)
                this.figs.depthIdx = 1;
            else
                this.figs.depthIdx = this.uVars.depthIdx;
            end
            
            this.figs.depthPos = this.figs.depthAxis(this.figs.depthIdx);
            
            this.figs.validStruct = this.uVars.validStruct;
            
            if this.figs.repeats == 1
                this.figs.validStruct.curMainAxisAvg  = false;
                this.figs.validStruct.curMainPlaneAvg = false;
            end
            
            this.figs.curRep       = 1;
            this.figs.curScanIdx   = 1;
            this.figs.curSecondIdx = 1;
            
            this.updateGraphicsConstruction();
            this.resetPlots();         
 
            this.figs.fonts = this.uVars.fonts;
        end
        
        function updateCurPosAndIdx(this, pns)
            this.figs.curRep        = pns.curPosIdx(1);
            this.figs.curScanFirst  = pns.curPosIdx(2);
            this.figs.curSecondIdx  = 1;
            this.figs.curScanPos    = pns.curPos(1);
            this.figs.curSecondPos  = pns.curPos(2);
            
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
            if strcmp(ax, this.figs.scanAxLabel)
                this.figs.scanAxType = type;
                switch type
                    case 'Normal'
                        this.figs.scanAxis = this.figs.scanAxisNorm;
                    case 'Center'
                        this.figs.scanAxis = this.figs.scanAxisCntr;
                    case 'Zero'
                        this.figs.scanAxis = this.figs.scanAxisZero;
                    case 'Index'
                        this.figs.scanAxis = 1:length(this.figs.scanAxisNorm);
                end
                this.figs.dxFirstAx         = abs(this.figs.scanAxis(1) - this.figs.scanAxis(2));
                this.figs.scanAxLimsToPlot = [min(this.figs.scanAxis) - this.figs.dxFirstAx,...
                               max(this.figs.scanAxis) + this.figs.dxFirstAx];
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
        
        function setSepChIdx(this, ch)
            if this.figs.singleChannelAnalysis && ch<=this.figs.channelsInRes
                this.figs.sepChIdx = ch;
            else
                this.figs.sepChIdx = 1;
            end
        end
        
        %Set Data functions
        function resetPlots(this)
            this.figs.resetPlots = true;
            this.initDataArray();
            this.initClims();
            this.dispCurMainAxis();
            this.dispCurMainAxisLog();
            this.dispCurMainAxisAvg();
            this.dispCurMainPlane();
            this.dispCurMainPlaneLog()
            this.dispCurMainPlaneAvg();
            this.figs.resetPlots = false;
        end
        
        function resetGraphics(this)
            this.initDataArray()
        end
        
        function initClims(this)
            this.figs.scan2DClims  = zeros(this.figs.channelsInRes, 2);
            this.figs.scan2DClims(:, 1) = 0;
            this.figs.scan2DClims(:, 2) = 1;
            
            this.figs.scanAvg2DClims  = zeros(this.figs.channelsInRes, 2);
            this.figs.scanAvg2DClims(:, 1) = 0;
            this.figs.scanAvg2DClims(:, 2) = 1;
        end
        
        function initDataArray(this)
            rstPhi    = zeros(length(this.figs.scanAxis), 1, length(this.figs.depthAxis), this.uVars.repeats, this.uVars.channelsInRes);
            rstPhiAvg = zeros(length(this.figs.scanAxis), 1, length(this.figs.depthAxis), this.uVars.channelsInRes);
            
            this.data.phi       = rstPhi;
            this.data.phiLog    = rstPhi;
            this.data.phiAvg    = rstPhiAvg;
            this.data.phiAvgStd = rstPhiAvg;
        end
        
        function setLogData(this, phiLog)
            this.data.phiLog = phiLog;
        end
        
        function set1DData(this, data)
            if (this.figs.curScanFirst == 1)
               this.figs.scan2DClims  = zeros(this.figs.channelsInRes, 2);
               this.figs.scan2DClims(:, 1) = inf;
               this.figs.scan2DClims(:, 2) = 0;
            end
            
            for i=1:this.figs.channelsInRes
                this.data.phi(this.figs.curScanFirst, 1, :, this.figs.curRep, i)    = data.phi(:,i);

                curMin = min(this.data.phi(this.figs.curScanFirst, this.figs.curSecondIdx, :, this.figs.curRep, i));
                curMax = max(this.data.phi(this.figs.curScanFirst, this.figs.curSecondIdx, :, this.figs.curRep, i));
                this.figs.scan2DClims(i, 1) = min(this.figs.scan2DClims(i, 1), curMin);
                this.figs.scan2DClims(i, 2) = max(this.figs.scan2DClims(i, 2), curMax);

                if this.figs.scan2DClims(i, 1) < this.figs.extClims(1)
                    this.figs.extClims(1) = this.figs.scan2DClims(i, 1);
                end
                if this.figs.scan2DClims(i, 2) > this.figs.extClims(2)
                    this.figs.extClims(2) = this.figs.scan2DClims(i, 2);
                end
            end
            
            if this.validOwnerGraphics
                this.hOwnerGraphics.set1DData(data);
                this.hOwnerGraphics.updatePlots();
            end
            
        end

        function setAvg1DData(this, data)
            if (this.figs.curRep == 1) && (this.figs.curScanFirst == 1)
               this.figs.scanAvg2DClims  = zeros(this.figs.channelsInRes, 2);
               this.figs.scanAvg2DClims(:, 1) = inf;
               this.figs.scanAvg2DClims(:, 2) = 0;
            end
            
            for i=1:this.figs.channelsInRes
                this.data.phiAvg(this.figs.curScanFirst, this.figs.curSecondIdx, : , i)   = data.phiAvg(:,i);
                this.data.phiAvgStd(this.figs.curScanFirst, this.figs.curSecondIdx, :, i) = data.phiAvgStd(:,i);

                curMin = min(this.data.phiAvg(this.figs.curScanFirst, this.figs.curSecondIdx, :, i));
                curMax = max(this.data.phiAvg(this.figs.curScanFirst, this.figs.curSecondIdx, :, i));
                this.figs.scanAvg2DClims(i, 1) = min(this.figs.scanAvg2DClims(i, 1), curMin);
                this.figs.scanAvg2DClims(i, 2) = max(this.figs.scanAvg2DClims(i, 2), curMax);

                if this.figs.scanAvg2DClims(i, 1) < this.figs.extAvgClims(1)
                    this.figs.extAvgClims(1) = this.figs.scanAvg2DClims(i, 1);
                end
                if this.figs.scanAvg2DClims(i, 2) > this.figs.extAvgClims(2)
                    this.figs.extAvgClims(2) = this.figs.scanAvg2DClims(i, 2);
                end
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
        function data = extractMainAxisFromData(this, type)
            switch type
                case 'phi'
                    data.yData  = this.data.phi(:, this.figs.curSecondIdx, this.figs.depthIdx, this.figs.curRep, this.figs.sepChIdx);
                case 'phiAvg'
                    data.yData  = this.data.phiAvg   (:, this.figs.curSecondIdx, this.figs.depthIdx, this.figs.sepChIdx);
                    data.stdVec = this.data.phiAvgStd(:, this.figs.curSecondIdx, this.figs.depthIdx, this.figs.sepChIdx);    
                case 'phiLog'
                    data.yData = this.data.phiLog(:, this.figs.curSecondIdx, this.figs.depthIdx, this.figs.curRep, this.figs.sepChIdx);
            end
        end
        
        function data = extractMainPlaneFromData(this, type)
            switch type
                case 'phi'
                    if ~this.figs.useExtClims 
                        data.clims = this.figs.scan2DClims(this.figs.sepChIdx,:);
                    else
                        data.clims = this.figs.extClims;
                    end
    %                 data.cData = permute(this.data.phi(:, this.figs.curSecondIdx, :, this.figs.curRep), [3,1,2,4]);
                    data.cData = permute(this.data.phi(:, this.figs.curSecondIdx, :, this.figs.curRep, this.figs.sepChIdx), [1,3,2,4,5]);
                    if data.clims(1) == data.clims(2)
                        data.clims(2) = data.clims(1)+ 0.1;
                    end
                case 'phiAvg'
                    if ~this.figs.useExtClims 
                        data.clims = this.figs.scanAvg2DClims(this.figs.sepChIdx, :);
                    else
                        data.clims = this.figs.extAvgClims;
                    end
                    data.cData = permute(this.data.phiAvg(:, this.figs.curSecondIdx, :, this.figs.sepChIdx), [1,3,2,4]);
                    if data.clims(1) == data.clims(2)
                        data.clims(2) = data.clims(1)+ 0.1;
                    end
                case 'phiLog'
                    data.cData = permute(this.data.phiLog(:, this.figs.curSecondIdx, :, this.figs.curRep, this.figs.sepChIdx), [1,3,2,4,5]);
            end
        end
        
        % Display Functions
        function dispCurMainAxis(this)
            if ~isgraphics(this.figs.curMainAxis.handles.cur.ax)
               return
            end
            
            this.setTitleVariables('curMainAxis',...
                                   {{this.figs.scanAxLabel, this.figs.secondAxLabel, this.figs.depthAxLabel...
                                     this.figs.curRep, this.figs.curSecondPos, this.figs.depthPos, this.figs.sepChIdx}})
            
            %extract relevent data
            plotData = this.extractMainAxisFromData('phi');

            if ~isgraphics(this.figs.curMainAxis.handles.cur.plot) ||...
               this.figs.curMainAxis.update
                
                this.setAxesVar('curMainAxis', this.figs.scanAxLabel, 'Fluence');
                this.setLimits('curMainAxis', this.figs.scanAxLimsToPlot, []);
                
                switch this.figs.curMainAxis.type
                    case 'stem'
                        this.figs.curMainAxis.handles.cur.plot = ...
                            stem(this.figs.curMainAxis.handles.cur.ax,...
                            this.figs.scanAxis, plotData.yData); 
                    case 'errorbar'
                        this.figs.curMainAxis.handles.cur.plot = ...
                            errorbar(this.figs.curMainAxis.handles.cur.ax,...
                            this.figs.scanAxis, plotData.yData, plotData.stdVec);
                end
                
                this.setLimsToPlot('curMainAxis'); %apply limits to the plot
                this.setStringsToPlot('curMainAxis'); %allocate title and labels handles
                drawnow();
                this.figs.curMainAxis.update = false;
            else
                this.quickPlot('curMainAxis', this.figs.scanAxis, plotData.yData, [], [])
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
                                   {{this.figs.scanAxLabel, this.figs.secondAxLabel, this.figs.depthAxLabel,...
                                     this.figs.curRep,       this.figs.curSecondPos,  this.figs.depthPos, this.figs.sepChIdx}});            
            %extract relevent data
            plotData = this.extractMainAxisFromData('phiAvg');
            
            if ~isgraphics(this.figs.curMainAxisAvg.handles.cur.plot) ||...
                this.figs.curMainAxisAvg.update
                
                this.setAxesVar('curMainAxisAvg', this.figs.scanAxLabel, 'Fluence')
                this.setLimits('curMainAxisAvg', this.figs.scanAxLimsToPlot, [])
                switch this.figs.curMainAxisAvg.type
                    case 'stem'
                        this.figs.curMainAxisAvg.handles.cur.plot = ...
                            stem(this.figs.curMainAxisAvg.handles.cur.ax,...
                            this.figs.scanAxis, plotData.yData); 
                    case 'errorbar'
                        this.figs.curMainAxisAvg.handles.cur.plot = ...
                            errorbar(this.figs.curMainAxisAvg.handles.cur.ax,...
                            this.figs.scanAxis, plotData.yData, plotData.stdVec);
                end
                this.setLimsToPlot('curMainAxisAvg')
                this.setStringsToPlot('curMainAxisAvg');
                drawnow();
                this.figs.curMainAxisAvg.update = false;
            else
                this.quickPlot('curMainAxisAvg', this.figs.scanAxis, plotData.yData, [], plotData.stdVec)
            end
        end
        
        function dispCurMainAxisLog(this)
            if ~isgraphics(this.figs.curMainAxisLog.handles.cur.ax)
               return
            end
            
            this.setTitleVariables('curMainAxisLog',...
                                   {{this.figs.scanAxLabel, this.figs.secondAxLabel, this.figs.depthAxLabel...
                                     this.figs.curRep, this.figs.curSecondPos, this.figs.depthPos, this.figs.sepChIdx}})
            
            %extract relevent data
            plotData = this.extractMainAxisFromData('phiLog');

            if ~isgraphics(this.figs.curMainAxisLog.handles.cur.plot) ||...
               this.figs.curMainAxisLog.update
                
                this.setAxesVar('curMainAxisLog', this.figs.scanAxLabel, 'Fluence');
                this.setLimits('curMainAxisLog', this.figs.scanAxLimsToPlot, []);
                
                switch this.figs.curMainAxisLog.type
                    case 'stem'
                        this.figs.curMainAxisLog.handles.cur.plot = ...
                            stem(this.figs.curMainAxisLog.handles.cur.ax,...
                            this.figs.scanAxis, plotData.yData); 
                    case 'errorbar'
                        this.figs.curMainAxisLog.handles.cur.plot = ...
                            errorbar(this.figs.curMainAxisLog.handles.cur.ax,...
                            this.figs.scanAxis, plotData.yData, plotData.stdVec);
                end
                
                this.setLimsToPlot('curMainAxisLog'); %apply limits to the plot
                this.setStringsToPlot('curMainAxisLog'); %allocate title and labels handles
                drawnow();
                this.figs.curMainAxisLog.update = false;
            else
                this.quickPlot('curMainAxisLog', this.figs.scanAxis, plotData.yData, [], [])
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
                                      this.figs.curRep,         this.figs.curSecondPos, this.figs.sepChIdx}});
            
            plotData = this.extractMainPlaneFromData('phi');
                                                  
            this.figs.curMainPlane.lims.clims = plotData.clims;

            if ~isgraphics(this.figs.curMainPlane.handles.cur.plot) ||...
                this.figs.curMainPlane.update
                
                cla(this.figs.curMainPlane.handles.cur.ax)
                this.setAxesVar('curMainPlane', this.figs.depthAxLabel, this.figs.scanAxLabel);
                
                this.figs.curMainPlane.handles.cur.plot = ...
                    imagesc(this.figs.curMainPlane.handles.cur.ax,...
                    'XData', this.figs.depthAxis, 'Ydata', this.figs.scanAxis, 'CData', plotData.cData,...
                    plotData.clims);
                
                %title
                axis(this.figs.curMainPlane.handles.cur.ax, 'tight')
                colorbar(this.figs.curMainPlane.handles.cur.ax);
                this.setLimsToPlot('curMainPlane');
                this.setStringsToPlot('curMainPlane');
                drawnow();
                this.figs.curMainPlane.update = false;
            else
                this.quickPlot('curMainPlane', this.figs.depthAxis, this.figs.scanAxis, plotData.cData)
            end
        end
        
        function dispCurMainPlaneLog(this)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.curMainPlaneLog.handles.cur.ax)
               return
            end
            % "Current Main Plane %s: (R, %s) = (%d, %.2f)"
            this.setTitleVariables('curMainPlaneLog',...
                                   {{ this.figs.mainPlaneLabel, this.figs.secondAxLabel,...
                                      this.figs.curRep,         this.figs.curSecondPos, this.figs.sepChIdx}});
            
            plotData = this.extractMainPlaneFromData('phiLog');

            if ~isgraphics(this.figs.curMainPlaneLog.handles.cur.plot) ||...
                this.figs.curMainPlaneLog.update
                
                cla(this.figs.curMainPlaneLog.handles.cur.ax)
                this.setAxesVar('curMainPlaneLog', this.figs.depthAxLabel, this.figs.scanAxLabel);
                
                this.figs.curMainPlaneLog.handles.cur.plot = ...
                    imagesc(this.figs.curMainPlaneLog.handles.cur.ax,...
                    'XData', this.figs.depthAxis, 'Ydata', this.figs.scanAxis, 'CData', plotData.cData);
                
                %title
                axis(this.figs.curMainPlaneLog.handles.cur.ax, 'tight')
                colorbar(this.figs.curMainPlaneLog.handles.cur.ax);
                this.setLimsToPlot('curMainPlaneLog');
                this.setStringsToPlot('curMainPlaneLog');
                drawnow();
                this.figs.curMainPlaneLog.update = false;
            else
                this.quickPlot('curMainPlaneLog', this.figs.depthAxis, this.figs.scanAxis, plotData.cData)
            end
        end
        
        function dispCurMainPlaneAvg(this)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.curMainPlaneAvg.handles.cur.ax)
               return
            end
            %"Current Main Plane %s (Averaged): (R, %s) = (%d, %.2f)"
            this.setTitleVariables('curMainPlaneAvg', ...
                                   {{ this.figs.mainPlaneLabel, this.figs.secondAxLabel,...
                                      this.figs.curRep,        this.figs.curSecondPos, this.figs.sepChIdx}});
                                                                 
            plotData = this.extractMainPlaneFromData('phiAvg');
            this.figs.curMainPlaneAvg.lims.clims = plotData.clims;

            if ~isgraphics(this.figs.curMainPlaneAvg.handles.cur.plot) ||...
                this.figs.curMainPlane.update
                
                cla(this.figs.curMainPlaneAvg.handles.cur.ax)
                this.setAxesVar('curMainPlaneAvg', this.figs.depthAxLabel, this.figs.scanAxLabel);
                
                this.figs.curMainPlaneAvg.handles.cur.plot = ...
                    imagesc(this.figs.curMainPlaneAvg.handles.cur.ax,...
                    'XData', this.figs.depthAxis, 'Ydata', this.figs.scanAxis , 'CData', plotData.cData,...
                    plotData.clims);
                
                axis(this.figs.curMainPlaneAvg.handles.cur.ax, 'tight')
                colorbar(this.figs.curMainPlaneAvg.handles.cur.ax);
                this.setLimsToPlot('curMainPlaneAvg')
                this.setStringsToPlot('curMainPlaneAvg');
                
                drawnow();
                this.figs.curMainPlaneAvg.update = false;
            else
                this.quickPlot('curMainPlaneAvg', this.figs.depthAxis, this.figs.scanAxis,  plotData.cData)
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