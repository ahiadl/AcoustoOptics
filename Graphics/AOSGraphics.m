classdef AOSGraphics < Graphics
    %SACN2DGRAPHICS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        data
        hOwnerGraphics
        validOwnerGraphics
    end

    methods (Static)
       function figsNames = getGraphicsNames()
            figsNames = {'scanPlane'; 'scanPlaneLog'; "navLine"; "navLineLog"; "navPlane"; "navPlaneLog"};
       end 
       
       function vars = getLineNavVarsStruct()
           vars.ax1Label = [];
           vars.ax1Idx   = [];
           vars.ax2Label = [];
           vars.ax2Idx   = []; 
       end
        
       function vars = getPlaneNavVarsStruct()
           vars.normAx = [];
           vars.axIdx  = []; 
       end

       function figs  = createGraphicsVars()    
            % General Scan Parameters 
            figs.scan1AxType       = 'Normal';
            figs.scan1Vec       = [];
            figs.scan1VecCntr   = [];
            figs.scan1VecZero   = [];
            figs.scan1VecNorm   = [];
            
            figs.scan2AxType       = 'Normal';
            figs.scan2Vec       = [];
            figs.scan2VecCntr   = [];
            figs.scan2VecZero   = [];
            figs.scan2VecNorm   = [];

            figs.depthAxType       = 'Normal';
            figs.depthVec       = [];
            figs.depthVecCntr   = [];
            figs.depthVecZero   = [];
            figs.depthVecNorm   = [];
            
            figs.thirdPos   = [];

            figs.scan1Label     = 'Y';
            figs.scan2Label     = 'X';
            figs.thirdLabel     = 'Z';
            figs.depthLabel     = 'Z';
            figs.mainPlaneLabel = 'YZ';

            figs.scanAxLimsToPlot = [];
            figs.dScan1           = [];
            figs.dScan2           = [];
            figs.dDDepth          = [];
 
            % Current Position
            figs.scan1Pos  = 0;
            figs.scan2Pos  = 0;
            figs.thirdPos  = 0;         
            
            figs.scan1Idx = 1;
            figs.scan2Idx = 1;

            % Line Navigator Position          
            figs.varsNavLine  = AOSGraphics.getLineNavVarsStruct;
            figs.varsNavPlane = AOSGraphics.getPlaneNavVarsStruct;
            
            figs.scan1LimsToPlot = [];
            figs.scan2LimsToPlot = [];
            figs.depthLimsToPlot = [];
            
            figs.normDispPlaneColors = false;

            %Specific Figures parameters
            figs.intExt  = 'int';
            figsNames    = AOSGraphics.getGraphicsNames();
            
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
            figs.scan1AxType        = 'Normal';
            figs.scan2AxType        = 'Normal';
            figs.depthAxType        = 'Normal';

            figs.reopenFigures       = false;
            figs.normDispPlaneColors = false;

            figsNames = AOSGraphics.getGraphicsNames();
            figs.intExt      = 'int';
            for i=1:length(figsNames)
                figs.validStruct.(figsNames{i}) = true;
            end
            
            for i=1:length(figsNames)
                figs.extH.(figsNames{i}) =  Graphics.createHandlesStruct();
            end
            
            figs.varsNavLine  = AOSGraphics.getLineNavVarsStruct;
            figs.varsNavPlane = AOSGraphics.getPlaneNavVarsStruct;

            figs.fonts.type       = [];
            figs.fonts.titleSize  = 14;
            figs.fonts.labelsSize = 14;
            figs.fonts.axisSize   = 14;
        end

       function uVars = createOwnerVars() %for Reduced Operation
            uVars = AOSGraphics.createUserVars();

            uVars.scan1VecNorm  = [];
            uVars.scan2VecNorm  = [];
            uVars.depthVecNorm  = [];
            
            uVars.thirdPos  = [];

            uVars.scan1Label   = 'Y';
            uVars.scan2Label   = 'X';
            uVars.thirdLabel   = 'Z';
            uVars.depthLabel   = 'Z';
            uVars.mainPlaneLabel = 'YZ';
       end
    end

    methods
        function this = AOSGraphics()
            this@Graphics()
            this.validOwnerGraphics = false;
            this.hOwnerGraphics     = [];
            
            this.figsNames = AOSGraphics.getGraphicsNames();
            this.figs      = AOSGraphics.createGraphicsVars();
            this.uVars     = AOSGraphics.createOwnerVars();
            
            this.setGraphicsStaticVars();
            this.numOfFigs = length(this.figsNames);
        end

        function setGraphicsStaticVars(this)
            % scanPlane
            this.setType('scanPlane', 'imagesc');
            this.setStrings('scanPlane', "Scan Plane %s: %s = %.2f [mm] - Ch: %d", "Depth (US) %s [mm]", "Scan %s[mm]", []);
            
            % scanPlaneLog
            this.setType('scanPlaneLog', 'imagesc');
            this.setStrings('scanPlaneLog', "Scan Plane (Log) %s: %s = %.2f - Ch: %d", "Depth (US) %s [mm] ", "Scan %s[mm]", []);
            
            % navLine
            this.setType('navLine', 'stem');
            this.setStrings('navLine', "Line Navigator:  (%s, %s) = (%.2f, %.2f) - Ch: %d", "%s [mm]", "\\phi[v]", []);
            
            % navLineLog
            this.setType('navLineLog', 'plot');
            this.setStrings('navLineLog', "Line Navigator (Log): (%s, %s) = (%.2f, %.2f) - Ch: %d", "%s [mm]",  "\\phi[Log]", []);
            
             % navPlane
            this.setType('navPlane', 'imagesc');
            this.setStrings('navPlane', "Plane Navigator %s: %s = %.2f - Ch: %d", "%s[mm]", "%s[mm]", []);
            
            % navPlaneLog
            this.setType('navPlaneLog', 'imagesc');
            this.setStrings('navPlaneLog', "Plane Navigator %s (Log):  %s =  %.2f - Ch: %d", "%s[mm]", "%s[mm]", []);
        end

        function setGraphicsScanVars(this)           
            this.figs.channelsInRes         = this.uVars.channelsInRes;
            this.figs.singleChannelAnalysis = this.uVars.singleChannelAnalysis;
            this.setSepChIdx(this.uVars.sepChIdx); 
            this.figs.normDispPlaneColors = this.uVars.normDispPlaneColors;

            this.figs.scan1Label   = this.uVars.scan1Label;
            this.figs.scan2Label   = this.uVars.scan2Label;
            this.figs.thirdLabel   = this.uVars.thirdLabel;
            this.figs.depthLabel   = this.uVars.depthLabel;
            this.figs.mainPlaneLabel = this.uVars.mainPlane;
            
            this.figs.scan1AxType   = this.uVars.scan1AxType;
            this.figs.scan2AxType   = this.uVars.scan1AxType;
            this.figs.depthAxType   = this.uVars.depthAxType;

            this.figs.scan1VecNorm = this.uVars.scan1VecNorm;
            this.figs.scan1VecCntr = this.figs.scan1VecNorm - mean(this.figs.scan1VecNorm);
            this.figs.scan1VecZero = abs(this.figs.scan1VecNorm - this.figs.scan1VecNorm(1));
            this.setAxisType(this.figs.scan1Label);
            
            this.figs.scan2VecNorm = this.uVars.scan2VecNorm;
            this.figs.scan2VecCntr = this.figs.scan2VecNorm - mean(this.figs.scan2VecNorm);
            this.figs.scan2VecZero = abs(this.figs.scan2VecNorm - this.figs.scan2VecNorm(1));
            this.setAxisType(this.figs.scan2Label);

            this.figs.depthVecNorm = this.uVars.depthVecNorm*1e3;
            this.figs.depthVecCntr = this.figs.depthVecNorm - mean(this.figs.depthVecNorm);
            this.figs.depthVecZero = abs(this.figs.depthVecNorm - this.figs.depthVecNorm(1));
            this.setAxisType(this.figs.depthLabel);
            
            this.figs.thirdPos  = this.uVars.thirdPos;            
            
            this.figs.scan1Idx = 1;
            this.figs.scan2Idx = 1;
            this.figs.scan1Pos = 0;
            this.figs.scan2Pos = 0;

            this.figs.varsNavLine.ax1Label = this.uVars.varsNavLine.ax1Label;
            this.figs.varsNavLine.ax1Idx   = this.uVars.varsNavLine.ax1Idx;
            this.figs.varsNavLine.ax2Label = this.uVars.varsNavLine.ax2Label;
            this.figs.varsNavLine.ax2Idx   = this.uVars.varsNavLine.ax2Idx; 

            this.figs.varsNavPlane.normAx = this.uVars.varsNavPlane.normAx;
            this.figs.varsNavPlane.axIdx  = this.uVars.varsNavPlane.axIdx; 

            this.figs.validStruct = this.uVars.validStruct;

            this.updateGraphicsConstruction();
            this.resetPlots();         
        
            this.figs.fonts = this.uVars.fonts;
        end

        function updateCurPosAndIdx(this, pns)
            this.figs.scan1Idx  = pns.curIdx(1);
            this.figs.scan2Idx  = pns.curIdx(2);
            this.figs.scan1Pos  = pns.curPos(1);
            this.figs.scan2Pos  = pns.curPos(2);
        end

        function setAxisType(this, ax)
            if strcmp(ax, this.figs.scan1Label)
                switch this.figs.scan1AxType
                    case 'Normal'
                        this.figs.scan1Vec = this.figs.scan1VecNorm;
                    case 'Center'
                        this.figs.scan1Vec = this.figs.scan1VecCntr;
                    case 'Zero'
                        this.figs.scan1Vec = this.figs.scan1VecZero;
                    case 'Index'
                        this.figs.scan1Vec = 1:length(this.figs.scan1VecNorm);
                end
                this.figs.dScan1         = abs(this.figs.scan1Vec(1) - this.figs.scan1Vec(2));
                this.figs.scan1LimsToPlot = [min(this.figs.scan1Vec) - this.figs.dScan1,...
                                             max(this.figs.scan1Vec) + this.figs.dScan1];
            elseif strcmp(ax, this.figs.scan2Label)
                switch this.figs.scan2AxType
                    case 'Normal'
                        this.figs.scan2Vec = this.figs.scan2VecNorm;
                    case 'Center'
                        this.figs.scan2Vec = this.figs.scan2VecCntr;
                    case 'Zero'
                        this.figs.scan2Vec = this.figs.scan2VecZero;
                    case 'Index'
                        this.figs.scan2Vec = 1:length(this.figs.scan2VecNorm);
                end
                if length(this.figs.scan2Vec(1)) == 1
                    this.figs.dScan2 = 1;
                else
                    this.figs.dScan2          = abs(this.figs.scan2Vec(1) - this.figs.scan2Vec(2));
                end
                this.figs.scan2LimsToPlot = [min(this.figs.scan2Vec) - this.figs.dScan2,...
                                             max(this.figs.scan2Vec) + this.figs.dScan2];
            elseif strcmp(ax, this.figs.depthLabel)
                switch this.figs.depthAxType
                    case 'Normal'
                        this.figs.depthVec = this.figs.depthVecNorm;
                    case 'Center'
                        this.figs.depthVec = this.figs.depthVecCntr;
                    case 'Zero'
                        this.figs.depthVec = this.figs.depthVecZero;
                    case 'Index' 
                        this.figs.depthVec = 1:length(this.figs.depthVecNorm);
                end
                this.figs.dDepth = abs(this.figs.depthVec(1) - this.figs.depthVec(2));
                this.figs.depthLimsToPlot = [min(this.figs.depthVec) - this.figs.dDepth,...
                                             max(this.figs.depthVec) + this.figs.dDepth];
            end  
        end

        function setSepChIdx(this, ch)
            if this.figs.singleChannelAnalysis && ch<=this.figs.channelsInRes
                this.figs.sepChIdx = ch;
            else
                this.figs.sepChIdx = 1;
            end
        end

        function setNavLineVars(this, vars)
            this.figs.varsNavLine.ax1Label = vars.ax1Label;
            this.figs.varsNavLine.ax1Idx   = vars.ax1Idx;
            this.figs.varsNavLine.ax2Label = vars.ax2Label;
            this.figs.varsNavLine.ax2Idx   = vars.ax2Idx;
        end

        function setNavPlaneVars(this, vars)
            this.figs.varsNavPlane.normAx  = vars.normAx;
            this.figs.varsNavPlane.axIdx   = vars.axIdx;
        end
    
        function setDispPlaneColorNorm(this, val)
            this.figs.normDispPlaneColors = val;
        end

        %Set Data functions
        function resetPlots(this)
            this.figs.resetPlots = true;
            this.initDataArray();
            this.dispAll();
            this.figs.resetPlots = false;
        end

        function resetLineNav(this)
            navVars.navAx = this.figs.depthAxLabel;
            navVars.navIdx = 1;
            navVars.navRep = 1;
            this.figs.navUpdate = true;
            this.setLineNavVars(navVars);
            this.dispLineNav("lineNav");
            this.dispLineNav("lineNavLog");
        end

        function resetGraphics(this)
            this.initDataArray()
        end

        function initDataArray(this)
            rstPhi    = zeros(length(this.figs.scan1Vec), length(this.figs.scan2Vec), length(this.figs.depthVec),  this.uVars.channelsInRes);
            this.data.phi       = rstPhi;
            this.data.phiLog    = rstPhi;
            this.data.phiNorm   = rstPhi;

            this.data.glbMax  = ones(1, this.uVars.channelsInRes);
            this.data.glbMin  = zeros(1, this.uVars.channelsInRes);
            this.data.glbSpan = ones(1, this.uVars.channelsInRes);

            this.data.glbMinLog = -Inf*ones(1, this.uVars.channelsInRes);;
        end

        function setData(this, res)
            this.data.phi     = res.phi;
            this.data.phiLog  = res.phiLog;
            this.data.phiNorm = res.phiNorm;

            this.data.glbMax  = res.globalMax;
            this.data.glbMin  = res.globalMin;
            this.data.glbSpan = res.globalSpan;

            this.data.glbMinLog = res.globalMinLog;
        end

        function setLoadedData(this, data)
            this.data.phi       = data.phi;
            this.data.phiAvg    = data.phiAvg;
            this.data.phiAvgStd = data.phiAvgStd;
            
            this.setLoadedDataClims(data);
        end

        % Extract data functions      
        function dim = extractDimFromLabel(this, label)
            if strcmp(label, this.figs.scan1Label)
                dim = 1;
            elseif strcmp(label, this.figs.scan2Label)
                dim = 2;
            elseif strcmp(label, this.figs.depthLabel)
                dim = 3;
            end
        end
        
        function [ax1, ax2, lineVec, lineLabel] = extractLineNavAxFromLabels(this, label1, label2)
            if strcmp(label1, this.figs.scan1Label)
                ax1 = this.figs.scan1Vec;
                if strcmp(label2, this.figs.scan1Label)
                    error("Identical Line Navigator Labels");
                elseif strcmp(label2, this.figs.scan2Label)
                    ax2 = this.figs.scan2Vec;
                    lineVec = this.figs.depthVec;
                    lineLabel = this.figs.depthLabel;
                elseif strcmp(label2, this.figs.depthLabel)
                    ax2 = this.figs.depthVec;
                    lineVec = this.figs.scan2Vec;
                    lineLabel = this.figs.scan2Label;
                end
            elseif strcmp(label1, this.figs.scan2Label)
                ax1 = this.figs.scan2Vec;
                if strcmp(label2, this.figs.scan1Label)
                    ax2 = this.figs.scan1Vec;
                    lineVec = this.figs.depthVec;
                    lineLabel = this.figs.depthLabel;
                elseif strcmp(label2, this.figs.scan2Label)
                    error("Identical Line Navigator Labels");
                elseif strcmp(label2, this.figs.depthLabel)
                    ax2 = this.figs.depthVec;
                    lineVec = this.figs.scan1Vec;
                    lineLabel = this.figs.scan1Label;
                end
            elseif strcmp(label1, this.figs.depthLabel)
                ax1 = this.figs.depthVec;
                if strcmp(label2, this.figs.scan1Label)
                    ax2 = this.figs.scan1Vec;
                    lineVec = this.figs.scan2Vec;
                    lineLabel = this.figs.scan2Label;
                elseif strcmp(label2, this.figs.scan2Label)
                    ax2 = this.figs.scan2Vec;
                    lineVec = this.figs.scan1Vec;
                    lineLabel = this.figs.scan1Label;
                elseif strcmp(label2, this.figs.depthLabel)
                    error("Identical Line Navigator Labels");
                end
            end
        end
        
        function [ax1, ax2, ax3, label1, label2, planeLabel] = extractPlaneNavAxes(this, ax3Label)
            if strcmp(ax3Label, this.figs.scan1Label)
                ax1 = this.figs.depthVec;
                ax2 = this.figs.scan2Vec;
                ax3 = this.figs.scan1Vec;
                label1 = this.figs.depthLabel;
                label2 = this.figs.scan2Label;
            elseif strcmp(ax3Label, this.figs.scan2Label)
                ax1 = this.figs.depthVec;
                ax2 = this.figs.scan1Vec;
                ax3 = this.figs.scan2Vec;
                label1 = this.figs.depthLabel;
                label2 = this.figs.scan1Label;
            elseif strcmp(ax3Label, this.figs.depthLabel)
                ax1 = this.figs.scan1Vec;
                ax2 = this.figs.scan2Vec;
                ax3 = this.figs.depthVec;
                label1 = this.figs.scan1Label;
                label2 = this.figs.scan2Label;
            end
            planeLabel = sprintf("%s%s", label1, label2);
        end

        function data = extractScanPlaneFromData(this, gName)
            data.xData = this.figs.depthVec;
            data.yData = this.figs.scan1Vec;

            switch gName
                case 'scanPlane'
                    data.cData = permute(this.data.phiNorm(:, this.figs.scan2Idx, :, this.figs.sepChIdx), [1,3,2,4]);;
                    if this.figs.normDispPlaneColors
                        data.clims = [min(data.cData(:)), max(data.cData(:))];
                    else
                        data.clims = [0, 1];
                    end
                case 'scanPlaneLog'
                    data.cData = permute(this.data.phiLog(:, this.figs.scan2Idx, :, this.figs.sepChIdx), [1,3,2,4]);
                    if this.figs.normDispPlaneColors
                        data.clims = [min(data.cData(:)), 0];
                    else
                        data.clims = [this.data.glbMinLog(this.figs.sepChIdx), 0];
                    end
            end
            if data.clims(1) == data.clims(2)
                data.clims(2)= data.clims(1)+1;
            end
            data.titleVars = {{ this.figs.mainPlaneLabel,...
                                this.figs.scan2Label,...
                                this.figs.scan2Pos,...
                                this.figs.sepChIdx}};
            data.xLabelVar = this.figs.depthLabel;
            data.yLabelVar = this.figs.scan1Label;
        end
        
        function data = extractNavPlaneFromData(this, gName)
            [ax1, ax2, axNorm, label1, label2, planeLabel] = this.extractPlaneNavAxes(this.figs.varsNavPlane.normAx);
            data.xData = ax2;
            data.yData = ax1;
            
            dimNorm = this.extractDimFromLabel(this.figs.varsNavPlane.normAx);
            dimOrder = 1:length(size(this.data.phi));
            dimOrder(dimNorm) = [];
            dimOrder = [dimNorm,dimOrder];
            
            switch gName
                case 'navPlane'
                    phi = permute(this.data.phiNorm, dimOrder); % Bring normal axis to dim 1;
                    phi = squeeze(phi(this.figs.varsNavPlane.axIdx, :, :, this.figs.sepChIdx));
                    data.cData = phi;
                    if this.figs.normDispPlaneColors
                        data.clims = [min(phi(:)), max(phi(:))];
                    else
                        data.clims = [0, 1];
                    end
                case 'navPlaneLog'
                    phiLog     = permute(this.data.phiLog, dimOrder);
                    data.cData = squeeze(phiLog(this.figs.varsNavPlane.axIdx, :, :, this.figs.sepChIdx));
                    if this.figs.normDispPlaneColors
                        data.clims = [min(data.cData(:)), 0];
                    else
                        data.clims = [this.data.glbMinLog(this.figs.sepChIdx), 0];
                    end
            end
            if data.clims(1) == data.clims(2)
                data.clims(2)= data.clims(1)+1;
            end
            data.titleVars = {{ planeLabel,...
                                this.figs.varsNavPlane.normAx,...
                                axNorm(this.figs.varsNavPlane.axIdx),...
                                this.figs.sepChIdx}};
            data.xLabelVar = label2;
            data.yLabelVar = label1;
        end

        function data = extractNavLineFromData(this, gName)
            dim1 = this.extractDimFromLabel(this.figs.varsNavLine.ax1Label);
            dim2 = this.extractDimFromLabel(this.figs.varsNavLine.ax2Label);
            
            dimOrder = 1:length(size(this.data.phi));
            dimOrder([dim1, dim2]) = [];
            dimOrder = [dim1, dim2, dimOrder];

            [ax1, ax2, data.xData, lineLabel] = this.extractLineNavAxFromLabels(this.figs.varsNavLine.ax1Label,...
                                                                                this.figs.varsNavLine.ax2Label);
            switch gName
                case "navLine"
                    phi = permute(this.data.phiNorm, dimOrder);
                case "navLineLog"
                    phi = permute(this.data.phiLog, dimOrder);
            end
            data.yData = squeeze(phi(this.figs.varsNavLine.ax1Idx, this.figs.varsNavLine.ax2Idx, :,  this.figs.sepChIdx));
            
            data.titleVars  =  {{this.figs.varsNavLine.ax1Label,...
                                 this.figs.varsNavLine.ax2Label,...
                                 ax1(this.figs.varsNavLine.ax1Idx),...
                                 ax2(this.figs.varsNavLine.ax2Idx),...
                                 this.figs.sepChIdx}};

            data.xLabelVar = lineLabel;
            data.yLabelVar = [];
        end

        % Display Functions
        function dispScanPlanes(this)
            this.dispScanPlane("scanPlane");
            this.dispScanPlane("scanPlaneLog");
        end

        function dispLineNavAll(this)
            this.dispNavLine("navLine");
            this.dispNavLine("navLineLog");
        end

        function dispPlaneNavAll(this)
            this.dispNavPlane("navPlane");
            this.dispNavPlane("navPlaneLog");
        end

        function dispAll(this)
            this.dispScanPlanes();
            this.dispLineNavAll();
            this.dispPlaneNavAll();
        end

        function dispScanPlane(this, gName)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.(gName).handles.cur.ax)
               return
            end
            % "Current Main Plane %s: (R, %s) = (%d, %.2f)"
            plotData = this.extractScanPlaneFromData(gName);

            this.setTitleVariables(gName, plotData.titleVars)...                    

            if ~isgraphics(this.figs.(gName).handles.cur.plot) ||...
                this.figs.(gName).update
                
                cla(this.figs.(gName).handles.cur.ax)
                this.setAxesVar(gName, plotData.xLabelVar, plotData.yLabelVar);
                
                this.figs.(gName).handles.cur.plot = ...
                    imagesc(this.figs.(gName).handles.cur.ax,...
                            'XData', plotData.xData,...
                            'YData', plotData.yData,...
                            'CData', plotData.cData,...
                            plotData.clims);
                
                %title
                axis(this.figs.(gName).handles.cur.ax, 'tight')
                colorbar(this.figs.(gName).handles.cur.ax);
                this.setLimsToPlot(gName);
                this.setStringsToPlot(gName);
                drawnow();
                this.figs.(gName).update = false;
            else
                this.quickPlot(gName, plotData.xData, plotData.yData, plotData.cData, plotData.clims)
            end
        end
        
        function dispNavPlane(this, gName)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.(gName).handles.cur.ax)
               return
            end
            % "Current Main Plane %s: (R, %s) = (%d, %.2f)"
            plotData = this.extractNavPlaneFromData(gName);

            this.setTitleVariables(gName, plotData.titleVars)                   

            if ~isgraphics(this.figs.(gName).handles.cur.plot) ||...
                this.figs.(gName).update
                
                cla(this.figs.(gName).handles.cur.ax)
                this.setAxesVar(gName, plotData.xLabelVar, plotData.yLabelVar);
                
                this.figs.(gName).handles.cur.plot = ...
                    imagesc(this.figs.(gName).handles.cur.ax,...
                    'XData', plotData.xData,...
                    'YData', plotData.yData,...
                    'CData', plotData.cData, ...
                    plotData.clims);
                
                %title
                axis(this.figs.(gName).handles.cur.ax, 'tight')
                colorbar(this.figs.(gName).handles.cur.ax);
                this.setLimsToPlot(gName);
                this.setStringsToPlot(gName);
                drawnow();
                this.figs.(gName).update = false;
            else
                this.quickPlot(gName, plotData.xData, plotData.yData, plotData.cData, plotData.clims)
            end
        end

        function dispNavLine(this, gName)
            if ~isgraphics(this.figs.(gName).handles.cur.ax)
               return
            end
            
            %extract relevent data
            plotData = this.extractNavLineFromData(gName);

            this.setTitleVariables(gName,plotData.titleVars)

            if ~isgraphics(this.figs.(gName).handles.cur.plot) || this.figs.(gName).update
                this.setAxesVar(gName, plotData.xLabelVar, plotData.yLabelVar);

                switch this.figs.(gName).type
                    case 'stem'
                        this.figs.(gName).handles.cur.plot = ...
                            stem(this.figs.(gName).handles.cur.ax,...
                            plotData.xData, plotData.yData); 
                    case 'plot'
                        this.figs.(gName).handles.cur.plot = ...
                            plot(this.figs.(gName).handles.cur.ax,...
                            plotData.xData, plotData.yData);
                end
                
                this.setLimsToPlot(gName); %apply limits to the plot
                this.setStringsToPlot(gName); %allocate title and labels handles
                drawnow();
                this.figs.(gName).update = false;
            else
                this.quickPlot(gName, plotData.xData, plotData.yData, [], [])
            end
        end

        function quickPlot(this, gName, xData, yData, cData, vars) 
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
                        'YPositiveDelta', vars,...
                        'YNegativeDelta', vars);
                case 'imagesc'
                    set(this.figs.(gName).handles.cur.plot,...
                        'XData', xData,...
                        'YData', yData,...
                        'CData', cData);
                    set(this.figs.(gName).handles.cur.ax, 'CLim', vars)
                    this.setLimsToPlot(gName)

            end
            set(this.figs.(gName).handles.cur.title, 'String', this.figs.(gName).strings.title);
            pause(0.01);
            drawnow();
        end

    end
end