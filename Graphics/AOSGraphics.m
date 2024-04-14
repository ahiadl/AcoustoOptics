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
            figsNames = {'scanPlane'; 'scanPlaneLog'; 'navLine'; 'navLineLog'; 'navPlane'; 'navPlaneLog'};
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
            
            figs.scan3AxType       = 'Normal';
            figs.scan3Vec       = [];
            figs.scan3VecCntr   = [];
            figs.scan3VecZero   = [];
            figs.scan3VecNorm   = [];

            figs.depthAxType       = 'Normal';
            figs.depthVec       = [];
            figs.depthVecCntr   = [];
            figs.depthVecZero   = [];
            figs.depthVecNorm   = [];
            
            figs.scan1Label     = 'Y';
            figs.scan2Label     = 'X';
            figs.scan3Label     = 'Z';
            figs.depthLabel     = 'D';
            figs.mainPlaneLabel = 'YZ';

            figs.scanAxLimsToPlot = [];
            figs.dScan1           = [];
            figs.dScan2           = [];
            figs.dScan3           = [];
            figs.dDDepth          = [];
 
            % Current Position
            figs.scan1Pos  = 0;
            figs.scan2Pos  = 0;
            figs.scan3Pos  = 0;         
            
            figs.scan1Idx = 1;
            figs.scan2Idx = 1;
            figs.scan3Idx = 1;

            % Line Navigator Position          
            figs.varsNavLine  = AOSGraphics.getLineNavVarsStruct;
            figs.varsNavPlane = AOSGraphics.getPlaneNavVarsStruct;
            
            figs.scan1LimsToPlot = [];
            figs.scan2LimsToPlot = [];
            figs.scan3LimsToPlot = [];
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
            figs.scan3AxType        = 'Normal';
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
            this.setStrings('scanPlane', "Scan Plane %s: (%s, %s) = (%.2f, %.2f) [mm] - Ch: %d", "Depth (US) %s [mm]", "Scan %s[mm]", []);
            
            % scanPlaneLog
            this.setType('scanPlaneLog', 'imagesc');
            this.setStrings('scanPlaneLog', "Scan Plane %s (Log): (%s, %s) = (%.2f, %.2f) [mm] - Ch: %d", "Depth (US) %s [mm] ", "Scan %s[mm]", []);
            
            % navLine
            this.setType('navLine', 'stem');
            this.setStrings('navLine', "Line Navigator %s: (%s, %s, %s) = (%.2f, %.2f, %.2f) - Ch: %d", "%s [mm]", "\\phi[v]", []);
            
            % navLineLog
            this.setType('navLineLog', 'plot');
            this.setStrings('navLineLog', "Line Navigator %s (Log): (%s, %s, %s) = (%.2f, %.2f, %.2f) - Ch: %d", "%s [mm]",  "\\phi[Log]", []);
            
             % navPlane
            this.setType('navPlane', 'imagesc');
            this.setStrings('navPlane', "Plane Navigator %s: (%s, %s) = (%.2f, %.2f) [mm] - Ch: %d", "%s[mm]", "%s[mm]", []);
            
            % navPlaneLog
            this.setType('navPlaneLog', 'imagesc');
            this.setStrings('navPlaneLog', "Plane Navigator %s (Log): (%s, %s) = (%.2f, %.2f) [mm] - Ch: %d", "%s[mm]", "%s[mm]", []);
        end

        function setGraphicsScanVars(this)
            this.figs.channelsInRes         = this.uVars.channelsInRes;
            this.figs.singleChannelAnalysis = this.uVars.singleChannelAnalysis;
            this.setSepChIdx(this.uVars.sepChIdx); 
            this.figs.normDispPlaneColors = this.uVars.normDispPlaneColors;

            this.figs.scan1Label   = this.uVars.scan1Label;
            this.figs.scan2Label   = this.uVars.scan2Label;
            this.figs.scan3Label   = this.uVars.scan3Label;
            this.figs.thirdLabel   = this.uVars.thirdLabel;
            this.figs.depthLabel   = 'D';
%             this.figs.depthLabel   = this.uVars.depthLabel;
            this.figs.mainPlaneLabel = this.uVars.mainPlane;
            
            this.figs.scan1AxType   = this.uVars.scan1AxType;
            this.figs.scan2AxType   = this.uVars.scan2AxType;
            this.figs.scan3AxType   = this.uVars.scan3AxType;
            this.figs.depthAxType   = this.uVars.depthAxType;

            this.figs.scan1VecNorm = this.uVars.scan1VecNorm;
            this.figs.scan1VecCntr = this.figs.scan1VecNorm - mean(this.figs.scan1VecNorm);
            this.figs.scan1VecZero = abs(this.figs.scan1VecNorm - this.figs.scan1VecNorm(1));
            this.setAxisType(this.figs.scan1Label);
            
            this.figs.scan2VecNorm = this.uVars.scan2VecNorm;
            this.figs.scan2VecCntr = this.figs.scan2VecNorm - mean(this.figs.scan2VecNorm);
            this.figs.scan2VecZero = abs(this.figs.scan2VecNorm - this.figs.scan2VecNorm(1));
            this.setAxisType(this.figs.scan2Label);
            
            this.figs.scan3VecNorm = this.uVars.scan3VecNorm;
            this.figs.scan3VecCntr = this.figs.scan3VecNorm - mean(this.figs.scan3VecNorm);
            this.figs.scan3VecZero = abs(this.figs.scan3VecNorm - this.figs.scan3VecNorm(1));
            this.setAxisType(this.figs.scan3Label);

            this.figs.depthVecNorm = this.uVars.depthVecNorm*1e3;
            this.figs.depthVecCntr = this.figs.depthVecNorm - mean(this.figs.depthVecNorm);
            this.figs.depthVecZero = abs(this.figs.depthVecNorm - this.figs.depthVecNorm(1));
            this.setAxisType(this.figs.depthLabel);
            
            this.figs.thirdPos  = this.uVars.thirdPos;            
            
            this.figs.scan1Idx = 1;
            this.figs.scan2Idx = 1;
            this.figs.scan3Idx = 1;
            this.figs.scan1Pos = 0;
            this.figs.scan2Pos = 0;
            this.figs.scan3Pos = 0;

            this.figs.varsNavLine.labels = this.uVars.varsNavLine.labels;
            this.figs.varsNavLine.idx    = this.uVars.varsNavLine.idx; 

            this.figs.varsNavPlane.labels = this.uVars.varsNavPlane.labels;
            this.figs.varsNavPlane.idx    = this.uVars.varsNavPlane.idx; 

            this.figs.validStruct = this.uVars.validStruct;

            this.updateGraphicsConstruction();
            this.resetPlots();         
        
            this.figs.fonts = this.uVars.fonts;
        end

        function updateCurPosAndIdx(this, pns)
            this.figs.scan1Idx  = pns.curIdx(1);
            this.figs.scan2Idx  = pns.curIdx(2);
            this.figs.scan3Idx  = pns.curIdx(3);
            this.figs.scan1Pos  = pns.curPos(1);
            this.figs.scan2Pos  = pns.curPos(2);
            this.figs.scan3Pos  = pns.curPos(3);
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
                    this.figs.dScan2 = abs(this.figs.scan2Vec(1) - this.figs.scan2Vec(2));
                end
                this.figs.scan2LimsToPlot = [min(this.figs.scan2Vec) - this.figs.dScan2,...
                                             max(this.figs.scan2Vec) + this.figs.dScan2];
            elseif strcmp(ax, this.figs.scan3Label)
                switch this.figs.scan3AxType
                    case 'Normal'
                        this.figs.scan3Vec = this.figs.scan3VecNorm;
                    case 'Center'
                        this.figs.scan3Vec = this.figs.scan3VecCntr;
                    case 'Zero'
                        this.figs.scan3Vec = this.figs.scan3VecZero;
                    case 'Index'
                        this.figs.scan3Vec = 1:length(this.figs.scan3VecNorm);
                end
                if length(this.figs.scan3Vec(1)) == 1
                    this.figs.dScan3 = 1;
                else
                    this.figs.dScan3 = abs(this.figs.scan3Vec(1) - this.figs.scan3Vec(2));
                end
                this.figs.scan3LimsToPlot = [min(this.figs.scan3Vec) - this.figs.dScan3,...
                                             max(this.figs.scan3Vec) + this.figs.dScan3];
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
            this.figs.varsNavLine.labels = vars.labels;
            this.figs.varsNavLine.idx    = vars.idx;
        end

        function setNavPlaneVars(this, vars)
            this.figs.varsNavPlane.labels = vars.normAx;
            this.figs.varsNavPlane.idx    = vars.axIdx;
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
            rstPhi    = zeros(length(this.figs.scan1Vec), length(this.figs.scan2Vec), length(this.figs.scan3Vec), length(this.figs.depthVec),  this.uVars.channelsInRes);
            this.data.phi       = rstPhi;
            this.data.phiLog    = rstPhi;
            this.data.phiNorm   = rstPhi;

            this.data.glbMax  = ones(1, this.uVars.channelsInRes);
            this.data.glbMin  = zeros(1, this.uVars.channelsInRes);
            this.data.glbSpan = ones(1, this.uVars.channelsInRes);

            this.data.glbMinLog = -Inf*ones(1, this.uVars.channelsInRes);
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
            elseif strcmp(label, this.figs.scan3Label)
                dim = 3;
            elseif strcmp(label, this.figs.depthLabel)
                dim = 4;
            end
        end
        
        function lineNavData = extractLineNavAxFromLabels(this)
            
            labels = this.figs.varsNavLine.labels;
            idx    = this.figs.varsNavLine.idx;

            if strcmp(labels(1), this.figs.scan1Label)
                lineNavData.lineAxLabel = this.figs.scan1Label;
                lineNavData.lineAxVec = this.figs.scan1Vec;
            elseif strcmp(labels(1), this.figs.scan2Label)
                lineNavData.lineAxLabel = this.figs.scan2Label;
                lineNavData.lineAxVec = this.figs.scan2Vec;
            elseif strcmp(labels(1), this.figs.scan3Label)
                lineNavData.lineAxLabel = this.figs.scan3Label;
                lineNavData.lineAxVec = this.figs.scan3Vec;
            elseif strcmp(labels(1), this.figs.depthLabel)
                lineNavData.lineAxLabel = this.figs.depthLabel;
                lineNavData.lineAxVec = this.figs.depthVec;
            end

            for i=2:4
                axLabelField = sprintf("ax%dLabel", i);
                axVecField   = sprintf("ax%dVec", i);
                axValField   = sprintf("ax%dVal", i);
                axIdxField   = sprintf("ax%dIdx", i);
                
                if strcmp(labels(i), this.figs.scan1Label)
                    lineNavData.(axLabelField) = this.figs.scan1Label;
                    lineNavData.(axVecField) = this.figs.scan1Vec;
                    lineNavData.(axValField) = this.figs.scan1Vec(idx(i-1));
                elseif strcmp(labels(i), this.figs.scan2Label)
                    lineNavData.(axLabelField) = this.figs.scan2Label;
                    lineNavData.(axVecField) = this.figs.scan2Vec;
                    lineNavData.(axValField) = this.figs.scan2Vec(idx(i-1));
                elseif strcmp(labels(i), this.figs.scan3Label)
                    lineNavData.(axLabelField) = this.figs.scan3Label;
                    lineNavData.(axVecField) = this.figs.scan3Vec;
                    lineNavData.(axValField) = this.figs.scan3Vec(idx(i-1));
                elseif strcmp(labels(i), this.figs.depthLabel)
                    lineNavData.(axLabelField) = this.figs.depthLabel;
                    lineNavData.(axVecField) = this.figs.depthVec;
                    lineNavData.(axValField) = this.figs.depthVec(idx(i-1));
                end
                lineNavData.(axIdxField) = idx(i-1);
            end
        end

        function planeData = extractPlaneNavAxesFromLabel(this)
            labels = this.figs.varsNavPlane.labels;
            idx    = this.figs.varsNavPlane.idx;

            for i=1:2
                axLabelField = sprintf("line%dAxLabel", i);
                axVecField   = sprintf("line%dAxVec", i);
                if strcmp(labels(i), this.figs.scan1Label)
                    planeData.(axLabelField) = this.figs.scan1Label;
                    planeData.(axVecField) = this.figs.scan1Vec;
                elseif strcmp(labels(i), this.figs.scan2Label)
                    planeData.(axLabelField) = this.figs.scan2Label;
                    planeData.(axVecField) = this.figs.scan2Vec;
                elseif strcmp(labels(i), this.figs.scan3Label)
                    planeData.(axLabelField) = this.figs.scan3Label;
                    planeData.(axVecField) = this.figs.scan3Vec;
                elseif strcmp(labels(i), this.figs.depthLabel)
                    planeData.(axLabelField) = this.figs.depthLabel;
                    planeData.(axVecField) = this.figs.depthVec;
                end
            end

            for i=3:4
                axLabelField = sprintf("ax%dLabel", i);
                axVecField   = sprintf("ax%dVec", i);
                axValField   = sprintf("ax%dVal", i);
                axIdxField   = sprintf("ax%dIdx", i);
                
                if strcmp(labels(i), this.figs.scan1Label)
                    planeData.(axLabelField) = this.figs.scan1Label;
                    planeData.(axVecField) = this.figs.scan1Vec;
                    planeData.(axValField) = this.figs.scan1Vec(idx(i-2));
                elseif strcmp(labels(i), this.figs.scan2Label)
                    planeData.(axLabelField) = this.figs.scan2Label;
                    planeData.(axVecField) = this.figs.scan2Vec;
                    planeData.(axValField) = this.figs.scan2Vec(idx(i-2));
                elseif strcmp(labels(i), this.figs.scan3Label)
                    planeData.(axLabelField) = this.figs.scan3Label;
                    planeData.(axVecField) = this.figs.scan3Vec;
                    planeData.(axValField) = this.figs.scan3Vec(idx(i-2));
                elseif strcmp(labels(i), this.figs.depthLabel)
                    planeData.(axLabelField) = this.figs.depthLabel;
                    planeData.(axVecField) = this.figs.depthVec;
                    planeData.(axValField) = this.figs.depthVec(idx(i-2));
                end
                planeData.(axIdxField) = idx(i-2);
            end

            planeData.planeLabel = sprintf("%s%s", planeData.line1AxLabel, planeData.line2AxLabel);
        end

        function data = extractScanPlaneFromData(this, gName)
            data.xData = this.figs.depthVec;
            data.yData = this.figs.scan1Vec;

            switch gName
                case 'scanPlane'
                    data.cData = permute(this.data.phiNorm(:, this.figs.scan2Idx, this.figs.scan3Idx, :, this.figs.sepChIdx), [1,4,2,3,5]);
                    if this.figs.normDispPlaneColors
                        data.clims = [min(data.cData(:)), max(data.cData(:))];
                    else
                        data.clims = [0, 1];
                    end
                case 'scanPlaneLog'
                    data.cData = permute(this.data.phiLog(:, this.figs.scan2Idx, this.figs.scan3Idx, :, this.figs.sepChIdx), [1,4,2,3,5]);
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
                                this.figs.scan3Label,...
                                this.figs.scan2Pos,...
                                this.figs.scan3Pos,...
                                this.figs.sepChIdx}};
            data.xLabelVar = this.figs.depthLabel;
            data.yLabelVar = this.figs.scan1Label;
        end

        function data = extractNavPlaneFromData(this, gName)
            dimOrder = zeros(1,5);
            for i=1:4
                dimOrder(i) = this.extractDimFromLabel(this.figs.varsNavLine.labels(i));
            end
            dimOrder(5) = 5;
            
            planeData = this.extractPlaneNavAxesFromLabel();
            data.xData = planeData.line2AxVec;
            data.yData = planeData.line1AxVec;
            
            switch gName
                case 'navPlane'
                    phi = permute(this.data.phiNorm, dimOrder);
                    phi = squeeze(phi(:, :, planeData.ax3Idx,  planeData.ax4Idx, this.figs.sepChIdx));
                    data.cData = phi;
                    if this.figs.normDispPlaneColors
                        data.clims = [min(phi(:)), max(phi(:))];
                    else
                        data.clims = [0, 1];
                    end
                case 'navPlaneLog'
                    phiLog = permute(this.data.phiLog, dimOrder);
                    phiLog = squeeze(phiLog(:, :, planeData.ax3Idx,  planeData.ax4Idx, this.figs.sepChIdx));
                    data.cData = phiLog;
                    if this.figs.normDispPlaneColors
                        data.clims = [min(phiLog(:)), max(phiLog(:))];
                    else
                        data.clims = [this.data.glbMinLog(this.figs.sepChIdx), 0];
                    end
            end
            if data.clims(1) == data.clims(2)
                data.clims(2)= data.clims(1)+1;
            end
            data.titleVars = {{ planeData.planeLabel,...
                                planeData.ax3Label,...
                                planeData.ax4Label,...
                                planeData.ax3Val,...
                                planeData.ax4Val,...
                                this.figs.sepChIdx}};

            data.xLabelVar = planeData.line2AxLabel;
            data.yLabelVar = planeData.line1AxLabel;
        end

        function data = extractNavLineFromData(this, gName)
            dimOrder = zeros(1,5);
            for i=1:4
                dimOrder(i) = this.extractDimFromLabel(this.figs.varsNavLine.labels(i));
            end
            dimOrder(5) = 5;

            lineData = this.extractLineNavAxFromLabels();
            data.xData = lineData.lineAxVec;

            switch gName
                case "navLine"
                    phi = permute(this.data.phiNorm, dimOrder);
                case "navLineLog"
                    phi = permute(this.data.phiLog, dimOrder);
            end
            data.yData = squeeze(phi(:, lineData.ax2Idx, lineData.ax3Idx, lineData.ax4Idx, this.figs.sepChIdx));
                        
            data.titleVars  =  {{lineData.lineAxLabel,...
                                 lineData.ax2Label,...
                                 lineData.ax3Label,...
                                 lineData.ax4Label,...
                                 lineData.ax2Val,...
                                 lineData.ax3Val,...
                                 lineData.ax4Val,...
                                 this.figs.sepChIdx}};

            data.xLabelVar = lineData.lineAxLabel;
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