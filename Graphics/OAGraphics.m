classdef OAGraphics < Graphics
    %AOGRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data
        dataScan
    end
    
    methods (Static)
        function gNames = getGraphicsNames()
            gNames = {'sinogram'; 'recon'; 'navSinogram'; 'navRecon'};
        end
        
        function figs  = createGraphicsVars()
            figs.firstAxLabel  = 'X';
            figs.secondAxLabel = 'Y';
            figs.scanAxLabel   = 'Z';
            
            figs.firstAxType  = 'Normal';
            figs.secondAxType = 'Normal';
            figs.scanAxType   = 'Normal';
            figs.timeAxType   = 'Normal';
            
            figs.axesVec.firstAxisNormal  = [];
            figs.axesVec.firstAxisZero    = [];
            figs.axesVec.firstAxisCenter  = [];
            figs.axesVec.firstAxisIdx     = [];
            
            figs.axesVec.secondAxisNormal = [];
            figs.axesVec.secondAxisZero   = [];
            figs.axesVec.secondAxisCenter = [];
            figs.axesVec.secondAxisIdx    = [];
            
            figs.axesVec.scanAxisNormal   = [];
            figs.axesVec.scanAxisZero     = [];
            figs.axesVec.scanAxisCenter   = [];
            figs.axesVec.scanAxisIdx      = [];
            
            figs.axesVec.timeAxNormal     = [];
            figs.axesVec.timeAxSamp       = [];

            figs.chAxVec       = [];
            figs.firstAxVec    = [];
            figs.secondAxVec   = [];
            figs.scanAxVec     = [];
            figs.timeAxVec     = [];

            figs.curPosIdx = [];
            figs.curPos    = [];
            figs.navIdx    = [];
            figs.navPos    = [];
            
            figs.numOfChannels = 256;
            figs.numOfSamples  = 2030;
            figs.reconSize     = 400;
            figs.imageWidth    = 10e-3;
            
            figs.extClims      = [];
            figs.useExtClims   = false;
            figs.reopenFigures = false;
            figs.intExt        = 'int';
            figsNames          = OAGraphics.getGraphicsNames();

            figs.fonts.type       = [];
            figs.fonts.titleSize  = 14;
            figs.fonts.labelsSize = 14;
            figs.fonts.axisSize   = 14;
            
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
                figs.(figsNames{i}).strings.updateLegend = true;

                figs.(figsNames{i}).dims.zDim  = [];
                figs.(figsNames{i}).dims.chDim   = [];
                figs.(figsNames{i}).dims.dataDim = [];

                figs.(figsNames{i}).markers.plotMark = false;

                figs.(figsNames{i}).lims.xlims  = [];
                figs.(figsNames{i}).lims.ylims  = [];
                figs.(figsNames{i}).lims.clims = [];
            end
        end 
        
        function uVars = createUserVars()
            %Variables that comes from the user
            figsNames = OAGraphics.getGraphicsNames();
                       
            uVars.firstAxType  = 'Normal';
            uVars.secondAxType = 'Normal';
            uVars.scanAxType   = 'Normal';
            uVars.timeAxType   = 'Normal';
            
            uVars.extClims      = [];
            uVars.useExtClims   = false;
            uVars.reopenFigures = false;
            
            uVars.intExt = [];
            for i=1:length(figsNames)
                uVars.validStruct.(figsNames{i}) = false;
            end
            
            for i=1:length(figsNames)
                uVars.extH.(figsNames{i}) =  Graphics.createHandlesStruct();
            end
        end
        
        function vars  = createOwnerVars()
            %Variables that comes from the AO Object
            vars = OAGraphics.createUserVars();
            vars.grid = [];
            
            vars.firstAxLabel  = 'X';
            vars.secondAxLabel = 'Y';
            vars.scanAxLabel   = 'Z';

            vars.numOfChannels = 256;
            vars.numOfSamples  = 2030;
            vars.reconSize     = 400;
            vars.imageWidth    = 10e-3;
            vars.tVec          = [];
            vars.reconAxis     = [];
        end 
    end
    
    methods
        function this = OAGraphics()
            this@Graphics()
            
            this.figsNames = this.getGraphicsNames();
            this.numOfFigs = length(this.figsNames);
            
            this.figs     = OAGraphics.createGraphicsVars();
            
            this.uVars    = OAGraphics.createOwnerVars();
            
            this.setGraphicsStaticVars();
        end
        
        % Set Functions 
        function setGraphicsDynamicVars(this)           
            %Generic functionality (written to uVars to sync with
            %construction function).
            this.figs.validStruct     = this.uVars.validStruct;
            this.figs.extH            = this.uVars.extH;
            this.figs.reopenFigures   = this.uVars.reopenFigures;
            
            %Measurement Parameters
            this.figs.numOfChannels = this.uVars.numOfChannels;
            this.figs.numOfSamples  = this.uVars.numOfSamples;
            this.figs.reconSize     = this.uVars.reconSize;
            this.figs.imageWidth    = this.uVars.imageWidth;

            %Axes Parameters
            this.figs.firstAxLabel  = this.uVars.firstAxLabel;
            this.figs.secondAxLabel = this.uVars.secondAxLabel;
            this.figs.scanAxLabel   = this.uVars.scanAxLabel;
            
            this.figs.firstAxType  = this.uVars.firstAxType;
            this.figs.secondAxType = this.uVars.secondAxType;
            this.figs.scanAxType   = this.uVars.scanAxType;
            this.figs.timeAxType   = this.uVars.timeAxType;
            
            this.figs.axesVec.firstAxisNormal  = this.uVars.grid.firstNormal;
            this.figs.axesVec.firstAxisZero    = this.uVars.grid.firstZero;
            this.figs.axesVec.firstAxisCenter  = this.uVars.grid.firstCntr;
            this.figs.axesVec.firstAxisIdx     = this.uVars.grid.firstIdx;
            
            this.figs.axesVec.secondAxisNormal = this.uVars.grid.secondNormal;
            this.figs.axesVec.secondAxisZero   = this.uVars.grid.secondZero;
            this.figs.axesVec.secondAxisCenter = this.uVars.grid.secondCntr;
            this.figs.axesVec.secondAxisIdx    = this.uVars.grid.secondIdx;
            
            this.figs.axesVec.scanAxisNormal   = this.uVars.grid.scanNormal;
            this.figs.axesVec.scanAxisZero     = this.uVars.grid.scanZero;
            this.figs.axesVec.scanAxisCenter   = this.uVars.grid.scanCntr;
            this.figs.axesVec.scanAxisIdx      = this.uVars.grid.scanIdx;
            
            this.figs.axesVec.timeAxNormal     = this.uVars.tVec;
            this.figs.axesVec.timeAxSamp       = 1:this.figs.numOfSamples;
            
%             this.figs.reconAxis     = this.uVars.reconAxis;
%             this.figs.reconIndAxis  = -this.figs.reconSize/2 : 1 : this.figs.reconSize/2;
            this.figs.chAxVec       = 1:this.figs.numOfChannels;
            this.figs.firstAxVec    = this.getAxByType('first');
            this.figs.secondAxVec   = this.getAxByType('second');
            this.figs.scanAxVec     = this.getAxByType('scan');
            this.figs.timeAxVec     = this.getAxByType('time');
            
            this.figs.useExtClims = this.uVars.useExtClims;
            
            this.figs.curPosIdx = 1;
            this.figs.curPos    = this.figs.scanAxVec(1);
            this.updateNavigatorIdxAndPos(1);
            
            this.updateGraphicsConstruction()
            this.resetGraphics();
        end
        
        function setGraphicsStaticVars(this)
%             this.graphics.setType(type);
%             this.graphics.setStrings(title, xlabel, ylabel, legend);
%             this.graphics.setDims(gName, zIdxDim, chDim, dataDim)
%             this.graphics.setMarkersEnable(enable);          
           
            this.setType('sinogram', 'imagesc');
            this.setStrings('sinogram',"Sinogram %s = %d",  "Samples[#]", "Channels[#]",[]);
            
            this.setType('recon', 'imagesc');
            this.setStrings('recon', "Reconstruction %s = %d", "%s[mm]", "%s[mm]", []);
            
            this.setType('navSinogram', 'imagesc');
            this.setStrings('navSinogram',"Navigator Sinogram %s = %d",  "Samples[#]", "Channels[#]",[]);
            
            this.setType('navRecon', 'imagesc');
            this.setStrings('navRecon', "Navigator Reconstruction %s = %d", "%s[mm]", "%s[mm]", []);
        end

        function updateCurPos(this, pos, idx)
            this.figs.curPos    = pos;
            this.figs.curPosIdx = idx;
        end
        
        function setData(this, data)
            this.data.sigMat = data.sigMat;
            this.data.recon  = data.recon;
        end
        
        function setDataScan(this, data)
            this.dataScan.sigMat = data.sigMat;
            this.dataScan.recon  = data.recon;
            
            this.data.sigMat = data.sigMat(:,:, this.figs.curPosIdx);
            this.data.recon  = data.recon(:,:, this.figs.curPosIdx);
        end
        
        function updateNavigatorIdxAndPos(this, idx)
            if idx > this.figs.curPosIdx
                return;
            end
            this.figs.navIdx = idx;
            this.figs.navPos = this.figs.scanAxVec(idx);
        end
        
        function updateNavigator(this, idx)
            this.updateNavigatorIdxAndPos(idx);
            this.dispNavSinogram();
            this.dispNavRecon();
        end
        
        function resetGraphics(this)
            this.data.recon  = zeros(this.uVars.grid.firstIdxLen, this.uVars.grid.secondIdxLen);
            this.data.sigMat = zeros(this.figs.numOfChannels, this.figs.numOfSamples);
            
            this.dataScan.recon  = zeros(this.uVars.grid.firstIdxLen, this.uVars.grid.secondIdxLen, this.uVars.grid.scanIdxLen);
            this.dataScan.sigMat = zeros(this.figs.numOfChannels, this.figs.numOfSamples, this.uVars.grid.scanIdxLen);
            
            this.dispSinogram();
            this.dispRecon();
            this.dispNavSinogram();
            this.dispNavRecon();
        end
        % Get functions
        function axVec = getAxByType(this, ax)
            switch ax
                case 'first'
                    switch this.figs.firstAxType
                        case 'Normal'
                            axVec = this.figs.axesVec.firstAxisNormal;
                        case 'Center'
                            axVec = this.figs.axesVec.firstAxisCenter;
                        case 'Zero'
                            axVec = this.figs.axesVec.firstAxisZero;
                        case 'Index'
                            axVec = this.figs.axesVec.firstAxisIdx;
                    end
                case 'second'
                    switch this.figs.secondAxType
                        case 'Normal'
                            axVec = this.figs.axesVec.secondAxisNormal;
                        case 'Center'
                            axVec = this.figs.axesVec.secondAxisCenter;
                        case 'Zero'
                            axVec = this.figs.axesVec.secondAxisZero;
                        case 'Index'
                            axVec = this.figs.axesVec.secondAxisIdx;
                    end
                case 'scan'
                    switch this.figs.scanAxType
                        case 'Normal'
                            axVec = this.figs.axesVec.scanAxisNormal;
                        case 'Center'
                            axVec = this.figs.axesVec.scanAxisCenter;
                        case 'Zero'
                            axVec = this.figs.axesVec.scanAxisZero;
                        case 'Index'
                            axVec = this.figs.axesVec.scanAxisIdx;
                    end
                case 'time'
                    switch this.figs.scanAxType
                        case 'Normal'
                            axVec = this.figs.axesVec.timeAxNormal;
                        case 'Sample'
                            axVec = this.figs.axesVec.timeAxSamp;
                    end
            end
            
        end
        
        function signal = extractData(this, figName)
            switch figName
                case 'sinogram'
                    signal.xData = this.figs.chAxVec;  
                    signal.yData = this.figs.timeAxVec;
                    signal.cData = this.data.sigMat;
                    
                    if this.figs.useExtClims
                        signal.cLims = this.figs.extClims.sinogram;
                    else
                        signal.clims = [min(min(signal.cData)), max(max(signal.cData))];
                    end
                case 'recon'
                    signal.xData = this.figs.firstAxVec;
                    signal.yData = this.figs.secondAxVec;
                    signal.cData = this.data.recon;
                    
                    if this.figs.useExtClims
                        signal.cLims = this.figs.extClims.recon;
                    else
                        signal.clims = [min(min(signal.cData)), max(max(signal.cData))];
                    end
                case 'navSinogram'
                    signal.xData = this.figs.chAxVec;
                    signal.yData = this.figs.timeAxVec;
                    
                    signal.cData = this.dataScan.sigMat(:,:, this.figs.navIdx);
                    
                    if this.figs.useExtClims
                        signal.cLims = this.figs.extClims.sinogram;
                    else
                        signal.clims = [min(min(signal.cData)), max(max(signal.cData))];
                    end
                case 'navRecon'
                    signal.xData = this.figs.firstAxVec;
                    signal.yData = this.figs.secondAxVec;
                    signal.cData = this.dataScan.recon(:,:,this.figs.navIdx);
                    
                    if this.figs.useExtClims
                        signal.cLims = this.figs.extClims.recon;
                    else
                        signal.clims = [min(min(signal.cData)), max(max(signal.cData))];
                    end   
            end
            if signal.clims(1) == signal.clims(2)
                signal.clims(2) = signal.clims(1)+ 0.1;
            end
        end

        % Disply Functions
        function quickPlot(this, gName, xData, yData, cData) 
            switch this.figs.(gName).type
                case 'plot'    
                    set(this.figs.(gName).handles.cur.plot,...
                        'XData', xData,...
                        'YData', yData);
                case 'imagesc'
                    set(this.figs.(gName).handles.cur.plot,...
                        'XData', xData,...
                        'YData', yData,...
                        'CData', cData);
                    this.setLimsToPlot(gName)
            end 
            set(this.figs.(gName).handles.cur.title, 'String', this.figs.(gName).strings.title);
            pause(0.001);
            drawnow();
         end 
        
        function dispSinogram(this)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.sinogram.handles.cur.ax)
               return
            end
            % "Current Main Plane %s: (R, %s) = (%d, %.2f)"
            this.setTitleVariables('sinogram',...
                                   { {this.figs.scanAxLabel, this.figs.curPos}});
            
            plotData = this.extractData('sinogram');
                                                  
            this.figs.sinogram.lims.clims = plotData.clims;

            if ~isgraphics(this.figs.sinogram.handles.cur.plot) ||...
                this.figs.sinogram.update
                
                cla(this.figs.sinogram.handles.cur.ax)
                this.setAxesVar('sinogram', this.figs.firstAxLabel, this.figs.secondAxLabel);
                
                this.figs.sinogram.handles.cur.plot = ...
                    imagesc(this.figs.sinogram.handles.cur.ax,...
                    'XData', plotData.xData, 'Ydata', plotData.yData,...
                    'CData', plotData.cData, plotData.clims);
                
                %title
                axis(this.figs.sinogram.handles.cur.ax, 'tight')
                colorbar(this.figs.sinogram.handles.cur.ax);
                this.setLimsToPlot('sinogram');
                this.setStringsToPlot('sinogram');
                drawnow();
                this.figs.sinogram.update = false;
            else
                quickPlot(this, 'sinogram',  plotData.xData, plotData.yData, plotData.cData)
            end
        end
        
        function dispRecon(this)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.recon.handles.cur.ax)
               return
            end
            % "Current Main Plane %s: (R, %s) = (%d, %.2f)"
            this.setTitleVariables('recon',...
                                   {{this.figs.scanAxLabel, this.figs.curPos}});
            
            plotData = this.extractData('recon');
                                                  
            this.figs.recon.lims.clims = plotData.clims;

            if ~isgraphics(this.figs.recon.handles.cur.plot) ||...
                this.figs.recon.update
                
                cla(this.figs.recon.handles.cur.ax)
                this.setAxesVar('recon', this.figs.firstAxLabel, this.figs.secondAxLabel);
                
                this.figs.recon.handles.cur.plot = ...
                    imagesc(this.figs.recon.handles.cur.ax,...
                    'XData', plotData.xData, 'Ydata', plotData.yData,...
                    'CData', plotData.cData, plotData.clims);
                
                %title
                axis(this.figs.recon.handles.cur.ax, 'equal', 'tight')
                colorbar(this.figs.recon.handles.cur.ax);
                this.setLimsToPlot('recon');
                this.setStringsToPlot('recon');
                drawnow();
                this.figs.recon.update = false;
            else
                quickPlot(this, 'recon', plotData.xData, plotData.yData, plotData.cData)
            end
        end
        
        function dispNavSinogram(this)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.navSinogram.handles.cur.ax)
               return
            end
            % "Current Main Plane %s: (R, %s) = (%d, %.2f)"
            this.setTitleVariables('navSinogram',...
                                   { {this.figs.scanAxLabel, this.figs.navPos}});
            
            plotData = this.extractData('navSinogram');
                                                  
            this.figs.navSinogram.lims.clims = plotData.clims;

            if ~isgraphics(this.figs.navSinogram.handles.cur.plot) ||...
                this.figs.navSinogram.update
                
                cla(this.figs.navSinogram.handles.cur.ax)
                this.setAxesVar('navSinogram', this.figs.firstAxLabel, this.figs.secondAxLabel);
                
                this.figs.navSinogram.handles.cur.plot = ...
                    imagesc(this.figs.navSinogram.handles.cur.ax,...
                    'XData', plotData.xData, 'Ydata', plotData.yData,...
                    'CData', plotData.cData, plotData.clims);
                
                %title
                axis(this.figs.navSinogram.handles.cur.ax, 'tight')
                colorbar(this.figs.navSinogram.handles.cur.ax);
                this.setLimsToPlot('navSinogram');
                this.setStringsToPlot('navSinogram');
                drawnow();
                this.figs.navSinogram.update = false;
            else
                quickPlot(this, 'navSinogram',  plotData.xData, plotData.yData, plotData.cData)
            end
        end
        
        function dispNavRecon(this)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.navRecon.handles.cur.ax)
               return
            end
            % "Current Main Plane %s: (R, %s) = (%d, %.2f)"
            this.setTitleVariables('navRecon',...
                                   {{this.figs.scanAxLabel, this.figs.navPos}});
            
            plotData = this.extractData('navRecon');
                                                  
            this.figs.navRecon.lims.clims = plotData.clims;

            if ~isgraphics(this.figs.navRecon.handles.cur.plot) ||...
                this.figs.navRecon.update
                
                cla(this.figs.navRecon.handles.cur.ax)
                this.setAxesVar('navRecon', this.figs.firstAxLabel, this.figs.secondAxLabel);
                
                this.figs.navRecon.handles.cur.plot = ...
                    imagesc(this.figs.navRecon.handles.cur.ax,...
                    'XData', plotData.xData, 'Ydata', plotData.yData,...
                    'CData', plotData.cData, plotData.clims);
                
                %title
                axis(this.figs.navRecon.handles.cur.ax, 'equal', 'tight')
                colorbar(this.figs.navRecon.handles.cur.ax);
                this.setLimsToPlot('navRecon');
                this.setStringsToPlot('navRecon');
                drawnow();
                this.figs.navRecon.update = false;
            else
                quickPlot(this, 'navRecon', plotData.xData, plotData.yData, plotData.cData)
            end
        end
        
    end
end
