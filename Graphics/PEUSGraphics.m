classdef PEUSGraphics < Graphics
    %AOGRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data
    end
    
    methods (Static)
        function gNames = getGraphicsNames()
            gNames = {'rawData'; 'recon'};
        end
        
        function figs  = createGraphicsVars()
            figs.intExt           = 'int';
            figs.reopenFigures    = false;
            figs.useExtClims = [];
            
            figs.numOfSamples  = [];
            figs.imageWidth    = [];
            
            figs.tVec          = [];
            
            figs.depthAxis     = [];
            figs.depthIndAxis  = [];
            figs.depthAxisCntr = [];
            
            figs.scanAxis      = [];
            figs.scanIndAxis   = [];
            figs.scanAxisCntr  = [];

            figs.depthAxLabel  = 'Z';
            figs.scanAxLabel   = 'X';
            figs.thirdAxLabel  = 'Y';
            
            figs.depthAxType   = [];
            figs.scanAxType    = [];
            
            figs.thirdAxCoor   = [];
            
            figs.depthAxisLen  = []; 
            figs.scanAxisLen   = [];
            
            figsNames = PEUSGraphics.getGraphicsNames();

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
            figsNames = PEUSGraphics.getGraphicsNames();
            
            uVars.intExt        = [];
            uVars.useExtClims   = false;
            uVars.reopenFigures = false;

            uVars.depthAxType   = [];
            uVars.scanAxType    = [];

            for i=1:length(figsNames)
                uVars.validStruct.(figsNames{i}) = false;
            end
            
            for i=1:length(figsNames)
                uVars.extH.(figsNames{i}) =  Graphics.createHandlesStruct();
            end
            
            uVars.fonts.type       = [];
            uVars.fonts.titleSize  = 14;
            uVars.fonts.labelsSize = 14;
            uVars.fonts.axisSize   = 14;
        end
        
        function vars  = createOwnerVars()
            vars = PEUSGraphics.createUserVars();
            
            vars.numOfSamples  = [];
            vars.imageWidth    = [];
            
            vars.tVec          = [];
            
            vars.depthAxis     = [];
            vars.depthAxisCntr = [];
            
            vars.scanAxis      = [];
            vars.scanAxisCntr  = [];
            
            vars.depthAxLabel = 'X';
            vars.scanAxLabel  = 'Y';
            vars.thirdAxLabel = 'Z';
            
            vars.thirdAxCoor = [];
        end 
    end
    
    methods
        function this = PEUSGraphics()
            this@Graphics()
            
            this.figsNames = this.getGraphicsNames();
            this.numOfFigs = length(this.figsNames);
            
            this.figs     = PEUSGraphics.createGraphicsVars();
            
            this.uVars    = PEUSGraphics.createOwnerVars();
            
            this.setGraphicsStaticVars();
        end
        
        % Set Functions     
        function setGraphicsScanVars(this)
            % Generic functionality
            this.figs.reopenFigures   = this.uVars.reopenFigures;
            this.figs.useExtClims     = this.uVars.useExtClims;
            
            this.figs.validStruct     = this.uVars.validStruct;
            this.figs.extH            = this.uVars.extH;     
            
            % Measurement Parameters
            this.figs.numOfSamples  = this.uVars.numOfSamples;
            this.figs.imageWidth    = this.uVars.imageWidth;
                     
            % Axes Parameters
            this.figs.tVec          = this.uVars.tVec*1e6;
            
            this.figs.depthAxis     = this.uVars.depthAxis *1e3;
            this.figs.depthIndAxis  = 1 : length(this.figs.depthAxis);
            this.figs.depthAxisCntr = this.uVars.depthAxisCntr*1e3;
            
            this.figs.scanAxis      = this.uVars.scanAxis;
            this.figs.scanIndAxis   = 1 : length(this.figs.scanAxis);
            this.figs.scanAxisCntr  = this.uVars.scanAxisCntr;
            
            this.figs.depthAxType   = this.uVars.depthAxType;
            this.figs.scanAxType    = this.uVars.scanAxType;
            
            this.figs.depthAxLabel  = this.uVars.depthAxLabel;
            this.figs.scanAxLabel   = this.uVars.scanAxLabel;
            this.figs.thirdAxLabel  = this.uVars.thirdAxLabel;
            
            this.figs.depthAxisLen  = length(this.figs.depthAxis);
            this.figs.scanAxisLen   = length(this.figs.scanAxis);
            
            % Title Parameters
            this.figs.thirdAxCoor  = this.uVars.thirdAxCoor;
            
            % Fonts
            this.figs.fonts = this.uVars.fonts;
            
            % Open/Close Figures
            this.updateGraphicsConstruction()
            this.resetPlots();
        end
        
        function setGraphicsStaticVars(this)
%             this.graphics.setType(type);
%             this.graphics.setStrings(title, xlabel, ylabel, legend);
%             this.graphics.setDims(gName, zIdxDim, chDim, dataDim)
%             this.graphics.setMarkersEnable(enable);          
           
            % phi [numOfPos]
            this.setType('rawData', 'imagesc');
            this.setStrings('rawData',"Raw Data %s = %d",  "t[\\mus]", "Scan Axis[mm]",[]);
            
            % phiRaw [numOfPos]
            this.setType('recon', 'imagesc');
            this.setStrings('recon', "Reconstruction %s = %d", "%s[mm] (Depth)", "%s[mm] (Scan)", []);
        end
        
        function resetPlots(this)
            this.initDataArray();
            this.plotRawData();
            this.plotRecon();
        end
        
        function resetGraphics(this)
            this.initDataArray()
        end
        
        function initDataArray(this)
            this.data.rawData = zeros(this.figs.depthAxisLen, this.figs.scanAxisLen);
            this.data.recon   = zeros(this.figs.depthAxisLen, this.figs.scanAxisLen);
        end
        
        function setData(this, data)
            this.data = data;
        end
         
        function signal = extractData(this, figName)
            switch figName
                case 'rawData'
                    switch this.figs.scanAxType
                        case 'Index'
                            signal.yData = this.figs.scanIndAxis;
                        case 'Normal'
                            signal.yData = this.figs.scanAxis;
                        case 'Center'
                            signal.yData = this.figs.scanAxisCntr;  
                    end
                    signal.xData = this.figs.tVec;
                    
                    signal.cData = this.data.rawData; %TODO: check if need to be transposed
                    if this.figs.useExtClims
                        signal.cLims = this.figs.extClims.rawData;
                    else
                        signal.clims = [min(min(signal.cData)), max(max(signal.cData))];
                        if signal.clims(1) == signal.clims(2)
                            signal.clims(2) = signal.clims(1) +1;
                        end
                    end
                case 'recon'
                     switch this.figs.scanAxType
                        case 'Index'
                            signal.yData = this.figs.scanIndAxis;
                        case 'Normal'
                            signal.yData = this.figs.scanAxis;
                        case 'Center'
                            signal.yData = this.figs.scanAxisCntr;  
                     end
                    
                     switch this.figs.depthAxType
                        case 'Index'
                            signal.xData = this.figs.depthIndAxis;
                        case 'Normal'
                            signal.xData = this.figs.depthAxis;
                        case 'Center'
                            signal.xData = this.figs.depthAxisCntr;  
                     end
                    
                    signal.cData = this.data.recon;
                    
                    if this.figs.useExtClims
                        signal.cLims = this.figs.extClims.recon;
                    else
                        signal.clims = [min(min(signal.cData)), max(max(signal.cData))];
                        if signal.clims(1) == signal.clims(2)
                            signal.clims(2) = signal.clims(1) +1;
                        end
                    end
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
        
        function plotRawData(this)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.rawData.handles.cur.ax)
               return
            end

            plotData = this.extractData('rawData');
                                                  
            this.figs.rawData.lims.clims = plotData.clims;

            if ~isgraphics(this.figs.rawData.handles.cur.plot) ||...
                this.figs.rawData.update
                
                this.setTitleVariables('rawData',...
                	{ {this.figs.thirdAxLabel, this.figs.thirdAxCoor}});
                this.setAxesVar('rawData', this.figs.depthAxLabel, this.figs.scanAxLabel);
                
                cla(this.figs.rawData.handles.cur.ax)
                
                this.figs.rawData.handles.cur.plot = ...
                    imagesc(this.figs.rawData.handles.cur.ax,...
                    'XData', plotData.xData, 'Ydata', plotData.yData,...
                    'CData', plotData.cData, plotData.clims);
                
                axis(this.figs.rawData.handles.cur.ax, 'tight')
                colorbar(this.figs.rawData.handles.cur.ax);
                this.setLimsToPlot('rawData');
                this.setStringsToPlot('rawData');
                drawnow();
                this.figs.rawData.update = false;
            else
                quickPlot(this, 'rawData',  plotData.xData, plotData.yData, plotData.cData)
            end
        end
        
        function plotRecon(this)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.recon.handles.cur.ax)
               return
            end
   
            plotData = this.extractData('recon');
                                                  
            this.figs.recon.lims.clims = plotData.clims;

            if ~isgraphics(this.figs.recon.handles.cur.plot) ||...
                this.figs.recon.update
                
                this.setTitleVariables('recon',...
                	{{this.figs.thirdAxLabel, this.figs.thirdAxCoor}});
                this.setAxesVar('recon', this.figs.depthAxLabel, this.figs.scanAxLabel);

                cla(this.figs.recon.handles.cur.ax)
                
                this.figs.recon.handles.cur.plot = ...
                    imagesc(this.figs.recon.handles.cur.ax,...
                    'XData', plotData.xData, 'Ydata', plotData.yData,...
                    'CData', plotData.cData, plotData.clims);
                
                axis(this.figs.recon.handles.cur.ax, 'tight')
                colorbar(this.figs.recon.handles.cur.ax);
                this.setLimsToPlot('recon');
                this.setStringsToPlot('recon');
                drawnow();
                this.figs.recon.update = false;
            else
                quickPlot(this, 'recon', plotData.xData, plotData.yData, plotData.cData)
            end
        end
        
    end
end
