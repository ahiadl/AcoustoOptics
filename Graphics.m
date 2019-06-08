classdef Graphics < handle
    %GRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        requests
        globalReq
        globalReqOld
        numOfGraphs
        graphNames
    end
    
    methods(Static)
        
        function graphReq = getGraphicsRequest(plotNames)
            graphReq.pos     = 1;
            graphReq.ch      = 1;
            graphReq.intExt  = 'internal';
%             graphReq.numOfCh = 4;
%             graphReq.z = [];
%             graphReq.fSin = [];
            
            for i = 1:length(plotNames)  

%                 graphReq.(plotNames{i}).valid     = false;
                graphReq.(plotNames{i}).type      = []; 
 
                graphReq.(plotNames{i}).handles = Graphics.createHandlesStruct();

                graphReq.(plotNames{i}).strings.title       = [];
                graphReq.(plotNames{i}).strings.titleModel  = [];
                graphReq.(plotNames{i}).strings.xlabel      = [];
                graphReq.(plotNames{i}).strings.ylabel      = [];
                graphReq.(plotNames{i}).strings.legend      = [];
                graphReq.(plotNames{i}).strings.legendModel = [];

                graphReq.(plotNames{i}).dims.posDim  = [];
                graphReq.(plotNames{i}).dims.chDim   = [];
                graphReq.(plotNames{i}).dims.dataDim = [];

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
        
        function hS = createHandlesStruct()
            hS.fig        = NaN;
            hS.ax         = NaN;
            hS.plot       = NaN;
            hS.plotMarker = NaN;
            hS.title      = NaN;
            hS.leg        = NaN;
        end
        
        function gReq = createGraphicsRunVars(gNames) %for Reduced Operation
            gReq.ch  = 1;
            gReq.pos = 1;
            
            gReq.intExt = 'int';
            for i=1:length(gNames)
                gReq.validStruct.(gNames{i}) = false;
            end
            
            for i=1:length(gNames)
                gReq.extH.(gNames{i}) =  Graphics.createHandlesStruct();
            end
        end
        
    end
    
    methods
        function this = Graphics(plotNames)
            this.requests = Graphics.getGraphicsRequest(plotNames);
            this.graphNames = plotNames;
            this.numOfGraphs = length(this.graphNames);
        end
        
        % Static Vars set functions
        function setType(this, gName, type)
            this.requests.(gName).type = type;
        end
        
        function setStrings(this, gName, title, xlabel, ylabel, legend)
            this.requests.(gName).strings.titleModel  = title;
            this.requests.(gName).strings.xlabel = xlabel;
            this.requests.(gName).strings.ylabel = ylabel;
            this.requests.(gName).strings.legendModel = legend;
        end
        
        function setDims(this, gName, posDim, chDim, dataDim)
            this.requests.(gName).dims.posDim  = posDim;
            this.requests.(gName).dims.chDim   = chDim;
            this.requests.(gName).dims.dataDim = dataDim;
        end
        
        function setChAndPos(this, ch, pos)
            this.requests.pos     = ch;
            this.requests.ch      = pos;
        end
        
        function setMarkersEnable(this, gName, enDis)
             this.requests.(gName).markers.plotMark = enDis;
        end
        
        % Dynamic Vars set functions
        function setTitleVariables(this, gName, titleVars)
        	this.requests.(gName).strings.title  =...
                sprintf(this.requests.(gName).strings.titleModel, titleVars);
        end
        
        function setLimits(this, gName, xlims, ylims)
            this.requests.(gName).lims.xlims  = xlims;
            this.requests.(gName).lims.ylims  = ylims;
%             graphReq.(plotNames{i}).lims.colors = [];
        end
        
        function setLegendVariables(this, gName, legendVars)
            for i=1:length(legendVars)
                this.requests.(gName).strings.legend{i}  =...
                    sprintf(this.requests.(gName).strings.legendModel, legendVars(i));
            end
        end
        
        % Set Global Request functions
        function setGlobalReq(this, newGlobalReq)
           this.globalReqOld = this.globalReq;
           this.globalReq = newGlobalReq;
           
           this.updateGraphicsConstruction()
        end
        
        function updateGraphicsConstruction(this)
%             if isequal(globalReq, this.globalReq); return; end
            if (strcmp(this.globalReq.intExt, 'int')) %if current is internal
                    for i = 1:this.numOfGraphs
                        if this.globalReq.validStruct.(this.graphNames{i}) %this plot was asked
                            if ~isgraphics(this.requests.(this.graphNames{i}).handles.ax) % but its graphic is closed, reopen
                                this.requests.(this.graphNames{i}).handles.fig = figure();
                                this.requests.(this.graphNames{i}).handles.ax = axes();
                                
%                                 this.requests.(this.graphNames{i}).handles.plot
%                                 = plot(); %not necessary

%                                 this.requests.(this.graphNames{i}).handles.title =...
%                                     title(this.requests.(this.graphNames{i}).handles.ax,...
%                                           this.requests.(this.graphNames{i}).strings.title);
%                                 set(this.requests.(this.graphNames{i}).handles.title, 'FontSize', this.requests.(this.graphNames{i}).fonts.titleSize)      
%                                 this.requests.(this.graphNames{i}).handles.legend =...
%                                     legend(this.requests.(this.graphNames{i}).handles.ax,...
%                                           this.requests.(this.graphNames{i}).strings.legend);
                            end
                        elseif isgraphics(this.requests.(this.graphNames{i}).handles.ax) && ~(strcmp(this.globalReqOld.intExt, 'ext')) 
                               % this plot was not asked but the figure is still open and previos was not external
                                close(this.requests.(this.graphNames{i}).handles.fig)    % close the figure
                        end
                    end
            else %if current request is external
                for i = 1:this.numOfGraphs
                    if strcmp(this.globalReqOld.intExt, 'int') && isgraphics(this.requests.(this.graphNames{i}).handles.ax) %if previous was internal, close it
                            close(this.requests.(this.graphNames{i}).handles.fig) % close the figure
                    end
                    if this.globalReq.validStruct.(this.graphNames{i}) %this plot was asked
                            this.requests.(this.graphNames{i}).handles.fig    = this.globalReq.handles.fig;
                            this.requests.(this.graphNames{i}).handles.ax     = this.globalReq.handles.ax;
%                             this.requests.(this.graphNames{i}).handles.plot   = this.globalReq.handles.plot;
%                             this.requests.(this.graphNames{i}).handles.title  = this.globalReq.handles.title;
%                             this.requests.(this.graphNames{i}).handles.legend = this.globalReq.handles.legend;
                    end
                end
            end
        end
        
        % Plot Functions
        function plotSingleSignal(this, gName, xData, yData)
            if ~isgraphics(this.requests.(gName).handles.ax)
               return
            end
           
            if ~isgraphics(this.requests.(gName).handles.plot)
                switch this.requests.(gName).type
                   case 'stem'
                       this.requests.(gName).handles.plot =...
                           stem(this.requests.(gName).handles.ax,...
                                xData, yData);    
                   case 'plot'
                        this.requests.(gName).handles.plot =...
                           plot(this.requests.(gName).handles.ax,...
                                xData, yData); 
                end
                set(this.requests.(gName).handles.ax, 'FontSize', this.requests.(gName).fonts.axisSize)
                this.requests.(gName).handles.title =...
                    title(this.requests.(gName).handles.ax, this.requests.(gName).strings.title, ...
                    'FontSize', this.requests.(gName).fonts.titleSize);
                xlabel(this.requests.(gName).handles.ax, this.requests.(gName).strings.xlabel,...
                     'FontSize', this.requests.(gName).fonts.labelsSize);
                ylabel(this.requests.(gName).handles.ax, this.requests.(gName).strings.ylabel,...
                     'FontSize', this.requests.(gName).fonts.labelsSize);
            else
                set(this.requests.(gName).handles.plot, 'XData', xData, 'YData', gather(yData));
            end
        end
        
        function plotMultipleSignals(this, gName, xData, yData)
            if ~isgraphics(this.requests.(gName).handles.ax)
               return
            end
            
            chDim   = this.requests.(gName).dims.chDim;
            dataDim = this.requests.(gName).dims.dataDim;
            yData = permute(yData, [chDim, dataDim]);
            
            if ~isgraphics(this.requests.(gName).handles.plot)
                for i=1:size(yData, 1)
                    switch this.requests.(gName).type
                       case 'stem'
                           this.requests.(gName).handles.plot(i) =...
                               stem(this.requests.(gName).handles.ax,...
                                    xData, yData(i,:)); 
                       case 'plot'
                            this.requests.(gName).handles.plot(i) =...
                               plot(this.requests.(gName).handles.ax,...
                                    xData, yData(i,:)); 
                    end
                    hold(this.requests.(gName).handles.ax, 'on');
                end
                hold(this.requests.(gName).handles.ax, 'off');
                
                set(this.requests.(gName).handles.ax, 'FontSize', this.requests.(gName).fonts.axisSize)
                this.requests.(gName).handles.title =...
                    title(this.requests.(gName).handles.ax, this.requests.(gName).strings.title, ...
                    'FontSize', this.requests.(gName).fonts.titleSize);
                xlabel(this.requests.(gName).handles.ax, this.requests.(gName).strings.xlabel,...
                     'FontSize', this.requests.(gName).fonts.labelsSize);
                ylabel(this.requests.(gName).handles.ax, this.requests.(gName).strings.ylabel,...
                     'FontSize', this.requests.(gName).fonts.labelsSize);
                this.requests.(gName).handles.legend =...
                    legend(this.requests.(gName).handles.ax, this.requests.(gName).strings.legend);
            else
                for i=1:size(yData, 1)
                    set(this.requests.(gName).handles.plot(i), 'XData', xData, 'YData', gather(yData(i,:)) );
                end
            end
        end
        
        function plotSelectedSignals(this, gName, xData, yData)
            if ~isgraphics(this.requests.(gName).handles.ax)
               return
            end

            chDim   = this.requests.(gName).dims.chDim;
            posDim  = this.requests.(gName).dims.posDim;
            dataDim = this.requests.(gName).dims.dataDim;
            yData = permute(yData(), [chDim, dataDim, posDim]);

            pos = this.globalReq.pos;
            ch  = this.globalReq.ch;

            if ~isgraphics(this.requests.(gName).handles.plot)
                switch this.requests.(gName).type
                   case 'stem'
                       this.requests.(gName).handles.plot =...
                           stem(this.requests.(gName).handles.ax,...
                                xData, yData(ch, :, pos));    
                   case 'plot'
                        this.requests.(gName).handles.plot =...
                           plot(this.requests.(gName).handles.ax,...
                                xData, yData(ch, :, pos)); 
                end
                set(this.requests.(gName).handles.ax, 'FontSize', this.requests.(gName).fonts.axisSize)
                this.requests.(gName).handles.title =...
                    title(this.requests.(gName).handles.ax, this.requests.(gName).strings.title, ...
                    'FontSize', this.requests.(gName).fonts.titleSize);
                xlabel(this.requests.(gName).handles.ax, this.requests.(gName).strings.xlabel,...
                     'FontSize', this.requests.(gName).fonts.labelsSize);
                ylabel(this.requests.(gName).handles.ax, this.requests.(gName).strings.ylabel,...
                     'FontSize', this.requests.(gName).fonts.labelsSize);
            else
                set(this.requests.(gName).handles.plot, 'XData', xData, 'YData', gather(yData(ch, :, pos)));
            end
        end
        
        function plotFFT(this, gName, xData, yData, markerX, markerY)

            if ~isgraphics(this.requests.(gName).handles.ax)
               return
            end

            chDim   = this.requests.(gName).dims.chDim;
            posDim  = this.requests.(gName).dims.posDim;
            dataDim = this.requests.(gName).dims.dataDim;
            yData = permute(yData(), [chDim, dataDim, posDim]);

            pos = this.globalReq.pos;
            chNum = size(yData, chDim);
            if ~isgraphics(this.requests.(gName).handles.plot(1))
                for i=1:chNum
                   switch this.requests.(gName).type
                       case 'stem'
                           this.requests.(gName).handles.plot(i) =...
                               stem(this.requests.(gName).handles.ax,...
                                    xData, yData(i, :, pos));    
                       case 'plot'
                            this.requests.(gName).handles.plot(i) =...
                               plot(this.requests.(gName).handles.ax,...
                                    xData, yData(i, :, pos)); 
                   end
                   hold(this.requests.(gName).handles.ax, 'on');
                end
                
                if this.requests.(gName).markers.plotMark
                    for i=1:chNum
                       this.requests.(gName).handles.plotMarker(i) = ...
                            plot(this.requests.(gName).handles.ax,...
                                 markerX, markerY(i,:,pos), 'g+');
                       hold(this.requests.(gName).handles.ax, 'on');
                    end
                end
                hold(this.requests.(gName).handles.ax, 'off');
                set(this.requests.(gName).handles.ax, 'FontSize', this.requests.(gName).fonts.axisSize)
                this.requests.(gName).handles.title =...
                    title(this.requests.(gName).handles.ax, this.requests.(gName).strings.title, ...
                    'FontSize', this.requests.(gName).fonts.titleSize);
                xlabel(this.requests.(gName).handles.ax, this.requests.(gName).strings.xlabel,...
                     'FontSize', this.requests.(gName).fonts.labelsSize);
                ylabel(this.requests.(gName).handles.ax, this.requests.(gName).strings.ylabel,...
                     'FontSize', this.requests.(gName).fonts.labelsSize);
                 
                this.requests.(gName).handles.legend =...
                    legend(this.requests.(gName).handles.ax, this.requests.(gName).strings.legend);
                xlim(this.requests.(gName).handles.ax, this.requests.(gName).lims.xlims)

            else
                for i=1:size(yData,this.requests.(gName).dims.chDim) 
                    set(this.requests.(gName).handles.plot(i), 'XData', xData, 'YData', gather(yData(i, :, pos)));
                    if this.requests.(gName).markers.plotMark
                        set(this.requests.(gName).handles.plotMarker(i), 'XData', markerX, 'YData', gather(markerY(i,:,pos)));
                    end
                end
            end
        end
        
    end
end