classdef Graphics < handle
    %GRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        requests
        uVars
        uVarsOld
    end
    
    methods(Static)
        
        function graphReq = getGraphicsRequest(plotNames)
            for i = 1:length(plotNames)  
                graphReq.(plotNames{i}).internal  = true;
                graphReq.(plotNames{i}).valid     = false;
                graphReq.(plotNames{i}).type      = []; 

                graphReq.(plotNames{i}).handles.fig    = [];
                graphReq.(plotNames{i}).handles.ax     = [];
                graphReq.(plotNames{i}).handles.plot   = [];
                graphReq.(plotNames{i}).handles.title  = [];
                graphReq.(plotNames{i}).handles.leg    = [];

                graphReq.(plotNames{i}).strings.title  = [];
                graphReq.(plotNames{i}).strings.xlabel = [];
                graphReq.(plotNames{i}).strings.ylabel = [];
                graphReq.(plotNames{i}).strings.legend = [];

                graphReq.(plotNames{i}).values.pos     = 1;
                graphReq.(plotNames{i}).values.ch      = 1;
                graphReq.(plotNames{i}).values.posDim  = [];
                graphReq.(plotNames{i}).values.chDim   = [];
                graphReq.(plotNames{i}).values.dataDim = [];

                graphReq.(plotNames{i}).markers.plotMark = false;
                graphReq.(plotNames{i}).markers.markX    = [];
                graphReq.(plotNames{i}).markers.markY    = []; 

                graphReq.(plotNames{i}).lims.xlims  = [];
                graphReq.(plotNames{i}).lims.ylims  = [];
                graphReq.(plotNames{i}).lims.colors = [];

                graphReq.(plotNames{i}).fonts.type       = [];
                graphReq.(plotNames{i}).fonts.titleSize  = 18;
                graphReq.(plotNames{i}).fonts.labelsSize = 18;
                graphReq.(plotNames{i}).fonts.axisSize   = 18;
            end
        end
    end
    
    methods
        function this = Graphics(plotNames)
            this.requests = Graphics.getGraphicsRequest(plotNames);
        end
        
        function setGraphicsRequests(this, gReq)
           this.uVarsOld = this.uVars;
           this.uVars = gReq.vars;
           this.setPlotRequests(gReq.plotRequests)  
        end
        
        function setPlotRequests(this, newPlotReq)
            if ~isempty(this.plotRequests)
                names = fieldnames(this.plotRequests);
                for i=1:length(names)
                    this.plotRequests.(names(i)).close();
                    this.plotRequests = rmfield(this.plotRequests,names(i));
                end
            end
%             names = fieldnames(newPlotReq);
%             for i=1:length(names)
            this.plotRequests =  newPlotReq;
%             end
        end
        
        function valid = isRequest(this, prName)
            valid = false;
            if ~isfield(this.plotRequests, prName); return; end
            valid = this.plotRequests.(prName).getIsValid();
        end
        
        function directGraphics(this)
            fprintf("------- Preparing Algorithm Graphics Requests ----------\n")
            if this.uVars.internal % if the current request is internal - 
                ch    = this.uVars.ch;
                pos   = this.uVars.pos;
                names = fieldnames(this.plotRequests);
                for i=1:length(names)
                    if this.uVarsOld.internal && this.uVarsOld.mask(i) % if previous was internal
                        close(this.plotRequests.(names{i}).getFigHandle()); %close old figs
                    end
                    if this.uVars.mask(i)
                        fh = figure(); ax = axes(); %open new figs
                        this.plotRequests.(names{i}).setParams(...
                                             this.uVars.mask(i), fh, ax, ch, pos)
                    end
                end
            else % if the current request is external
                names = fieldnames(this.plotRequests);
                for i=1:length(names)
                    if this.uVarsOld.internal && this.uVarsOld.mask(i) % if previous was internal
                        close(this.graphics.plotRequests.getFigHandle()); %close old figs
                    end
                end
                this.plotRequests = this.graphics.uVars.externalPlotRequests;
            end
            pause(0.1);
%             this.algo.setPlotRequests(this.graphics.plotRequests)
        end
        
        

        
        
        
        function displayResults(this, vars)
            % x should be a vector
            % y should be [pos,data,ch]
            
            % vars:
            % prName
            % yData
            % xData
            % titleVars



            for i = 1:size(vars.yData, pr.values.chDim)
                switch pr.getType()
                    case 'errBar'
                        set(h(i),...
                            'XData', vars.xData,...
                            'YData', vars. yData(pr.pos, :, i),...
                            'YNegativeDelta', vars.yStd(pr.pos, :, i),...
                            'YPositiveDelta', vars.yStd(pr.pos, :, i) );
                    case 'stem'
                        set(h(i),...
                            'XData', vars.xData,...
                            'YData', vars.yData(pr.pos, :, i));
                    case 'plot'
                        set(h(i),...
                            'XData', vars.xData,...
                            'YData', vars.yData(pr.pos, :, i));
                end
            end
            
            titleStr = sprintf(strings.title, vars.titleVars);
            title(ax, 'Title', titleStr);
            xlabel(ax, strings.xStr);
            ylabel(ax, strings.yStr);
            
            if size(vars.yData,3) > 1
                for i = 1:size(vars.legDim)
                    legStr{i} = sprintf(pr.strings.leg, pr.legendValues(i));
                end 
            end
            
            legend(ax, legStr);
            pause(0.01);
        end
        
        function plotSingleSignal(this, prName, xData, yData, titleArgs)
            pr = this.plotRequests.(prName);
            if ~pr.request; return; end
            h = pr.getHandles();
            if ~isgraphics(h.ax); return; end
            
            string = pr.getStrings();

            set(h.plot, 'XData', xData, 'YData', yData);
            set(h.title, 'String', sprintf(strings.title, titleArgs));
          
        end
        
        
        function plotFFT(this, xData, yData, xlims, ylims)
            
            
        end
        
        
        function plotChSignal(this)
        
            
        end
        
        
        function plotStatistics(this)
            
            
        end
        
        function plot2DColor(this)
             
        
        end
        
    end
end