classdef Graphics < handle
    %GRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        requests
        globalReq
        globalReqOld
        extUpdate
        numOfGraphs
        graphNames
    end
    
    methods(Static)
        function hS = createHandlesStruct()
            hS.fig        = NaN;
            hS.ax         = NaN;
            hS.plot       = NaN;
            hS.plotMarker = NaN;
            hS.title      = NaN;
            hS.leg        = NaN;
        end
    end
    
    methods
        function this = Graphics(plotNames)
            this.graphNames = plotNames;
            this.numOfGraphs = length(this.graphNames);
            this.extUpdate = true;
        end
        
        % Static Vars set functions
        function setType(this, gName, type)
            this.requests.(gName).type = type;
        end
        
        function setStrings(this, gName, title, xlabel, ylabel, legend)
            this.requests.(gName).strings.titleModel  = title;
            this.requests.(gName).strings.xlabel = xlabel;
            this.requests.(gName).strings.ylabel = ylabel;
            this.requests.(gName).strings.xlabelModel = xlabel;
            this.requests.(gName).strings.ylabelModel = ylabel;
            this.requests.(gName).strings.legendModel = legend;
        end        
        
        function setChAndPos(this, ch, zIdx, quant)
            this.requests.zIdx  = zIdx;
            this.requests.ch    = ch;
            this.requests.quant = quant;
        end
        
        function setMarkersEnable(this, gName, enDis)
             this.requests.(gName).markers.plotMark = enDis;
        end
        
        % Dynamic Vars set functions
        function setTitleVariables(this, gName, titleVars)
            for i = 1:length(this.requests.(gName).strings.titleModel)
                str{i} = sprintf(this.requests.(gName).strings.titleModel{i}, titleVars{i});
            end
            this.requests.(gName).strings.title  = str;
        end
        
        function setLimits(this, gName, xlims, ylims)
            this.requests.(gName).lims.xlims  = xlims;
            this.requests.(gName).lims.ylims  = ylims;
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
           
           this.requests.zIdx = this.globalReq.zIdx;
           this.requests.ch  = this.globalReq.ch;
           this.requests.quant = this.globalReq.quant;
        end
        
        function updateGraphicsConstruction(this)
            if (strcmp(this.globalReq.intExt, 'int')) 
                for i = 1:this.numOfGraphs
                    if strcmp(this.globalReqOld.intExt, 'ext') % backup current handles
                        this.requests.(this.graphNames{i}).handles.ext       = this.requests.(this.graphNames{i}).handles.cur;
                    elseif strcmp(this.globalReqOld.intExt, 'int')
                        this.requests.(this.graphNames{i}).handles.int       = this.requests.(this.graphNames{i}).handles.cur;
                    end
                    if this.globalReq.validStruct.(this.graphNames{i})
                        if ~isgraphics(this.requests.(this.graphNames{i}).handles.int.ax)  % check if you can use the old figures, if not, open a new one
                            this.requests.(this.graphNames{i}).handles.int.fig = figure();  
                            this.requests.(this.graphNames{i}).handles.int.ax  = axes();
                            this.requests.(this.graphNames{i}).handles.int.plot = NaN;
                        end
                        % and update cur with the new/old handles
                        this.requests.(this.graphNames{i}).handles.cur = this.requests.(this.graphNames{i}).handles.int;
                    else
                        this.requests.(this.graphNames{i}).handles.cur.fig  = NaN;
                        this.requests.(this.graphNames{i}).handles.cur.ax   = NaN;
                        this.requests.(this.graphNames{i}).handles.cur.plot = NaN;
                    end
                end
            else %if current is external
                for i = 1:this.numOfGraphs
                    if this.extUpdate % if its the first external copy the handles.
                        this.requests.(this.graphNames{i}).handles.ext.fig = this.globalReq.extH.(this.graphNames{i}).fig;
                        this.requests.(this.graphNames{i}).handles.ext.ax  = this.globalReq.extH.(this.graphNames{i}).ax;
                    end
                    if strcmp(this.globalReqOld.intExt, 'int') % if previous was internal, copy the handles
                        this.requests.(this.graphNames{i}).handles.int  = this.requests.(this.graphNames{i}).handles.cur;
                    elseif strcmp(this.globalReqOld.intExt, 'ext')
                        this.requests.(this.graphNames{i}).handles.ext  = this.requests.(this.graphNames{i}).handles.cur;
                    end
                    this.requests.(this.graphNames{i}).handles.cur = this.requests.(this.graphNames{i}).handles.ext;
                    this.requests.(this.graphNames{i}).fonts.titleSize  = this.globalReq.fonts.titleSize;
                    this.requests.(this.graphNames{i}).fonts.labelsSize = this.globalReq.fonts.labelsSize;
                    this.requests.(this.graphNames{i}).fonts.axisSize   = this.globalReq.fonts.axisSize;
                end
                this.extUpdate = false;
            end
        end   

        function setStringsToPlot(this, gName)
            set(this.requests.(gName).handles.cur.ax, 'FontSize', this.requests.(gName).fonts.axisSize)
            this.requests.(gName).handles.cur.title =...
                title(this.requests.(gName).handles.cur.ax, this.requests.(gName).strings.title, ...
                'FontSize', this.requests.(gName).fonts.titleSize);
            xlabel(this.requests.(gName).handles.cur.ax, this.requests.(gName).strings.xlabel,...
                 'FontSize', this.requests.(gName).fonts.labelsSize);
            ylabel(this.requests.(gName).handles.cur.ax, this.requests.(gName).strings.ylabel,...
                 'FontSize', this.requests.(gName).fonts.labelsSize);
            if length( this.requests.(gName).handles.cur.plot) > 1
                this.requests.(gName).handles.cur.legend =...
                    legend(this.requests.(gName).handles.cur.ax, this.requests.(gName).strings.legend);
            end
        end
        
        function setLimsToPlot(this, gName)
            if ~isempty(this.requests.(gName).lims.xlims)
                xlim(this.requests.(gName).handles.cur.ax, this.requests.(gName).lims.xlims)
            end
            if ~isempty(this.requests.(gName).lims.ylims)
                ylim(this.requests.(gName).handles.cur.ax, this.requests.(gName).lims.ylims)
            end
            if ~isempty(this.requests.(gName).lims.clims)
                caxis(this.requests.(gName).handles.cur.ax, this.requests.(gName).lims.clims)
            end
        end
    end
end