classdef Graphics < handle
    %GRAPHICS Summary of this class goes here
    % Detailed explanation goes here
    
    properties
        figs
        
        uVars
        uVarsOld
        
        extUpdate % what this is for?
        
        numOfFigs
        figsNames
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
        function this = Graphics()
%             this.graphNames = plotNames;
%             this.numOfGraphs = length(this.graphNames);
            this.extUpdate = true; %what this is for?
        end
        
        % Static Vars set functions
        function setType(this, gName, type)
            this.figs.(gName).type = type;
        end
        
        function setStrings(this, gName, title, xlabel, ylabel, legend)
            this.figs.(gName).strings.titleModel  = title;
            this.figs.(gName).strings.xlabel = xlabel;
            this.figs.(gName).strings.ylabel = ylabel;
            this.figs.(gName).strings.xlabelModel = xlabel;
            this.figs.(gName).strings.ylabelModel = ylabel;
            this.figs.(gName).strings.legendModel = legend;
        end        
        
        function setChAndPos(this, ch, zIdx, quant)
            this.figs.zIdx  = zIdx;
            this.figs.ch    = ch;
            this.figs.quant = quant;
        end
        
        function setMarkersEnable(this, gName, enDis)
             this.figs.(gName).markers.plotMark = enDis;
        end
        
        % Dynamic Vars set functions
        function setTitleVariables(this, gName, titleVars)
            str = cell(1,length(this.figs.(gName).strings.titleModel));
            for i = 1:length(this.figs.(gName).strings.titleModel)
                str{i} = sprintf(this.figs.(gName).strings.titleModel{i}, titleVars{i});
            end
            this.figs.(gName).strings.title  = str;
        end
        
        function setLimits(this, gName, xlims, ylims)
            this.figs.(gName).lims.xlims  = xlims;
            this.figs.(gName).lims.ylims  = ylims;
        end
        
        function setLegendVariables(this, gName, legendVars)
            for i=1:length(legendVars)
                this.figs.(gName).strings.legend{i}  =...
                    sprintf(this.figs.(gName).strings.legendModel, legendVars(i));
            end
        end
        
        % Set Global Request functions
        function setUserVars(this, uVars)
           this.uVarsOld = this.uVars;
           this.uVars    = uVars;
           
           this.setChAndPos(uVars.ch, uVars.zIdx, uVars.quant);
%            this.figs.zIdx     = this.globalReq.zIdx;
%            this.requests.ch    = this.globalReq.ch;
%            this.requests.quant = this.globalReq.quant;
        end
        
        function updateGraphicsConstruction(this)
            if (strcmp(this.uVars.intExt, 'int')) 
                for i = 1:this.numOfFigs
                    if strcmp(this.uVarsOld.intExt, 'ext') % backup current handles
                        this.figs.(this.figsNames{i}).handles.ext       = this.figs.(this.figsNames{i}).handles.cur;
                    elseif strcmp(this.uVarsOld.intExt, 'int')
                        this.figs.(this.figsNames{i}).handles.int       = this.figs.(this.figsNames{i}).handles.cur;
                    end
                    if this.uVars.validStruct.(this.figsNames{i})
                        if ~isgraphics(this.figs.(this.figsNames{i}).handles.int.ax)  % check if you can use the old figures, if not, open a new one
                            this.figs.(this.figsNames{i}).handles.int.fig = figure();  
                            this.figs.(this.figsNames{i}).handles.int.ax  = axes();
                            this.figs.(this.figsNames{i}).handles.int.plot = NaN;
                        end
                        % and update cur with the new/old handles
                        this.figs.(this.figsNames{i}).handles.cur = this.figs.(this.figsNames{i}).handles.int;
                    else
                        this.figs.(this.figsNames{i}).handles.cur.fig  = NaN;
                        this.figs.(this.figsNames{i}).handles.cur.ax   = NaN;
                        this.figs.(this.figsNames{i}).handles.cur.plot = NaN;
                    end
                end
            else %if current is external
                for i = 1:this.numOfFigs
                    if this.extUpdate % if its the first external copy the handles.
                        this.figs.(this.figsNames{i}).handles.ext.fig = this.uVars.extH.(this.figsNames{i}).fig;
                        this.figs.(this.figsNames{i}).handles.ext.ax  = this.uVars.extH.(this.figsNames{i}).ax;
                    end
                    if strcmp(this.uVarsOld.intExt, 'int') % if previous was internal, copy the handles
                        this.figs.(this.figsNames{i}).handles.int  = this.figs.(this.figsNames{i}).handles.cur;
                    elseif strcmp(this.uVarsOld.intExt, 'ext')
                        this.figs.(this.figsNames{i}).handles.ext  = this.figs.(this.figsNames{i}).handles.cur;
                    end
                    this.figs.(this.figsNames{i}).handles.cur = this.figs.(this.figsNames{i}).handles.ext;
                    this.figs.(this.figsNames{i}).fonts.titleSize  = this.uVars.fonts.titleSize;
                    this.figs.(this.figsNames{i}).fonts.labelsSize = this.uVars.fonts.labelsSize;
                    this.figs.(this.figsNames{i}).fonts.axisSize   = this.uVars.fonts.axisSize;
                end
                this.extUpdate = false;
            end
        end   

        function setStringsToPlot(this, gName)
            set(this.figs.(gName).handles.cur.ax, 'FontSize', this.figs.(gName).fonts.axisSize)
            this.figs.(gName).handles.cur.title =...
                title(this.figs.(gName).handles.cur.ax, this.figs.(gName).strings.title, ...
                'FontSize', this.figs.(gName).fonts.titleSize);
            xlabel(this.figs.(gName).handles.cur.ax, this.figs.(gName).strings.xlabel,...
                 'FontSize', this.figs.(gName).fonts.labelsSize);
            ylabel(this.figs.(gName).handles.cur.ax, this.figs.(gName).strings.ylabel,...
                 'FontSize', this.figs.(gName).fonts.labelsSize);
            if length( this.figs.(gName).handles.cur.plot) > 1
                this.figs.(gName).handles.cur.legend =...
                    legend(this.figs.(gName).handles.cur.ax, this.figs.(gName).strings.legend);
            end
        end
        
        function setLimsToPlot(this, gName)
            if ~isempty(this.figs.(gName).lims.xlims)
                xlim(this.figs.(gName).handles.cur.ax, this.figs.(gName).lims.xlims)
            end
            if ~isempty(this.figs.(gName).lims.ylims)
                ylim(this.figs.(gName).handles.cur.ax, this.figs.(gName).lims.ylims)
            end
            if ~isempty(this.figs.(gName).lims.clims)
                caxis(this.figs.(gName).handles.cur.ax, this.figs.(gName).lims.clims)
            end
        end
        
%         TODO: complete
%         function popOutFigures(this)

%         end
    end
end