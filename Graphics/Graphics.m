classdef Graphics < handle
    %GRAPHICS Summary of this class goes here
    % Detailed explanation goes here
    
    properties
        figs
        
        uVars
%         uVarsOld
        
        extUpdate % set external handles only once on c'tor
        
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
            hS.text       = NaN;
        end
    end
    
    methods
        function this = Graphics()
%             this.graphNames = plotNames;
%             this.numOfGraphs = length(this.graphNames);
            this.extUpdate = true;
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
            this.figs.(gName).strings.updateLegend = true;
        end        
        
        function setChAndPos(this, ch, depthIdx, quant)
            this.figs.depthIdx  = depthIdx;
            this.figs.ch    = ch;
            this.figs.quant = quant;
        end
        
        function setMarkersEnable(this, gName, enDis)
             this.figs.(gName).markers.plotMark = enDis;
        end
        
        function setFonts(this, fonts)
            this.figs.fonts.type       = fonts.type;
            this.figs.fonts.titleSize  = fonts.titleSize;
            this.figs.fonts.labelsSize = fonts.labelsSize;
            this.figs.fonts.axisSize   = fonts.axisSize; 
        end
        
        % Dynamic Vars set functions
        function setTitleVariables(this, gName, titleVars)
            str = cell(1,length(this.figs.(gName).strings.titleModel));
            for i = 1:length(this.figs.(gName).strings.titleModel)
                str{i} = this.figs.(gName).strings.titleModel(i);
%                 for j=1:length(titleVars{i})
                    str{i} = sprintf(str{i}, titleVars{i}{:});
%                 end
            end
            this.figs.(gName).strings.title  = str;
        end
        
        function setAxesVar(this, gName, xName, yName)
           if ~isempty(xName)
           this.figs.(gName).strings.xlabel = ...
               sprintf(this.figs.(gName).strings.xlabelModel, xName);
           end
           if ~isempty(yName)
               this.figs.(gName).strings.ylabel = ...
                   sprintf(this.figs.(gName).strings.ylabelModel, yName);
           end
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
           this.uVars    = uVars;
        end
        
        function updateGraphicsConstruction(this)
            % NOTICE: in internal operation, this function should be called
            % in every run. in external operation, this function can be
            % called only once, and if called in every run will not take
            % any effect
            
            % 1. backup current handles
            for i = 1:this.numOfFigs
                if strcmp(this.figs.intExt, 'ext') 
                    this.figs.(this.figsNames{i}).handles.ext       = this.figs.(this.figsNames{i}).handles.cur;
                elseif strcmp(this.figs.intExt, 'int')
                    this.figs.(this.figsNames{i}).handles.int       = this.figs.(this.figsNames{i}).handles.cur;
                end
            end

            % 2. Manage construction
            %if requested is internal
            if (strcmp(this.uVars.intExt, 'int')) 
                for i = 1:this.numOfFigs
                    if this.uVars.validStruct.(this.figsNames{i})
                        if ~isgraphics(this.figs.(this.figsNames{i}).handles.int.ax) || this.uVars.reopenFigures % check if you can use the old figures, if not, open a new one
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
            else %if requested is external
                for i = 1:this.numOfFigs
                    if this.extUpdate % if it's the first external, copy the handles.
                        this.figs.(this.figsNames{i}).handles.ext.fig = this.uVars.extH.(this.figsNames{i}).fig;
                        this.figs.(this.figsNames{i}).handles.ext.ax  = this.uVars.extH.(this.figsNames{i}).ax;
%                         this.figs.(this.figsNames{i}).handles.cur     = this.figs.(this.figsNames{i}).handles.ext;
                    end
                    if ~strcmp(this.figs.intExt, 'ext') %if current is not external
                        this.figs.(this.figsNames{i}).handles.cur = this.figs.(this.figsNames{i}).handles.ext;
                    end
                end
                this.extUpdate = false;
                % NOTICE: External operation is independant in the valide
                % struct.
            end
            %3. update control variable
            this.figs.intExt = this.uVars.intExt;
            
            %4. turn on update for valid graphics
            for i = 1:this.numOfFigs
                if this.uVars.validStruct.(this.figsNames{i})
                    this.figs.(this.figsNames{i}).update = true;
                end
            end
            
        end   

        function setStringsToPlot(this, gName)
            set(this.figs.(gName).handles.cur.ax, 'FontSize', this.figs.fonts.axisSize)
            this.figs.(gName).handles.cur.title =...
                title(this.figs.(gName).handles.cur.ax, this.figs.(gName).strings.title, ...
                'FontSize', this.figs.fonts.titleSize);
            
            
            xlabel(this.figs.(gName).handles.cur.ax, this.figs.(gName).strings.xlabel,...
                 'FontSize', this.figs.fonts.labelsSize);
            ylabel(this.figs.(gName).handles.cur.ax, this.figs.(gName).strings.ylabel, ...
                 'FontSize', this.figs.fonts.labelsSize);
             
            if length( this.figs.(gName).handles.cur.plot) > 1 && this.figs.(gName).strings.updateLegend
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
        
        

    end
end