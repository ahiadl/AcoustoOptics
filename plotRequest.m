classdef plotRequest < handle
    %PLOTREQUEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        request
        type
        pos
        ch
        
        values
        handles
        lims
        strings
    end
    
    methods
        
        function this = plotRequest(type, tit, xlabel, ylabel, legStr, fontSize)
            this.request = false;          
            this.type    = type; 
            
            this.values.pos     = 1;
            this.values.ch      = 1;
            this.values.chDim   = 3;
            this.values.legendValues = [];          
            this.values.posDim  = 1;
            
            this.strings.title  = tit;
            this.strings.xlabel = xlabel;
            this.strings.ylabel = ylabel;
            this.strings.legend = legStr;
            this.strings.titleFontSize  = fontSize(1);
            this.strings.labelsFontSize = fontSize(2);
            this.strings.axesFontSize   = fontSize(3);
           
            this.lims.xlims  = [];
            this.lims.ylims  = [];
            this.lims.colors = [];
            
            this.handles.fig  = NaN;
            this.handles.ax   = NaN;
            this.handles.plot = NaN;
            this.handles.tit  = NaN;
        end

%         function ax = getAx(this)
%             ax = this.handles.ax;
%         end
        
        function request = getIsValid(this)
            request = this.request;
        end
        
%         function setParams(this, request, fig, ax, ch, pos)
%             this.request = request;
%             this.handles.ax = ax;
%             this.ch = ch;
%             this.pos = pos;
%             this.handles.fig = fig;
%         end
        
%         function [ax, ch, pos, plot] = getRequest(this)
%             ax  = this.handles.ax;
%             ch  = this.ch;
%             pos = this.pos;
%             plot   = this.handles.plot;
%         end
        
        function h = getHandle(this)
           h = this.handles.h; 
        end
        
        function setHandle(this, plot)
           this.handles.plot = plot; 
        end
        
        function fig = getFigHandle(this)
           fig = this.handles.fig; 
        end
        
        function type = getType(this)
            type = this.type;
        end
        
        function h = getHandles(this)
           h = this.handles;
        end
        
        function str = getStrings(this)
            str = this.strings;
        end
        
            
    end
end

