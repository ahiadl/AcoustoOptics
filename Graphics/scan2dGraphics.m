classdef scan2DGraphics < Graphics
    %SACN2DGRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data
    end
    
    methods (Static)
  
       function figsNames = getGraphicsNames()
            figsNames = {'curMainAxis'; 'curMainAxisAvg'; 'curMainPlane'; 'curMainPlaneAvg'; 'nav'; 'navAvg'};
       end 
        
       function figs = getGraphicsVars()    
            % General Scan Parameters
            figs.xAxis   = [];
            figs.yAxis   = [];
            figs.zAxis   = [];
            figs.repeats = [];
            
            figs.scanFirstAx   = 'Y';
            figs.scanSecondAx  = 'X';
            figs.scanMainPlane = 'YZ';
            
            figs.firstAxis  = [];
            figs.secondAxis = [];
            figs.xlim       = [];
            figs.dx         = [];
            
            figs.zIdx    = 1; 
            figs.zPos    = 0; % the z value related to zIdx;
            
            % Current Position
            figs.curXpos  = 0;
            figs.curYpos  = 0;
            
            figs.curXidx  = 0;
            figs.curYidx  = 0;
            figs.curRep   = 0;
            figs.curFirstPos  = 0;
            figs.curSecondPos = 0;
            
            % Navigator Position
            figs.navAx    = [];  % navigating on this axes
            figs.navPlane = []; % displaying these planes

            figs.navPos  = 0;
            figs.navIdx  = 0;
            figs.navRep  = 0;
            
            figs.navXaxis = [];
            figs.navYaxis = [];
            
            figs.navXname = 'X';
            figs.navYname = 'Y';
            
            %Specific Figures parameters
            figs.intExt  = 'int';
            figsNames = scan2DGraphics.getGraphicsNames();
            
            figs.fonts.type       = [];
            figs.fonts.titleSize  = 18;
            figs.fonts.labelsSize = 18;
            figs.fonts.axisSize   = 18;
            
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
        
       function uVars = createGraphicsUserVars() %for Reduced Operation
            uVars.zIdx = 1;
            
            uVars.useQuant = true;
            uVars.repeats  = 1;
            uVars.xAxis    = [];
            uVars.yAxis    = [];
            uVars.zAxis    = [];
            
            uVars.scanFirstAx = 'Y';
            
            figsNames = scan2DGraphics.getGraphicsNames();
            
            uVars.intExt = [];
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
       
       function navVars = createNavigatorVars()
            navVars.navAx   = 'Z';
            navVars.navIdx  = 0;
            navVars.navRep  = 0;
       end
    end
    
    methods
        function this = scan2DGraphics()
            this@Graphics()
            this.figsNames = scan2DGraphics.getGraphicsNames();
            this.figs      = scan2DGraphics.getGraphicsVars();
            this.uVars     = scan2DGraphics.createGraphicsUserVars();
            this.setGraphicsStaticVars();
            this.numOfFigs = length(this.figsNames);
        end
 
        % Set vars Functions
        function setGraphicsStaticVars(this)
            % curMainAxis
            this.setType('curMainAxis', 'stem');
            this.setStrings('curMainAxis', "Current Main Axis %s: (R,%s,Z) = (%d, %.2f, %.2f)", "%s[mm]", "\\phi[v]", []);
            
            % curMainAxisAvg
            this.setType('curMainAxisAvg', 'errorbar');
            this.setStrings('curMainAxisAvg', "Current Main Axis %s (Averaged): (R, %s, Z) = (%d, %.2f, %.2f)", "%s[mm]",  "\\phi[v]", []);
            
            % curMainPlane
            this.setType('curMainPlane', 'imagesc');
            this.setStrings('curMainPlane', "Current Main Plane %s: (R, %s) = (%d, %.2f)", "%s [mm]", "%s[mm]", []);
            
            % curMainPlaneAvg
            this.setType('curMainPlaneAvg', 'imagesc');
            this.setStrings('curMainPlaneAvg', "Current Main Plane %s (Averaged): (R, %s) = (%d, %.2f)", "%s [mm]", "%s[mm]", []);
           
            % navPlane
            this.setType('nav', 'imagesc');
            this.setStrings('nav', "Navigator Plane %s: (R, %s) = (%d, %.2f)", "%s[mm]", "%s[mm]", []);
            
            % navPlaneRep
            this.setType('navAvg', 'imagesc');
            this.setStrings('navAvg', "Navigator Plane %s (Averaged): (R, %s) = (%d, %.2f)", "%s[mm]", "%s[mm]", []);
            
            
%             this.figs.intExt      = staticVars.intExt;
%             this.figs.validStruct = staticVars.validStruct;
%             this.figs.extH        = staticVars.extH;
        end
                 
        function setGraphicsScanVars(this)
            this.figs.useQuant = this.uVars.useQuant;
            this.figs.repeats  = this.uVars.repeats;
            this.figs.xAxis    = this.uVars.xAxis;
            this.figs.yAxis    = this.uVars.yAxis;
            this.figs.zAxis    = this.uVars.zAxis*1e3;
            
            % reset data arrays:
            rstPhi    = zeros(length(this.figs.xAxis), length(this.figs.yAxis), length(this.figs.zAxis), this.uVars.repeats);
            rstPhiAvg = zeros(length(this.figs.xAxis), length(this.figs.yAxis), length(this.figs.zAxis), this.uVars.repeats);
            
            this.data.phi       = rstPhi;
            this.data.phiStd    = rstPhi;
            this.data.phiAvg    = rstPhiAvg;
            this.data.phiAVgStd = rstPhiAvg;
            
            if this.uVars.zIdx > length(this.figs.zAxis)
                this.figs.zIdx = 1;
            else
                this.figs.zIdx = this.uVars.zIdx;
            end
            this.figs.zPos = this.figs.zAxis(this.figs.zIdx);
            
            this.figs.scanFirstAx = this.uVars.scanFirstAx;
            
            this.figs.validStruct = this.uVars.validStruct;
            
            switch this.figs.scanFirstAx
                case 'X'
                    this.figs.scanSecondAx  = 'Y';
                    this.figs.scanMainPlane = 'XZ';
                    this.figs.firstAxis     = this.figs.xAxis;
                    this.figs.secondAxis    = this.figs.yAxis;
                    this.figs.dx            = abs(this.figs.xAxis(1) - this.figs.xAxis(2));
                    this.figs.xlim          = [this.figs.xAxis(end)-this.figs.dx, this.figs.xAxis(1)+this.figs.dx];
                case 'Y'
                    this.figs.scanSecondAx  = 'X';
                    this.figs.scanMainPlane = 'YZ';
                    this.figs.firstAxis     = this.figs.yAxis;
                    this.figs.secondAxis    = this.figs.xAxis;
                    this.figs.dx            = abs(this.figs.yAxis(1) - this.figs.yAxis(2));
                    this.figs.xlim          = [this.figs.yAxis(1)- this.figs.dx, this.figs.yAxis(end)+this.figs.dx];
            end
            
            if this.figs.useQuant
                this.setType('curMainAxis', 'errorbar');
            else
                this.setType('curMainAxis', 'stem');
            end
            
            if this.figs.repeats > 1
                this.setType('curMainAxisAvg', 'errorbar');
            else
                this.setType('curMainAxisAvg', 'stem');
            end

            this.figs.fonts = this.uVars.fonts;
            
            this.resetNavigator();
        end

        function updateCurPosAndIdx(this, pns)
            this.figs.curRep       = pns.curPosIdx(1);
            this.figs.curXidx      = pns.curPosIdx(2);
            this.figs.curYidx      = pns.curPosIdx(3);
            this.figs.curXpos      = pns.curPos(1);
            this.figs.curYpos      = pns.curPos(2);
            this.figs.curFirstPos  = pns.curScan(1);
            this.figs.curSecondPos = pns.curScan(2);
        end

        function setNavVars(this, navVars)
            this.figs.navIdx =  navVars.navIdx;
            this.figs.navRep =  navVars.navRep;
            if ~strcmp( this.figs.navAx, navVars.navAx) || this.figs.navUpdate
                this.figs.nav.update = true;
                this.figs.navAvg.update = true;
                this.figs.navAx = navVars.navAx;
                switch this.figs.navAx
                    case 'X'
                        this.figs.navPlane = 'YZ';
                        this.figs.navAxis  = this.figs.xAxis;
                        this.figs.navXaxis = this.figs.yAxis;
                        this.figs.navYaxis = this.figs.zAxis;
                        this.figs.navXname = 'Y';
                        this.figs.navYname = 'Z';
                    case 'Y'
                        this.figs.navPlane = 'XZ';
                        this.figs.navAxis  = this.figs.yAxis;
                        this.figs.navXaxis = this.figs.xAxis;
                        this.figs.navYaxis = this.figs.zAxis;
                        this.figs.navXname = 'X';
                        this.figs.navYname = 'Z';
                    case 'Z'
                        this.figs.navPlane = 'XY';
                        this.figs.navAxis  = this.figs.zAxis;
                        this.figs.navXaxis = this.figs.xAxis;
                        this.figs.navYaxis = this.figs.yAxis;
                        this.figs.navXname = 'X';
                        this.figs.navYname = 'Y';
                end
            end
            this.figs.navPos   = this.figs.navAxis(this.figs.navIdx);
        end
        
        function resetNavigator(this)
            navVars.navAx = 'Z';
            navVars.navIdx = 1;
            navVars.navRep = 1;
            this.figs.navUpdate = true;
            this.setNavVars(navVars);
            this.dispNavPlane();
            this.dispNavPlaneAvg();
        end
        
        function setData(this, phi, phiStd, phiAvg, phiAvgStd)
           this.data.phi       = phi;
           this.data.phiStd    = phiStd;
           this.data.phiAvg    = phiAvg;
           this.data.phiAvgStd = phiAvgStd;           
        end
        
        % Extract data functions
        function data = extractMainAxisFromData(this, isAvg)
            if ~isAvg
                if strcmp(this.figs.scanFirstAx, 'Y')
                    data.yData  = this.data.phi(this.figs.curXidx, :, this.figs.zIdx, this.figs.curRep);
                    data.stdVec = this.data.phiStd(this.figs.curXidx, :, this.figs.zIdx, this.figs.curRep);
                elseif strcmp(this.figs.scanFirstAx, 'X')
                    data.yData  = this.data.phi(:, this.figs.curYidx, this.figs.zIdx, this.figs.curRep);
                    data.stdVec = this.data.phiStd(:, this.figs.curYidx, this.figs.zIdx, this.figs.curRep);
                end
            else
                if strcmp(this.figs.scanFirstAx, 'Y')
                    data.yData  = this.data.phiAvg (this.figs.curXidx, :, this.figs.zIdx);
                    data.stdVec = this.data.phiStd(this.figs.curXidx, :, this.figs.zIdx);
                elseif strcmp(this.figs.scanFirstAx, 'X')
                    data.yData  = this.data.phiAvg (:, this.figs.curYidx, this.figs.zIdx);
                    data.stdVec = this.data.phiStd(:, this.figs.curYidx, this.figs.zIdx);
                end
            end
        end
        
        function data = extractMainPlaneFromData(this, isAvg)
            if ~isAvg
                data.clims = [min(min(min(min(this.data.phi)))), max(max(max(max(this.data.phi))))];
                if strcmp(this.figs.scanFirstAx, 'Y')
                    data.cData  = permute(this.data.phi(this.figs.curXidx, :, :, this.figs.curRep), [3,2,1,4]);
                    data.xData  = this.figs.yAxis;
                elseif strcmp(this.figs.scanFirstAx, 'X')
                    data.cData  = permute(this.data.phi(:, this.figs.curYidx, :, this.figs.curRep), [3,1,2,4]);                
                    data.xData  = this.figs.xAxis;
                end
            else
                data.clims = [min(min(min(min(this.data.phiAvg)))), max(max(max(max(this.data.phiAvg))))];
                if strcmp(this.figs.scanFirstAx, 'Y')
                    data.cData  = permute(this.data.phiAvg(this.figs.curXidx, :, :), [3,2,1]);
                    data.xData  = this.figs.yAxis;
                elseif strcmp(this.figs.scanFirstAx, 'X')
                    data.cData  = permute(this.data.phiAvg(:, curYidx, :), [3,1,2]);                
                    data.xData  = this.figs.xAxis;
                end
            end
            
            if data.clims(1) == data.clims(2)
                data.clims(2) = data.clims(1)+ 0.1;
            end
        end
        
        function data = extractNavPlaneFromData(this, isAvg)
            if ~isAvg
                data.clims = [min(min(min(min(this.data.phi)))), max(max(max(max(this.data.phi))))];
                switch this.figs.navAx
                    case 'X'
                        data.cData  = permute(this.data.phi(this.figs.navIdx, :, :, this.figs.navRep), [3,2,1,4]);
                    case 'Y'     
                        data.cData  = permute(this.data.phi(:, this.figs.navIdx, :, this.figs.navRep), [3,1,2,4]);
                    case 'Z' 
                        data.cData  = permute(this.data.phi(:, :, this.figs.navIdx, this.figs.navRep), [2,1,3,4]); 
                end
            else
                data.clims = [min(min(min(min(this.data.phiAvg)))), max(max(max(max(this.data.phiAvg))))];
                switch this.figs.navAx
                    case 'X'
                        data.cData  = permute(this.data.phiAvg(this.figs.navIdx, :, :), [3,2,1]);
                    case 'Y'     
                        data.cData  = permute(this.data.phiAvg(:, this.figs.navIdx, :), [3,1,2]);
                    case 'Z' 
                        data.cData  = permute(this.data.phiAvg(:, :, this.figs.navIdx), [2,1,3]); 
                end
            end
            if (data.clims(1) == data.clims(2)) 
                data.clims(2) = data.clims(1)+0.1; 
            end
        end
        
        % Display Functions
        function dispCurMainAxis(this)
            % yData - (xLen, yLen, zLen, r)
            % stdMat - (xLen, yLen, zLen, r)
            if ~isgraphics(this.figs.curMainAxis.handles.cur.ax)
               return
            end
            
            this.setTitleVariables('curMainAxis',...
                                   {{this.figs.scanFirstAx, this.figs.scanSecondAx, this.figs.curRep, ...
                                     this.figs.curSecondPos, this.figs.zPos}})
            
            %extract relevent data
            plotData = this.extractMainAxisFromData(false);

            if ~isgraphics(this.figs.curMainAxis.handles.cur.plot) ||...
               this.figs.curMainAxis.update
                
                this.setAxesVar('curMainAxis', this.figs.scanFirstAx, 'Fluence');
                this.setLimits('curMainAxis', this.figs.xlim, []);
                
                switch this.figs.curMainAxis.type
                    case 'stem'
                        this.figs.curMainAxis.handles.cur.plot = ...
                            stem(this.figs.curMainAxis.handles.cur.ax,...
                            this.figs.firstAxis, plotData.yData); 
                    case 'errorbar'
                        this.figs.curMainAxis.handles.cur.plot = ...
                            errorbar(this.figs.curMainAxis.handles.cur.ax,...
                            this.figs.firstAxis, plotData.yData, plotData.stdVec);
                end
                
                this.setLimsToPlot('curMainAxis'); %apply limits to the plot
                this.setStringsToPlot('curMainAxis'); %allocate title and labels handles
                drawnow();
                this.figs.curMainAxis.update = false;
            else
                quickPlot(this, 'curMainAxis', this.figs.firstAxis, plotData.yData, [], plotData.stdVec)
            end
        end
        
        function dispCurMainAxisAvg(this)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.curMainAxisAvg.handles.cur.ax)
               return
            end
            this.setTitleVariables('curMainAxisAvg', {{this.figs.scanFirstAx, ...
                                                       this.figs.scanSecondAx, this.figs.curRep,...
                                                       this.figs.curSecondPos, this.figs.zPos}});            
            %extract relevent data
            plotData = this.extractMainAxisFromData(true);
            
            if ~isgraphics(this.figs.curMainAxisAvg.handles.cur.plot) ||...
                this.figs.curMainAxisAvg.update
                
                this.setAxesVar('curMainAxisAvg', this.figs.scanFirstAx, 'Fluence')
                this.setLimits('curMainAxisAvg', this.figs.xlim, [])
                switch this.figs.curMainAxisAvg.type
                    case 'stem'
                        this.figs.curMainAxisAvg.handles.cur.plot = ...
                            stem(this.figs.curMainAxisAvg.handles.cur.ax,...
                            this.figs.firstAxis, plotData.yData); 
                    case 'errorbar'
                        this.figs.curMainAxisAvg.handles.cur.plot = ...
                            errorbar(this.figs.curMainAxisAvg.handles.cur.ax,...
                            this.figs.firstAxis, plotData.yData, plotData.stdVec);
                end
                this.setLimsToPlot('curMainAxisAvg')
                this.setStringsToPlot('curMainAxisAvg');
                drawnow();
                this.figs.curMainAxisAvg.update = false;
            else
                quickPlot(this, 'curMainAxisAvg', this.figs.firstAxis, plotData.yData, [], plotData.stdVec)
            end
        end
        
        function dispCurMainPlane(this)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.curMainPlane.handles.cur.ax)
               return
            end

            this.setTitleVariables('curMainPlane', {{ this.figs.scanMainPlane,...
                                                      this.figs.scanSecondAx, this.figs.curRep,...
                                                      this.figs.curSecondPos}});
            
            plotData = this.extractMainPlaneFromData(false);
                                                  
            this.figs.curMainPlane.lims.clims = plotData.clims;

            if ~isgraphics(this.figs.curMainPlane.handles.cur.plot) ||...
                this.figs.curMainPlane.update
                
                cla(this.figs.curMainPlane.handles.cur.ax)
                this.setAxesVar('curMainPlane', this.figs.scanFirstAx, 'Z');
                
                this.figs.curMainPlane.handles.cur.plot = ...
                    imagesc(this.figs.curMainPlane.handles.cur.ax,...
                    'XData', plotData.xData, 'Ydata', this.figs.zAxis, 'CData', plotData.cData,...
                    plotData.clims);
                
                %title
                axis(this.figs.curMainPlane.handles.cur.ax, 'tight')
                colorbar(this.figs.curMainPlane.handles.cur.ax);
                this.setLimsToPlot('curMainPlane');
                this.setStringsToPlot('curMainPlane');
                drawnow();
                this.figs.curMainPlane.update = false;
            else
                quickPlot(this, 'curMainPlane', plotData.xData, this.figs.zAxis, plotData.cData)
            end
        end
        
        function dispCurMainPlaneAvg(this)
            % yData - (xLen, yLen, zLen)
            % stdMat - (xLen, yLen, zLen)
            if ~isgraphics(this.figs.curMainPlaneAvg.handles.cur.ax)
               return
            end

            this.setTitleVariables('curMainPlane', {{ this.figs.scanMainPlane, this.figs.scanSecondAx,...
                                                      this.figs.curRep, this.figs.curSecondPos}});
                                                                 
            plotData = this.extractMainPlaneFromData(true);
            this.figs.curMainPlaneAvg.lims.clims = plotData.clims;

            if ~isgraphics(this.figs.curMainPlaneAvg.handles.cur.plot) ||...
                this.figs.curMainPlane.update
                
                cla(this.figs.curMainPlaneAvg.handles.cur.ax)
                this.setAxesVar('curMainPlaneAvg', this.figs.scanSecondAx, 'Z');
                
                this.figs.curMainPlaneAvg.handles.cur.plot = ...
                    imagesc(this.figs.curMainPlaneAvg.handles.cur.ax,...
                    'XData', plotData.xData, 'Ydata', this.figs.zAxis, 'CData', plotData.cData,...
                    plotData.clims);
                
                axis(this.figs.curMainPlaneAvg.handles.cur.ax, 'tight')
                colorbar(this.figs.curMainPlaneAvg.handles.cur.ax);
                this.setLimsToPlot('curMainPlaneAvg')
                this.setStringsToPlot('curMainPlaneAvg');
                
                drawnow();
                this.figs.curMainPlaneAvg.update = false;
            else
                quickPlot(this, 'curMainPlaneAvg', plotData.xData, this.figs.zAxis, plotData.cData)
            end
        end

        function dispNavPlane(this)
            % yData  - (xLen, yLen, zLen, rep)
            if ~isgraphics(this.figs.nav.handles.cur.ax)
               return
            end
            
            this.setTitleVariables('nav', {{ this.figs.navPlane, this.figs.navAx, ...
                                             this.figs.navRep, this.figs.navPos}});

            plotData = this.extractNavPlaneFromData(false);
            this.figs.nav.lims.clims = plotData.clims;
            
            if ~isgraphics(this.figs.nav.handles.cur.plot) ||...
                this.figs.nav.update
                
                cla(this.figs.nav.handles.cur.ax)
                this.setAxesVar('nav', this.figs.navXname, this.figs.navYname);
                
                this.figs.nav.handles.cur.plot = ...
                    imagesc(this.figs.nav.handles.cur.ax, 'XData', this.figs.navXaxis,...
                    'Ydata', this.figs.navYaxis, 'CData', plotData.cData,...
                    plotData.clims);
                
                axis(this.figs.nav.handles.cur.ax, 'tight')
                colorbar(this.figs.nav.handles.cur.ax);
                this.setLimsToPlot('nav')
                this.setStringsToPlot('nav');
                drawnow();
                this.figs.nav.update = false;
            else
                quickPlot(this, 'nav', this.figs.navXaxis, this.figs.navYaxis, plotData.cData)
            end
        end
        
        function dispNavPlaneAvg(this)
            % yData  - (xLen, yLen, zLen, rep)
            if ~isgraphics(this.figs.navAvg.handles.cur.ax)
               return
            end
            this.setTitleVariables('navAvg', {{ this.figs.navPlane, this.figs.navAx,...
                                                this.figs.curRep,   this.figs.navPos}});

            plotData = this.extractNavPlaneFromData(true);
            this.figs.nav.lims.clims = plotData.clims;
                                        
            if ~isgraphics(this.figs.navAvg.handles.cur.plot) ||...
                this.figs.navAvg.update
                
                cla(this.figs.navAvg.handles.cur.ax)
                this.setAxesVar('navAvg', this.figs.navXname, this.figs.navYname);
                
                this.figs.navAvg.handles.cur.plot = ...
                    imagesc(this.figs.navAvg.handles.cur.ax, 'XData', this.figs.navXaxis,...
                    'Ydata', this.figs.navYaxis, 'CData', plotData.cData, plotData.clims);
                
                axis(this.figs.navAvg.handles.cur.ax, 'tight')
                colorbar(this.figs.navAvg.handles.cur.ax);
                this.setLimsToPlot('navAvg')
                this.setStringsToPlot('navAvg');
                drawnow();
                this.figs.navAvg.update = false;
            else
                quickPlot(this, 'navAvg', this.figs.navXaxis, this.figs.navYaxis, plotData.cData)
            end    
        end
        
        function quickPlot(this, gName, xData, yData, cData, stdDev) 
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
                        'YPositiveDelta', stdDev,...
                        'YNegativeDelta', stdDev);
                case 'imagesc'
                    set(this.figs.(gName).handles.cur.plot,...
                        'XData', xData,...
                        'YData', yData,...
                        'CData', cData);
                    this.setLimsToPlot(gName)
            end 
            set(this.figs.(gName).handles.cur.title, 'String', this.figs.(gName).strings.title);
            pause(0.01);
            drawnow();
        end

    end
    
end