classdef dvGraphics < Graphics
    %SACN2DGRAPHICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data
    end
    
    methods (Static)
  
       function figsNames = getGraphicsNames()
            figsNames = {'deepView'; 'deepViewLog'};
       end 
        
       function figs  = createGraphicsVars()    
            % General Scan Parameters
            figs.depthIdx   = 1;
            figs.firstAxis  = [];
            figs.depthAxis  = [];
            figs.repeats    = [];

            figs.firstAxLabel   = 'Y';
            figs.dxFirstAx      = [];
            
            figs.numOfTimeFrames = 0;
            figs.timeFrames      = [];
            
            % Current Position
            figs.curPosFirst  = 0;
            figs.curIdxFirst  = 0;
            figs.curTime      = 0;
            figs.curIdxTime   = 0;
            figs.curRep       = 0;            
            
            %Specific Figures parameters
            figs.intExt  = 'int';
            figsNames    = dvGraphics.getGraphicsNames();
            
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
           figs.zIdx   = 0;
           
           figs.intExt      = 'int';
            
           figsNames = scan3DGraphics.getGraphicsNames();

           for i=1:length(figsNames)
               figs.validStruct.(figsNames{i}) = false;
           end

           for i=1:length(figsNames)
               figs.extH.(figsNames{i}) =  Graphics.createHandlesStruct();
           end
           
           figs.fonts.type       = [];
           figs.fonts.titleSize  = 14;
           figs.fonts.labelsSize = 14;
           figs.fonts.axisSize   = 14;
       end
       
       function uVars = createOwnerVars() %for Reduced Operation
            uVars = scan3DGraphics.createUserVars();
            
            uVars.firstAxis  = [];
            uVars.depthAxis  = [];
            uVars.timeFrames = [];
            uVars.repeats    = [];
            
            uVars.numOfTimeFrames = 0;
            
            uVars.firstAxLabel   = 'Y';

            % Current Position
            uVars.curPosFirst  = 0;
            uVars.curIdxFirst  = 0;
            uVars.curIdxTime   = 0;
            uVars.curRep       = 0; 

       end
    end
    
    methods
        function this = dvGraphics()
            this@Graphics()
            this.figsNames = dvGraphics.getGraphicsNames();
            this.figs      = dvGraphics.createGraphicsVars();
            this.uVars     = dvGraphics.createOwnerVars();
            this.setGraphicsStaticVars();
            this.numOfFigs = length(this.figsNames);
        end
 
        % Set vars Functions
        function setGraphicsStaticVars(this)
            this.setType('deepView', 'plot');
            this.setStrings('deepView', "Deep View Measurement", "%s[mm]", "\\phi", "T=%d");
            
            this.setType('deepViewLog', 'plot');
            this.setStrings('deepViewLog', "Deep View Measurement (dB)", "%s[mm]", "\\phi(dB)", "T=%d");
        end
                 
        function setGraphicsScanVars(this)
            this.figs.repeats    = this.uVars.repeats;
            this.figs.firstAxis  = this.uVars.firstAxis;
            this.figs.depthIdx   = this.uVars.depthIdx;
            this.figs.depthAxis  = this.uVars.depthAxis;
            
            this.figs.timeFrames      = this.uVars.timeFrames;
            this.figs.numOfTimeFrames = this.uVars.numOfTimeFrames;
            
            % reset data arrays:
            this.initDataArray();
            
            this.figs.firstAxLabel   = this.uVars.firstAxLabel;
            this.figs.dxFirstAx      = abs(this.figs.firstAxis(1) - this.figs.firstAxis(2));
            
            this.figs.validStruct = this.uVars.validStruct;
            %extH and intExt flag should not be transferred at this point.
            %it is updated by the updateGraphicsConstruction function
            
            if this.figs.repeats > 1
                this.setType('deepView', 'errorbar');
            else
                this.setType('deepView', 'plot');
            end
            this.updateGraphicsConstruction();
            this.figs.deepView.handles.cur.plot = zeros(1,this.figs.numOfTimeFrames);
        end
        
        function updateTimeAndInd(this, curTime, curIdx)
            this.figs.curIdxTime = curIdx;
            this.figs.curTime    = curTime;
        end
        
        function updateS2DCurPosAndIdx(this, pns)
            this.figs.curRep       = pns.curPosIdx(1);
            this.figs.curIdxFirst  = pns.curPosIdx(2);
            this.figs.curPosFirst  = pns.curPos(1);
        end
        
        function setDepthIdx(this, idx)
           this.figs.depthIdx= idx;
           this.updatePlots()
        end
        
        %Set Data functions
        function initDataArray(this)
            rstPhi    = zeros(length(this.figs.firstAxis), this.figs.numOfTimeFrames, length(this.figs.depthAxis), this.uVars.repeats);
            rstPhiAvg = zeros(length(this.figs.firstAxis), this.figs.numOfTimeFrames, length(this.figs.depthAxis));
            
            this.data.phi       = rstPhi;
            this.data.phiAvg    = rstPhiAvg;
            this.data.phiAvgStd = rstPhiAvg;
        end
        
        function set1DData(this, data)
            this.data.phi(this.figs.curIdxFirst, this.figs.curIdxTime, :, this.figs.curRep) = data.phi;
        end
        
        function setAvg1DData (this, data)
            this.data.phiAvg(this.figs.curIdxFirst, this.figs.curIdxTime, :)    = data.phiAvg;
            this.data.phiAvgStd(this.figs.curIdxFirst, this.figs.curIdxTime, :) = data.phiAvgStd;
        end
        
        function setLoadedData(this, data)
            this.data.phi = data.phi;
            this.data.phiStd = data.phiStd;
            this.data.phiAvg = data.phiAvg;
            this.data.phiAvgStd = data.phiAvgStd;
            
            this.setLoadedDataClims(data);
        end
        
        % Extract data functions
        function data = extractDeepViewFromData(this, figName)
            data.xData = squeeze(this.figs.firstAxis(1:this.figs.curIdxFirst));
            switch figName
                case 'deepView'
                    data.yData = squeeze(this.data.phiAvg(1:this.figs.curIdxFirst, this.figs.curIdxTime, this.figs.depthIdx));
                case 'deepViewLog'
                    data.yData = db(squeeze(this.data.phiAvg(1:this.figs.curIdxFirst, this.figs.curIdxTime, this.figs.depthIdx)));
            end
        end
        
        function notifyNewTimeFrame(this)
           this.figs.newTimeFrame = true; 
        end
        
        % Display Functions
        function updatePlots(this)
            this.dispDeepView('deepView');
            this.dispDeepView('deepViewLog');
        end
        
        function dispDeepView(this, figName)
            % yData  - (xLen, yLen, zLen, rep)
            if ~isgraphics(this.figs.(figName).handles.cur.ax)
               return
            end

            plotData = this.extractDeepViewFromData();
            
            if ~isgraphics(this.figs.(figName).handles.cur.plot(this.figs.curIdxTime)) ||...
                this.figs.newTimeFrame
                
                this.setAxesVar(figName, this.figs.firstAxLabel, []);
                
                switch this.figs.(figName).type
                    case 'plot'
                        this.figs.(figName).handles.cur.plot(this.figs.curIdxTime) = ...
                            plot(this.figs.(figName).handles.cur.ax,...
                            plotData.xData, plotData.yData, '-+'); 
                end
                hold(this.figs.(figName).handles.cur.ax, 'on');
                this.setLegendVariables(figName, this.figs.timeFrames(1:this.figs.curIdxTime));
                this.setLimsToPlot(figName)
                this.setStringsToPlot(figName);
                drawnow();
                this.figs.newTimeFrame = false;
            else
                quickPlot(this, figName, plotData.xData, plotData.yData)
            end
        end
        
        function quickPlot(this, gName, xData, yData) 
            switch this.figs.(gName).type
                case 'plot'    
                    set(this.figs.(gName).handles.cur.plot(this.figs.curIdxTime),...
                        'XData', xData,...
                        'YData', yData);                    
            end 
            this.setLimsToPlot(gName)
            set(this.figs.(gName).handles.cur.title, 'String', this.figs.(gName).strings.title);
            pause(0.01);
            drawnow();
        end

    end
    
end