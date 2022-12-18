classdef MTUA < handle
    % MTUA Summary of this class goes here
    % Measure Optical Transmission User Analysis
    
    properties
        i;      % disc2Idx
        j;      % disc1Idx
        k;      % contIdx
        hCs;    % handle to cs object
        aux;    % Auxilary variable for UA
        plots;  % figures handles
        uaData; % post processing data
        csVars;
    end
    
    methods
        function obj = MTUA()

        end
        
        function setCSHandle(this, hCs)
            this.hCs = hCs;
        end
        
        function setVars(this, vars)
            this.csVars = vars;
        end
        
        function updateScanIdx(this) 
            this.k = this.hCs.vars.idxs(1);
            this.j = this.hCs.vars.idxs(2);
            this.i = this.hCs.vars.idxs(3);
        end
        
        function initUAAux(this)
            nK = this.hCs.vars.scanSizeBin(1);
            nJ = this.hCs.vars.scanSizeBin(2);
            nI = this.hCs.vars.scanSizeBin(3);
            this.uaData.dcMatAvg = zeros(nK, nJ, nI);
            this.uaData.dcMatMax = zeros(nK, nJ, nI);
        end
        
        function init(this)
            this.updateScanIdx();
            this.initUAAux();
            this.initPlots();
        end
        
        function initPlots(this)
            this.plots.hFig = figure();
            this.plots.ax1 = subplot(2,2,1);
            this.plots.ax2 = subplot(2,2,2);
            this.plots.ax3 = subplot(2,2,3);
            
            xAxis = this.hCs.vars.disc1Vec;
            yAxis = this.hCs.vars.scanVecBin;
            
            this.plots.hP1 = imagesc(this.plots.ax1, xAxis, yAxis, this.uaData.dcMatAvg);
            xlabel(this.plots.ax1, "Scan Axis [mm]");
            ylabel (this.plots.ax1, "Disc Axis [mm]");
            colorbar (this.plots.ax1)
            axis(this.plots.ax1, 'equal', 'tight')
            title(this.plots.ax1, "2D DC");
            
            this.plots.hP2 = imagesc(this.plots.ax2, xAxis, yAxis, this.uaData.dcMatMax);
            xlabel(this.plots.ax2, "Scan Axis [mm]");
            ylabel (this.plots.ax2, "Disc Axis [mm]");
            axis(this.plots.ax2, 'equal', 'tight')
            colorbar(this.plots.ax2)
            title(this.plots.ax2, "2D Max");
            
            this.plots.hP3 = plot(this.plots.ax3, yAxis, zeros (1,length(yAxis)));
            xlabel(this.plots.ax3, "Scan Axis [mm]");
            title(this.plots.ax3, "1D DC");
            
        end
        
        function postProc(this, data)
            this.updateScanIdx();
            this.uaData.dcMatAvg(this.k, this.j, this.i) = squeeze(mean(data.curData, 5));
            this.uaData.dcMatMax(this.k, this.j, this.i) = squeeze(max(data.curData, [], 5));
        end
        
        function plotResults(this)
            this.updateScanIdx();
            set(this.plots.hP1, 'CData', this.uaData.dcMatAvg);
            set(this.plots.hP2, 'CData', this.uaData.dcMatMax);
            set(this.plots.hP3, 'YData', this.uaData.dcMatAvg(:,this.j, this.i));
            drawnow()
        end
        
    end
end

