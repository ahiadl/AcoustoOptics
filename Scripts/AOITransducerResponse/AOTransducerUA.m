classdef AOTransducerUA < handle
    % MTUA Summary of this class goes here
    % Measure Optical Transmission User Analysis
    
    properties
        i;      % disc2Idx
        j;      % disc1Idx
        k;      % contIdx
        hCs;    % handle to cs object
        plots;  % figures handles
        uaData; % post processing data
        vars;   
        csVars;
    end
    
    methods
        function obj = AOTransducerUA()

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
%             nK = this.hCs.vars.scanSizeBin(1);
%             nJ = this.hCs.vars.scanSizeBin(2);
%             nI = this.hCs.vars.scanSizeBin(3);
            
            this.vars.tVec = this.csVars.tVec;
            N = length(this.vars.tVec);
            fs = this.csVars.fs;
            this.vars.fBar = ((fs/N) *  ( (-N/2) : 1 : (N/2)-1 )) *1e-6;
            this.vars.fs = fs;
            this.vars.N = N;

            this.vars.fUSIdx = (this.vars.fBar == 1.25);
        end
        
        function init(this)
            this.updateScanIdx();
            this.initUAAux();
            this.initPlots();
        end
        
        function initPlots(this)
            this.plots.hFig = figure();
            this.plots.ax1 = subplot(1,2,1);
            this.plots.ax2 = subplot(1,2,2);

            hold(this.plots.ax1, 'on')

            xAxis = this.hCs.vars.tVec*1e6;
            yAxis = zeros(1, length(xAxis));

            this.plots.hP11 = plot(this.plots.ax1, xAxis, yAxis);
            this.plots.hP12 = plot(this.plots.ax1, xAxis, yAxis);
            xlabel(this.plots.ax1, "t [\mu s]");
            ylabel(this.plots.ax1, "Signal [KPa]");
            title(this.plots.ax1, "Temporal Signal");
            
            hold(this.plots.ax2, 'on')
            this.plots.hP21 = plot(this.plots.ax2, this.vars.fBar, yAxis);
            this.plots.hP22 = plot(this.plots.ax2, this.vars.fBar(this.vars.fUSIdx), 0, '+g');
            xlabel(this.plots.ax2, "frequency [MHz]");
            ylabel (this.plots.ax2, "Power Spectrum [AU]");
            title(this.plots.ax2, "Fourier Transform");
        end
        
        function postProc(this, data)
            this.updateScanIdx();
            this.uaData.sig  = (squeeze(data.curData)*1e3 /837 )* 1e3;
            this.uaData.acSig =  this.uaData.sig - mean(this.uaData.sig);
            this.uaData.envelope = envelope(this.uaData.acSig, 10, 'peak');
            this.uaData.fft = abs(fftshift(fft(this.uaData.acSig))).^2;
            this.uaData.fftNorm = this.uaData.fft./max(this.uaData.fft);
        end
        
        function plotResults(this)
            this.updateScanIdx();
            set(this.plots.hP11, 'YData', this.uaData.acSig);
            set(this.plots.hP12, 'YData', this.uaData.envelope);
            set(this.plots.hP21, 'YData', this.uaData.fftNorm);
            set(this.plots.hP22, 'YData', this.uaData.fftNorm(this.vars.fUSIdx));
            drawnow()
        end
        
    end
end

