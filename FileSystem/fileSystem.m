classdef fileSystem < handle
    
    properties
        uVars
        
        % Paths
        dirPath    % directory in the disk in which this scan directory is opened
        resDir     % this scan directory (full path) 
        rawDataDir % resDir\rawData - full path
        figsDir    % resDir\Figures - full path
        scanName   % as given from user (without date string)
        
        % Acousto Optics reconstruction saving indication 
        saveRawData
        saveNetSignal
        saveDemultiplexed
        saveReshapedSignal
        saveFFT
        savePhiChCmplx
        
        % Full scan saving indication
        saveFigs
        saveResults
        
        % Control Vars
        saveAny
    end
    
    methods
        function this = fileSystem()
            
        end
        
        function setUserVars(this, uVars)
            this.uVars = uVars;
            
            this.saveRawData        = uVars.saveRawData;
            this.saveNetSignal      = uVars.saveNetSignal;
            this.saveDemultiplexed  = uVars.saveDemultiplexed;
            this.saveReshapedSignal = uVars.saveReshapedSignal;
            this.saveFFT            = uVars.saveFFT;
            this.savePhiChCmplx     = uVars.savePhiChCmplx;
            
            this.saveFigs           = uVars.saveFigs;
            this.saveResults        = uVars.saveResults;
            
            this.saveAny = this.saveRawData        || this.saveNetSignal || this.saveDemultiplexed || ...
                           this.saveReshapedSignal || this.saveFFT       || this.savePhiChCmplx    || ...
                           this.saveFigs           || this.saveResults;
            
            if this.saveAny
                this.scanName   = uVars.scanName;
                dateStr         = strrep(datestr(datetime('now')), ':', '-');
                this.dirPath    = uVars.resDirPath;
                this.resDir     = sprintf("%s/%s-%s",this.dirPath, dateStr, this.scanName);
                this.rawDataDir = sprintf("%s/rawData", this.resDir);
                this.figsDir    = sprintf("%s/Figures", this.resDir);
                
                mkdir(this.resDir);
                mkdir(this.rawDataDir);
                mkdir(this.figsDir);
            end
        end
  
        function saveData(this, res, vars)
            if this.saveAny
                if this.saveRawData
                    data.rawData = res.rawData;    
                end
                if this.saveNetSignal
                    
                end
                if this.saveDemultiplexed
                    
                end
                if this.saveReshapedSignal
                    data.reshapedSignal = res.reshapedSignal;
                end
                if this.saveFFT 
                    data.FFT = res.FFT;
                end
                if this.savePhiChCmplx
                    data.phiChCmplx = data.phiChCmplx;
                end
                saveRawDataStructure(data, vars);
            end
        end
        
        function saveResultsToDisk(this, res, vars)
            if this.saveResults
                    save(sprintf("%s/%s", this.resDir, 'Vars.mat'),    'vars', '-v7.3');
                    save(sprintf("%s/%s", this.resDir, 'Results.mat'), 'res',  '-v7.3');    
            end   
        end
        
        function saveFigsToDisk(this)
            % TODO: complete
        end

    end
end

