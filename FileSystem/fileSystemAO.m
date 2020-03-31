classdef fileSystemAO < handle
    
    properties
        uVars
        
        % Paths
        scanName    % As given from user (without date string)
        dirPath     % Directory in the disk in which this scan directory is opened
        
        projPath    % this scan directory (full path)
        resultsPath % relative to resPath
        figsPath    % relative to resPath
        
        extProjPath         % absolute path
        extProjResultsPath  % relative to projPath
        extProjFigsPath     % relative to projPath
        
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
        saveVars
        
        % Control Vars
        saveAny
        extProject
        
        dataFilename
    end
    
    methods (Static)
        function vars = uVarsCreate()
            vars.saveRawData        = false;
            vars.saveNetSignal      = false;
            vars.saveDemultiplexed  = false;
            vars.saveReshapedSignal = false;
            vars.saveFFT            = false;
            vars.savePhiChCmplx     = false;
            
            vars.saveResults        = false;
            vars.saveFigs           = false;
            vars.saveVars           = false; 
            
            vars.saveAny            = [];
            
            vars.scanName           = [];
            vars.dirPath            = [];
            
            vars.extProject         = false;
            vars.extProjPath        = [];
            vars.extProjResultsPath = [];
            vars.extProjFigsPath    = [];
        end
    end
    
    methods
        function this = fileSystemAO()
        end
        
        function setUserVars(this, uVars)
            this.uVars = uVars;
            
            this.saveRawData        = uVars.saveRawData;
            this.saveNetSignal      = uVars.saveNetSignal;
            this.saveDemultiplexed  = uVars.saveDemultiplexed;
            this.saveReshapedSignal = uVars.saveReshapedSignal;
            this.saveFFT            = uVars.saveFFT;
            this.savePhiChCmplx     = uVars.savePhiChCmplx;
            
            this.saveResults = uVars.saveResults;
            this.saveFigs    = uVars.saveFigs;
            this.saveVars    = uVars.saveVars;
            
            this.saveAny     = uVars.saveAny;
            
            this.scanName    = uVars.scanName;
            this.dirPath     = uVars.dirPath;
            
            this.extProject            = uVars.extProject;            
            this.extProjPath           = uVars.extProjPath;
            this.extProjResultsPath    = uVars.extProjResultsPath;
            this.extProjFigsPath       = uVars.extProjFigsPath;
            
           
        end
        
        function  vars = configFileSystem(this)
            if this.saveAny || this.saveVars
                if this.extProject
                    %in case acousto optics is called by another project
                    %that already created results directory, AO should
                    %save its results into the rawData directory of this
                    %project. In that case, the supervising project should
                    %turn the existingProject flag, and fill its directory
                    %in the projDirPath variable
                    
                    this.projPath       = this.extProjPath;
                    this.resultsPath    = this.extProjResultsPath;
                    this.figsPath       = this.extProjFigsPath;
                    
                else
                    %in case acousto optics is independant, it creates its
                    %own directories.
                    
                    dateStr          = strrep(datestr(datetime('now')), ':', '-');
                    
                    this.resultsPath = "Results";
                    this.figsPath    = "Figures";
                    
                    this.projPath   = sprintf("%s/%s-%s", this.dirPath,  dateStr, this.scanName);
                    resultsDir      = sprintf("%s/%s",    this.projPath, this.resultsPath);
                    figsDir         = sprintf("%s/%s",    this.projPath, this.figsPath);
                    
                    mkdir(this.projPath);
                    mkdir(resultsDir);
                    mkdir(figsDir);
                end
            end
            this.dataFilename = "Results.mat";
            vars = this.getFilesystemVars();
        end
        
        function vars = getFilesystemVars(this)
            vars.saveRawData        = this.saveRawData;
            vars.saveNetSignal      = this.saveNetSignal;
            vars.saveDemultiplexed  = this.saveDemultiplexed;
            vars.saveReshapedSignal = this.saveReshapedSignal;
            vars.saveFFT            = this.saveFFT;
            vars.savePhiChCmplx     = this.savePhiChCmplx;
            
            vars.saveResults = this.saveResults;
            vars.saveFigs    = this.saveFigs;            
            vars.saveVars    = this.saveVars;
            
            vars.saveAny = this.saveAny;
            
            vars.scanName           = this.scanName;
            vars.dirPath            = this.dirPath;
            
            vars.projPath           = this.projPath;
            vars.resultsPath        = this.resultsPath;
            vars.figsPath           = this.figsPath;
            
            vars.extProject         = this.extProject;
            vars.extProjPath        = this.extProjPath;
            vars.extProjResultsPath = this.extProjResultsPath;
            vars.extProjFigsPath    = this.extProjFigsPath;
        end
        
        function saveData(this, res)
            if this.saveAny
                if this.saveRawData
                    data.rawData = res.rawData;    
                end
                if this.saveNetSignal
                     data.netSignal = res.netSignal;
                end
                if this.saveDemultiplexed
                    data.deMultiplexed = res.deMultiplexed;
                end
                if this.saveReshapedSignal
                    data.reshapedSignal = res.reshapedSignal;
                end
                if this.saveFFT 
                    data.FFT = res.fftRes;
                end
                if this.savePhiChCmplx
                    data.phiChCmplx = res.phiChCmplx;
                    data.phiCh      = res.phiCh;
                end
                
                data.phiQuant = res.phiQuant;
                data.phi      = res.phi;
                data.phiStd   = res.phiStd;
                
                filename = sprintf("%s/%s/%s", this.projPath, this.resultsPath, this.dataFilename);
                save(filename, '-struct', 'data', '-v7.3');                
            end
        end
        
        function updateFileName(this, str)
            this.dataFilename = sprintf("%s.mat", str);
        end 
        
        function saveVarsToDisk(this, vars, path)
            % the path is relative to resDir given earlier.
            if this.saveVars
                fprintf("AOI: FS: Saving variables.\n");
                filename = sprintf("%s/%sVars.mat", this.projPath, path);
                save(filename, '-struct', 'vars', '-v7.3');
            end
        end

    end
end

