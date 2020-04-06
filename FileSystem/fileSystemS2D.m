classdef fileSystemS2D < handle
    
    properties
        uVars
        
        % Paths
        scanName    % As given from user (without date string)
        dirPath     % Directory in the disk in which this scan directory is opened
        
        projPath    % this scan directory (full path)
        rawDataPath % AO rersults.  relative to resPath
        resultsPath % Scan results. relative to resPath
        figsPath    % relative to resPath
        
        extProjPath         % absolute path
        extProjRawDataPath  % relative to projPath
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
        
        dataFileNameModel
        varsFileNameModel
        dataFilename
        varsFilename
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
            vars.dataFileNameModel  = [];
            vars.varsFileNameModel  = [];
            
            vars.extProject         = false;
            vars.extProjPath        = [];
            vars.extProjRawDataPath = [];
            vars.extProjResultsPath = []; %relative to extProjPath
            vars.extProjFigsPath    = []; %relative to extProjPath
        end
    end
    
    methods
        function this = fileSystemS2D()
            this.uVars = fileSystemS2D.uVarsCreate();
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
            
            this.scanName          = uVars.scanName;
            this.dirPath           = uVars.dirPath;
            this.dataFileNameModel = uVars.dataFileNameModel;
            this.varsFileNameModel = uVars.varsFileNameModel;
            
            this.extProject            = uVars.extProject;            
            this.extProjPath           = uVars.extProjPath;
            this.extProjRawDataPath    = uVars.extProjRawDataPath;
            this.extProjResultsPath    = uVars.extProjResultsPath; 
            this.extProjFigsPath       = uVars.extProjFigsPath;    %relative to extProjPath
        end
        
        function vars = configFileSystem(this)
            if this.saveAny || this.saveVars
                if this.extProject
                    %in case acousto optics is called by another project
                    %that already created results directory, AO should
                    %save its results into the rawData directory of this
                    %project. In that case, the supervising project should
                    %turn the existingProject flag, and fill its directory
                    %in the projDirPath variable
                    
                    this.projPath       = this.extProjPath;
                    this.rawDataPath    = this.extProjRawDataPath;
                    this.resultsPath    = this.extProjResultsPath;
                    this.figsPath       = this.extProjFigsPath;      
                else
                    %in case acousto optics is independant, it creates its
                    %own directories.
                    
                    dateStr          = strrep(datestr(datetime('now')), ':', '-');
                    this.rawDataPath = "Results/";
                    this.resultsPath = "";
                    this.figsPath    = "Figures/";
                    
                    this.projPath   = sprintf("%s/%s-%s", this.dirPath,  dateStr, this.scanName);
                    rawDataDir      = sprintf("%s/%s",    this.projPath, this.rawDataPath);
                    resultsDir      = sprintf("%s/%s",    this.projPath, this.resultsPath);
                    figsDir         = sprintf("%s/%s",    this.projPath, this.figsPath);
                    
                    this.dataFilename = "ScanResults.mat";
                    this.varsFilename = "ScanVars.mat";
                    
                    mkdir(this.projPath);
                    mkdir(rawDataDir);
                    mkdir(resultsDir);
                    mkdir(figsDir);
                    
                end
            end
            
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
            vars.dataFileNameModel  = this.dataFileNameModel;
            vars.dataFilename       = this.dataFilename;
            vars.varsFileNameModel  = this.varsFileNameModel;
            vars.varsFilename       = this.varsFilename;
            
            vars.projPath           = this.projPath;
            vars.rawDataPath        = this.rawDataPath;
            vars.resultsPath        = this.resultsPath;
            vars.figsPath           = this.figsPath;
            
            vars.extProject         = this.extProject;
            vars.extProjPath        = this.extProjPath;
            vars.extProjRawDataPath = this.extProjRawDataPath;
            vars.extProjResultsPath = this.extProjResultsPath;
            vars.extProjFigsPath    = this.extProjFigsPath;
        end
        
        function saveResultsToDisk(this, res)
            fprintf("S2D: FS: Saving results.\n");
            filename = sprintf("%s/%s%s", this.projPath, this.resultsPath, this.dataFilename);
            save(filename, '-struct', 'res', '-v7.3');                
        end

        function saveVarsToDisk(this, vars, path)
            % the path is relative to resDir given earlier and contain / at the end.
            fprintf("S2D: FS: Saving variables.\n");
            filename = sprintf("%s/%s%s.mat", this.projPath, path, this.varsFilename);
            save(filename, '-struct', 'vars', '-v7.3');
        end
        
        function setDataFilenameModel(this, model)
            this.dataFileNameModel = model;
        end
        
        function setVarsFilenameModel(this, model)
            this.varsFileNameModel = model;
        end
        
        function setDataFilenameVariables(this, stringVars)
            this.dataFilename = sprintf(this.dataFileNameModel, stringVars{:});
        end
        
        function setVarsFilenameVariables(this, stringVars)
            this.varsFilename = sprintf(this.varsFileNameModel, stringVars{:});
        end

        function turnLogFileOn(this)
            if this.saveAny
                diary(sprintf("%s/log.txt", this.projPath));
                diary on
            end
        end
        
        function turnLogFileOff(this)
            if this.saveAny
                diary(sprintf("%s/log.txt", this.projPath))
                diary off
            end
        end
    end
end

