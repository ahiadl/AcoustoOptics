classdef fileSystemS2D < fileSystem
    
    properties

        
    end
    
    methods (Static)
        function vars = uVarsCreate()
            vars.data2DFileNameModel  = [];
        end
    end
    
    methods
        function this = fileSystemS2D(subObjHandle)
            this@fileSystem(subObjHandle)
            this.uVars = fileSystemS2D.uVarsCreate();
            this.fsName = "2DScan";
        end
        
        function setUserVars(this, uVars)
            setUserVars@fileSystem(this, uVars)
            this.data2DFileNameModel = uVars.data2DFileNameModel;
            
            this.saveAny = this.saveAny || uVars.saveRawData        || uVars.saveNetSignal || uVars.saveDemultiplexed || ...
                           uVars.saveReshapedSignal || uVars.saveFFT       || uVars.savePhiChCmplx;  
            
        end
        
        function vars = configFileSystem(this)
            if this.saveAny || this.saveVars
                if this.extProject
                    
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
                    
                    this.data3DFilename = "Scan3DResults.mat";
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
%             vars.saveRawData        = this.saveRawData;
%             vars.saveNetSignal      = this.saveNetSignal;
%             vars.saveDemultiplexed  = this.saveDemultiplexed;
%             vars.saveReshapedSignal = this.saveReshapedSignal;
%             vars.saveFFT            = this.saveFFT;
%             vars.savePhiChCmplx     = this.savePhiChCmplx;
            
            vars.saveResults = this.saveResults;
            vars.saveFigs    = this.saveFigs;            
            vars.saveVars    = this.saveVars;
            
            vars.saveAny = this.saveAny;
            
            vars.scanName             = this.scanName;
            vars.dirPath              = this.dirPath;
            
            vars.data3DFileNameModel  = this.data3DFileNameModel;
            vars.data3DFilename       = this.data3DFilename;
            vars.data2DFileNameModel  = this.data2DFileNameModel;
            vars.data2DFilename       = this.data2DFilename;
            
            vars.varsFileNameModel    = this.varsFileNameModel;
            vars.varsFilename         = this.varsFilename;
            
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
        
        function save3DResultsToDisk(this, res)
            fprintf("S2D: FS: Saving results.\n");
            filename = sprintf("%s/%s%s", this.projPath, this.resultsPath, this.data3DFilename);
            save(filename, '-struct', 'res', '-v7.3');                
        end

        function save2DResultsToDisk(this, res)
            fprintf("S2D: FS: Saving results.\n");
            filename = sprintf("%s/%s%s", this.projPath, this.resultsPath, this.data2DFilename);
            save(filename, '-struct', 'res', '-v7.3');                
        end
        
        function saveVarsToDisk(this, vars, path)
            % the path is relative to resDir given earlier and contain / at the end.
            fprintf("S2D: FS: Saving variables.\n");
            filename = sprintf("%s/%s%s.mat", this.projPath, path, this.varsFilename);
            save(filename, '-struct', 'vars', '-v7.3');
        end
        
        function set3DDataFilenameModel(this, model)
            this.data3DFileNameModel = model;
        end
        
        function set2DDataFilenameModel(this, model)
            this.data2DFileNameModel = model;
        end
       
        function set3DDataFilenameVariables(this, stringVars)
            this.data3DFilename = sprintf(this.data3DFileNameModel, stringVars{:});
        end
        
        function set2DDataFilenameVariables(this, stringVars)
            this.data2DFilename = sprintf(this.data2DFileNameModel, stringVars{:});
        end
        
        function setVarsFilenameModel(this, model)
            this.varsFileNameModel = model;
        end
        
        function setVarsFilenameVariables(this, stringVars)
            this.varsFilename = sprintf(this.varsFileNameModel, stringVars{:});
        end

        function turnLogFileOn(this)
            if this.saveAny && ~this.extProject
                diary(sprintf("%s/log.txt", this.projPath));
                diary on
            end
        end
        
        function turnLogFileOff(this)
            if this.saveAny  && ~this.extProject
                diary(sprintf("%s/log.txt", this.projPath))
                diary off
            end
        end
    end
end

