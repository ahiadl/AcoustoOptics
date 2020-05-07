classdef fileSystem < handle
    
    properties
        %New Used Variables
        uVars
        hSubFS      % handles to sub filesystems
        
        fsName      % for logging prints
        objName  % for saving results 
        
        dirPath     % Directory in the disk in which this project directory is created
        projPath    % This project directory (full path) [resDir/projName]
        extProjPath % Project path given bu external FS object
        resPath     % path to the subObjResults directory   
        
        projName    % project Name, relevant in case of no owner FS.
        resDirName  % The folder under projPath in which the results should be saved
        
        saveFigs    % should figures be saved?
        saveResults % should results be saved?
        saveVars    % should variables be saved?
        
        saveAny           % is there any action of saving in this object?
        saveAnySub        % is there any action of saving in subFS?
        saveAnyTot        % is there any action of saving in total?
        extProject        % is this an external project?
        stackAllSubObjRes % should all subObj Results not be divided into folders (for direct owner of acostoOptics)?
        useExtVarsPath    % shold Vars be saved in this projPath or in user given path?
        extVarsPath       % external path for Vars
        extVarsPrefix     % name of the vars file
        
        scanIdentifierPrefix       % current scan of the owner object:              Z=1
        scanIdentifierSuffixModel  % current object addition to the scanIdentifier: T=%d
        scanIdentifierModel        % current scan model:                            Z=1-T=%d
        scanIdentifier             % current scan including vars:                   Z=1-T=1
        
        defaultScanIdentifierPrefix
        
        fileNameModel              % current scan model of the filename:            AOS2D-Z=1-T=1.mat
    end
    
    methods (Static)
        function vars = uVarsCreate()
            % Should be received from the fileSystem of the related object
            % e.g. fileSystemAO
            vars.saveResults = false;
            vars.saveFigs    = false; 
            
            vars.dirPath     = [];
            vars.projName    = [];
            vars.resDirName  = [];
            
            vars.extProject        = false;
            vars.stackAllSubObjRes = false;
            vars.useExtVarsPath    = false;
            vars.extVarsPath       = [];
        end
    end
    
    methods
        function this = fileSystem(subObjHandle)
            if nargin>1
                this.hSubFS = subObjHandle;
            else
                this.hSubFS = [];
            end
            this.projPath = [];
            this.resDirName = "Results";
            this.stackAllSubObjRes = false;
        end 
        
        function setUserVars(this, uVars)
            this.uVars = uVars;
            
            this.saveResults = uVars.saveResults;
            this.saveFigs    = uVars.saveFigs;
            this.saveVars    = this.saveResults || this.saveFigs; 
            
            this.saveAny = this.saveVars || this.saveResults || this.saveFigs;
            
            this.projName = uVars.projName;
            this.dirPath  = uVars.dirPath;

            this.extProject            = uVars.extProject; 
            
            uVars.stackAllSubObjRes = false; %not relevant for ao
            uVars.useExtVarsPath    = false;
            uVars.extVarsPath       = [];
        end
        
        function configFileSystem(this)  
            if this.calcSaveAnyTot() 
                % If it is an external project, the project path should be
                % given already at this point using the updateFileSystem 
                % method of the owner FS.
                if ~this.extProject
                    dateStr       = strrep(datestr(datetime('now')), ':', '-');
                    this.projPath = sprintf("%s/%s-%s", this.dirPath, dateStr, this.projName);
                    mkdir(this.projPath);
                    this.turnOnLogFile();
                end
                if ~isempty(this.hSubFS)
                    this.resPath = sprintf("%s/%s", this.projPath, this.resDirName);
                    mkdir(this.resPath);
                end
            end 
            
        end
        
        function updateFileSystem(this, update)
            this.scanIdentifier = sprintf(this.scanIdentifierSuffixModel, update.scanIdentifierVars);
            if ~isempty(this.hSubFS) && ~this.stackAllSubObjRes
                scanPath = sprintf("%s/%s", this.resPath, this.scanIdentifier);
                mkdir(scanPath);
                for i = 1:length(this.hSubFS)
                    this.hSubFS(i).updateExtProjPath(scanPath);
                    this.hSubFS(i).updateFileNameIdentifier(this.scanIdentifier)
                end
            end     
        end
        
        function updateIdentifier(this, scanId)
            this.scanIdentifierPrefix = scanId; %Chiled: Z=1;
            this.scanIdentifierModel  = sprintf("%s-%s", this.scanIdentiferPrefix, this.scanIdentifierSuffixModel); %Child: 
            this.fileNameModel        = sprintf("%s-%s", this.objName, this.scanIdentifierModel);
        end 
        
        function updateExtProjPath(this, path)
            this.extProjPath = path; 
        end
        
        function saveResultsToDisk(this, res)
            if this.saveResults
                fprintf("%s: Saving results.\n", this.fsName);
                if strcmp(this.scanIdentifierPrefix, "")
                    dataFileName = sprintf("%s/%s-Results.mat", this.projPath, this.objName);
                else
                    dataFileName = sprintf("%s/%s-%s-Results.mat", this.projPath, this.objName, this.scanIdentifierPrefix);
                end
                save(dataFileName, '-struct', 'res', '-v7.3');
            end
        end
        
        function saveVarsToDisk(this, vars)
            if this.saveVars
                fprintf("%s: Saving Variables.\n", this.fsName);                 
                if this.useExtVarsPath
                    varsFileName = sprintf("%s/%s%s-Vars.mat", this.extVarsPath, this.ObjName, this.extVarsPrefix);
                else
                    if this.extProject
                        varsFileName = sprintf("%s/%s-%s-Vars.mat", this.projPath, this.objName, this.scanIdentifierPrefix);
                    else
                        if strcmp(this.defaultScanIdentifierPrefix, "")
                            varsFileName = sprintf("%s/%s-Vars.mat",  this.projPath, this.objName);
                        else
                            varsFileName = sprintf("%s/%s-%s-Vars.mat",  this.projPath, this.objName, this.defaultScanIdentifierPrefix);
                        end
                    end
                end
                save(varsFileName, '-struct', 'vars', '-v7.3');
            end
        end
        
        function enableSaveVars(this)
           this.saveVars = true;
        end
        
        function disableSaveVars(this)
           this.saveVars = false;
        end
        
        function turnOnLogFile(this)
           if this.saveAny && ~this.extProject
               diary(sprintf("%s/log.txt", this.projPath));
               diary on
           end
        end
        
        function turnOffLogFile(this)
            if this.saveAnyTot && ~this.extProject
                diary(sprintf("%s/log.txt", this.projPath))
                diary off
            end
        end
        
        function closeFileSystem(this)
            if this.saveAnyTot && ~this.extProject()
                this.turnOffLogFile(); 
            end
        end
        
        function saveAntTotVal = calcSaveAnyTot(this)
           this.saveAnySub = false;
           if ~isempty(this.hSubFS)
               for i=1:length(this.hSubFS)
                this.saveAnySub = this.saveAnySub || getSaveAny(this.hSubFS(i));
               end
           end
           this.saveAnyTot = this.saveAny || this.saveAnySub;
           saveAntTotVal = this.saveAnyTot;
        end
            
    end
end

