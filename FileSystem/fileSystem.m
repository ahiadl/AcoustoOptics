classdef fileSystem < handle
    
    properties
        %New Used Variables
        uVars
        hSubFS      % handles to sub filesystems
        hOwner
        
        fsName   % for logging prints
        objName  % for saving results 
        
        dirPath     % Directory in the disk in which this project directory is created
        projPath    % This project directory (full path) [resDir/projName]
        extProjPath % Project path given bu external FS object
%         resPath     % path to the subObjResults directory   [resDir/projName/resDirName]
        
        projName    % project Name, relevant in case of no owner FS.
        resDirName  % The folder under projPath in which the results should be saved
        
        saveFigs    % should figures be saved?
        saveResults % should results be saved?
        saveVars    % should variables be saved? (if user asked to ssave results save vars as well unless specifically requested not to.
        dontSaveVars% this is configured by the user
        
        saveAny           % is there any action of saving in this object?
        saveAnySub        % is there any action of saving in subFS?
        saveAnyTot        % is there any action of saving in total?
        extProject        % is this an external project?
        stackAllSubObjRes % should all subObj Results not be divided into folders (for direct owner of acostoOptics)?
        
        useExtVarsPath    % shold Vars be saved in this projPath or in user given path?
        extVarsPath       % external path for Vars
        extVarsName       % name of the vars file
        
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
            vars.dontSaveVars = false;
            
            vars.dirPath     = [];
            vars.projName    = [];
            vars.resDirName  = [];
            
            vars.extProject        = false;
            vars.stackAllSubObjRes = false;
            vars.dontSaveVars      = false;
            
            vars.useExtVarsPath    = false;
            vars.extVarsPath       = [];
        end
    end
    
    methods
        function this = fileSystem(hOwner, subObjHandle)
            this.hOwner = hOwner;
            if nargin>1
                this.hSubFS = subObjHandle;
            else
                this.hSubFS = [];
            end
            this.projPath = [];
            this.resDirName = "Results";
            this.stackAllSubObjRes = false;
            this.defaultScanIdentifierPrefix = "";
            this.useExtVarsPath = false;
        end 
        
        function setUserVars(this, uVars)
            this.uVars = uVars;
            
            this.saveResults  = uVars.saveResults;
            this.saveFigs     = uVars.saveFigs;
            this.dontSaveVars = uVars.dontSaveVars;
            this.saveVars     = (this.saveResults || this.saveFigs) && (~this.dontSaveVars); 
            
            this.saveAny = this.saveVars || this.saveResults || this.saveFigs;
            
            this.projName = uVars.projName;
            this.dirPath  = uVars.dirPath;
            
            this.extProject            = uVars.extProject;
            if ~uVars.extProject
                this.useExtVarsPath = false;
%                 this.projName = uVars.scanName;
            end
            this.stackAllSubObjRes     = uVars.stackAllSubObjRes; %not relevant for ao
        end
        
        function configFileSystem(this)  
            if this.calcSaveAnyTot() 
                % If it is an external project, the project path should be
                % given already at this point using the updateFileSystem 
                % method of the owner FS.
                if ~this.extProject
                    this.scanIdentifierPrefix = this.defaultScanIdentifierPrefix;
                    dateStr       = strrep(datestr(datetime('now')), ':', '-');
                    this.projPath = sprintf("%s/%s-%s", this.dirPath, dateStr, this.projName);
                    mkdir(this.projPath);
%                     this.turnOnLogFile();
                end
            end 
        end
        
        function updateFileSystem(this, updateVars)
            this.scanIdentifier = sprintf(this.scanIdentifierSuffixModel, updateVars);
            if ~isempty(this.hSubFS)&& this.saveAnyTot
                if ~this.stackAllSubObjRes
                    scanPath = sprintf("%s/%s/%s", this.projPath, this.resDirName, this.scanIdentifier);
                    mkdir(scanPath);
                else
                    scanPath = sprintf("%s/%s", this.projPath, this.resDirName);
                    if ~exist(scanPath, 'dir')
                        mkdir(scanPath)
                    end
                end
                for i = 1:length(this.hSubFS)
                    this.hSubFS(i).updateExtProjPath(scanPath);
                    this.hSubFS(i).updateIdentifier(this.scanIdentifier)
                end
            end     
        end
        
        function updateIdentifier(this, scanId)
            this.scanIdentifierPrefix = scanId; %Chiled: Z=1;
            this.scanIdentifierModel  = sprintf("%s-%s", this.scanIdentifierPrefix, this.scanIdentifierSuffixModel); %Child: 
            this.fileNameModel        = sprintf("%s-%s", this.objName, this.scanIdentifierModel);
        end 
        
        function setExtVars(this, path, name)
            this.useExtVarsPath = true;
            this.extVarsPath    = path;
            this.extVarsName    = name;
            
            if ~isempty(this.hSubFS)
               for i=1:length(this.hSubFS)
                   this.hSubFS(i).setExtVars(path, name)
               end
            end
        end
        
        function updateExtProjPath(this, path)
            this.projPath = path; 
        end
        
        function saveResultsToDisk(this, res)
            if this.saveResults
                fprintf("%s: Saving results.\n", this.fsName);
                if strcmp(this.scanIdentifierPrefix, "")
                    dataFileName = sprintf("%s/%s-Results.mat", this.projPath, this.objName);
                else
                    dataFileName = sprintf("%s/%s-%s-Results.mat", this.projPath, this.objName, this.scanIdentifierPrefix);
                end
                if exist(dataFileName, 'file')
                    save(dataFileName, '-struct', 'res', '-append', '-v7.3');
                else
                    if length(res) > 1
                        save(dataFileName, 'res', '-v7.3');
                    else
                        save(dataFileName, '-struct', 'res', '-v7.3');
                    end
                end
            end
        end
        
        function saveVarsToDisk(this)
            if this.saveVars
                fprintf("%s: Saving Variables.\n", this.fsName);                 
                if this.useExtVarsPath
                    varsFileName = sprintf("%s/%s%s-Vars.mat", this.extVarsPath, this.objName, this.extVarsName);
                else
                    if this.extProject
                        varsFileName = sprintf("%s/%s-%s-Vars.mat", this.projPath, this.objName, this.scanIdentifierPrefix);
                    else
                        if strcmp(this.scanIdentifierPrefix, "")
                            varsFileName = sprintf("%s/%s-Vars.mat",  this.projPath, this.objName);
                        else
                            varsFileName = sprintf("%s/%s-%s-Vars.mat",  this.projPath, this.objName, this.scanIdentifierPrefix);
                        end
                    end
                end
                vars = this.hOwner.getVars();
                save(varsFileName, '-Struct', 'vars', '-v7.3');
            end
        end
        
        function saveSubFSVarsToDisk(this)
            if ~isempty(this.hSubFS)
                for i=1:length(this.hSubFS)
                    this.hSubFS(i).enableSaveVars();
                    if ~this.useExtVarsPath
                        this.hSubFS(i).setExtVars(this.projPath, "");
                    end
                    this.hSubFS(i).saveVarsToDisk();
                    this.hSubFS(i).disableSaveVars();
                end
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
        
        function saveAnyTotVal = calcSaveAnyTot(this)
           this.saveAnySub = false;
           if ~isempty(this.hSubFS)
               for i=1:length(this.hSubFS)
                this.saveAnySub = this.saveAnySub || this.hSubFS(i).calcSaveAnyTot();
               end
           end
           this.saveAnyTot = this.saveAny || this.saveAnySub;
           saveAnyTotVal = this.saveAnyTot;
        end
            
    end
end

