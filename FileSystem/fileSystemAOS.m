classdef fileSystemAOS < fileSystem
    
    properties

        
    end
    
    methods (Static)
    
    function uVars = uVarsFileSystemCreate()
            uVars = fileSystem.uVarsCreate();
        end
    end
    
    methods
        function this = fileSystemAOS(hOwner, subObjHandle)
            this@fileSystem(hOwner, subObjHandle)
            this.uVars = fileSystemAOS.uVarsCreate();
            this.objName = "AOScan";
            this.fsName  = "AOS";
            this.resDirName = "AOResults";
        end
        
        function setUserVars(this, uVars)
           uVars.stackAllSubObjRes = true;
           setUserVars@fileSystem(this, uVars);
        end
        
        function configFileSystem(this, scan1Label, scan2Label)
           configFileSystem@fileSystem(this);
           this.scanIdentifierSuffixModel =  sprintf("%s-%s-%s-%s", scan2Label, "%.2f", scan1Label, "%.2f");
        end
        
        function saveVarsToDisk(this)
            saveVarsToDisk@fileSystem(this);
            if this.saveVars
                this.saveSubFSVarsToDisk();
            end
        end

        function saveAOVars(this, aoVars)
            this.hSubFS.enableSaveVars();
            if ~this.useExtVarsPath
                this.hSubFS.setExtVars(this.projPath, "");
            end
            this.hSubFS.saveVarsToDisk(aoVars);
            this.hSubFS.disableSaveVars();
        end
    end
end

