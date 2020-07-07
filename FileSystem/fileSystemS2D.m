classdef fileSystemS2D < fileSystem
    
    properties

        
    end
    
    methods (Static)
    
    function uVars = uVarsFileSystemCreate()
            uVars = fileSystem.uVarsCreate();
        end
    end
    
    methods
        function this = fileSystemS2D(hOwner, subObjHandle)
            this@fileSystem(hOwner, subObjHandle)
            this.uVars = fileSystemS2D.uVarsCreate();
            this.objName = "2DScan";
            this.fsName = "S2D";
            this.resDirName = "AOResults";
        end
        
        function setUserVars(this, uVars)
           uVars.stackAllSubObjRes = true;
           setUserVars@fileSystem(this, uVars);
        end
        
        function configFileSystem(this, firstAxis)
           configFileSystem@fileSystem(this);
           this.scanIdentifierSuffixModel =  sprintf("R-%s-%s-%s", "%d", firstAxis, "%.2f");
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

