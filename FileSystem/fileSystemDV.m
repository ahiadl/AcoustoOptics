classdef fileSystemDV < fileSystem
    
    properties
        
    end
    
    methods (Static)
        function uVars = uVarsCreate()
            uVars = fileSystem.uVarsCreate();
        end
    end
    
    methods
        function this = fileSystemDV(hOwner, subObjHandle)
            this@fileSystem(hOwner, subObjHandle)
            this.uVars   = fileSystemS3D.uVarsCreate();
            this.fsName  = "DV";
            this.objName = "DeepView";
            this.resDirName = "TimeResults";
            this.scanIdentifierSuffixModel = "T=%d";
        end
        
        function setUserVars(this, uVars)
           uVars.stackAllSubObjRes = false;
           setUserVars@fileSystem(this, uVars);
        end
        
        function configFileSystem(this)
           configFileSystem@fileSystem(this);
%            this.hSubFS.setExtVars(this.projPath, "");
        end
        
        function saveVarsToDisk(this)
            saveVarsToDisk@fileSystem(this);
%             this.saveSubFSVarsToDisk(); 
        end
    end
end

