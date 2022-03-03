classdef fileSystemOA < fileSystem
    %FILESYSTEMAO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)   
        function uVars = uVarsCreate()
            uVars = fileSystem.uVarsCreate();            
        end
    end
    
    methods
        function this = fileSystemOA(hOwner)
            this@fileSystem(hOwner);
            this.objName = "OA";
            this.fsName  = "OAFS";
            this.defaultScanIdentifierPrefix = "";
        end

        function setUserVars(this, uVars)
            uVars.stackAllSubObjRes = false;
            setUserVars@fileSystem(this, uVars);
        end
        
        function configFileSystem(this)
           configFileSystem@fileSystem(this)
           this.scanIdentifierPrefix = this.defaultScanIdentifierPrefix;
        end

    end
end

