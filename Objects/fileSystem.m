classdef fileSystem < handle
    %FILESYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        uVars
        
        dirPath
        resDir
        rawDataDir
        figsDir
        scanName
        
        saveRawDataFlag
        saveFigsFlag
        saveResultsFlag
        
        repFlag
        rep
    end
    
    methods
        function this = fileSystem()
            
        end
        
        function setUserVars(this, uVars)
            this.uVars = uVars;
            saveData = uVars.saveFullData || uVars.saveFigs || uVars.saveResults;
            this.saveRawDataFlag = uVars.saveFullData;
            this.saveFigsFlag    = uVars.saveFigs;
            this.saveResultsFlag = uVars.saveResults;
            this.repFlag    = true;
            this.rep = 1;
            if saveData
                this.scanName = uVars.scanName;
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

        function saveResults(this, res, vars)
            if this.saveResultsFlag
                if this.repFlag
                	resultsName = sprintf("%d-Results.mat", this.rep); 
                    varsName = sprintf("%s/Vars.mat", this.resDir);
                    if this.rep == 1
                        this.saveVars(varsName, vars);
                    end
                    this.rep = this.rep+1;
                else
                    resultsName = sprintf("Results.mat");
                    this.saveVars([this.resDir, '\Vars.mat'], vars)
                end 
                save(sprintf("%s/%s", this.resDir, resultsName), 'res', '-v7.3');    
            end   
        end
        
        function saveVars(this, name, vars)
            save(name, 'vars', '-v7.3');
        end
        
        function saveFigs(this)
            % TODO: complete
        end
        
        function saveRawData(this, name, rawData, vars)
            if this.saveRawData
                varsName    = sprintf("%s-vars.mat", name);
                if this.repFlag
                    rawDataName = sprintf("%d-%s-rawData.mat", this.rep, name); 
                    if this.rep == 1
                        this.saveReducedVars(varsName, vars)
                    end
                    this.rep = this.rep+1;
                    
                else
                    rawDataName = sprintf("%s-rawData.mat", name); 
                    this.saveReducedVars(varsName, vars)
                end
                save(sprintf("%s%s%s", this.rawDataDir,'\', rawDataName),  'rawData', '-v7.3');                
            end
        end

        function saveReducedVars(this, name, vars)          
            vars.acoustoOptics.algoVars.timing.tSigVec = [];
            vars.acoustoOptics.algoVars.timing.tMeasVec = [];
            vars.acoustoOptics.algoVars.timing.tAcqVec = [];
            vars.acoustoOptics.algoVars.timing.tQuantVec = [];
            vars.acoustoOptics.algoVars.timing.tPosVec = [];
            save(sprintf("%s%s%s",this.rawDataDir,'\', name), 'vars', '-v7.3');
        end
    end
end

