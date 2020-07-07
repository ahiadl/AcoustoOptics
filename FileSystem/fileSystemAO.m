classdef fileSystemAO < fileSystem
    %FILESYSTEMAO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        saveRawData
        saveNetSignal
        saveDemultiplexed
        saveReshapedSignal
        saveFFT
        
        splitMeas
        splitNum
        splitInd
        
        liveAO
        liveAOResDir
        liveAOInd
    end
    
    methods (Static)
        function vars = aoFSVarsCreate()
            vars.saveRawData        = false;
            vars.saveNetSignal      = false;
            vars.saveDemultiplexed  = false;
            vars.saveReshapedSignal = false;
            vars.saveFFT            = false;
        end
        
        function vars = uVarsCreate()
            vars = fileSystem.uVarsCreate();
            
            vars.saveRawData        = false;
            vars.saveNetSignal      = false;
            vars.saveDemultiplexed  = false;
            vars.saveReshapedSignal = false;
            vars.saveFFT            = false;            
        end
    end
    
    methods
        function this = fileSystemAO(hOwner)
            this@fileSystem(hOwner);
            this.fsName  = "AOFS";
            this.objName = "AO";
            this.defaultScanIdentifierPrefix = "";
        end

        function setUserVars(this, uVars)
            uVars.stackAllSubObjRes = false;
            setUserVars@fileSystem(this, uVars)
            
            this.saveRawData        = uVars.saveRawData;
            this.saveNetSignal      = uVars.saveNetSignal;
            this.saveDemultiplexed  = uVars.saveDemultiplexed;
            this.saveReshapedSignal = uVars.saveReshapedSignal;
            this.saveFFT            = uVars.saveFFT;
            
            this.splitMeas          = uVars.splitMeas;
            this.splitNum           = uVars.splitNum;
            this.splitInd           = 1;
            
            saveAlgoData = this.saveRawData        || this.saveNetSignal || this.saveDemultiplexed || ...
                           this.saveReshapedSignal || this.saveFFT;
            this.saveResults = this.saveResults    || saveAlgoData;
            this.saveAny     = this.saveAny        || saveAlgoData;
            this.saveVars    = this.saveAny && ~this.dontSaveVars;

            %Save vars in case there is any saving option on. saving vars
            %can be externally disabled later according to user
            %requirements.
            
            this.liveAO = false; %should be manualy on by liveAO function
        end
        
        function configFileSystem(this)
           configFileSystem@fileSystem(this)
           if ~this.extProject
              this.scanIdentifierPrefix = this.defaultScanIdentifierPrefix;
           end
        end
        
        function saveLiveAOVarsToDisk(this, vars)
           this.enableSaveVars(); 
           this.saveVarsToDisk(vars);
           this.disableSaveVars();
        end
        
        function saveResultsToDisk(this, res)
            if this.saveAnyTot
                if this.saveRawData && ~this.splitMeas
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
                
                data.qAvgChFFT        = res.qAvgChFFT;
                data.unFittedFFT      = res.unFittedFFT;
                data.unFittedFFTShift = res.unFittedFFTShift;
                data.fitModel         = res.fitModel;
                data.fittedFFT        = res.fittedFFT;
                data.phi              = res.phi;

                if this.liveAO
                   if this.saveResults
                        fprintf("%s: Saving Live AO results.\n", this.fsName);
                        if strcmp(this.scanIdentifierPrefix, "")
                            dataFileName = sprintf("%s/%s/%s-Results-%d.mat", this.projPath, this.liveAOResDir, this.objName, this.liveAOInd);
                        else
                            dataFileName = sprintf("%s/%s/%s-%s-Results-%d.mat", this.projPath, this.liveAOResDir, this.objName, this.scanIdentifierPrefix, this.liveAOInd);
                        end
                        save(dataFileName, '-struct', 'data', '-v7.3');
                    end
                else
                   saveResultsToDisk@fileSystem(this, data); 
                end
            end
        end
        
        function saveRawDataToDisk(this, rawData)
            if this.saveResults && this.saveRawData
                fieldName = sprintf("rawData%d", this.splitInd);
                data.(fieldName) = rawData;
                fprintf("%s: Saving rawData.\n", this.fsName);
                if strcmp(this.scanIdentifierPrefix, "")
                    dataFileName = sprintf("%s/%s-Results.mat", this.projPath, this.objName);
                else
                    dataFileName = sprintf("%s/%s-%s-Results.mat", this.projPath, this.objName, this.scanIdentifierPrefix);
                end
                if exist(dataFileName, 'file')
                    save(dataFileName, '-struct', 'data', '-append', '-v7.3');
                else
                    save(dataFileName, '-struct', 'data', '-v7.3');
                end
            end
        end
        
        function initLiveAOFS(this, vars)
            if this.saveAny
                this.liveAO    = true;
                this.liveAOInd = 1;
                this.saveLiveAOVarsToDisk(vars);
                if strcmp(this.scanIdentifierPrefix, "")  
                    this.liveAOResDir = "LiveAOResults";
                else
                    this.liveAOResDir = sprintf("LiveAOResults-%s", this.scanIdentifierPrefix);
                end
                mkdir(sprintf("%s/%s", this.projPath, this.liveAOResDir));
            end
        end
        
        function setLiveAOInd(this, ind)
            this.liveAOInd = ind;
        end
        
        function closeFileSystem(this)
            if ~this.liveAO
                closeFileSystem@fileSystem(this)
            end
        end
        
        function closeLiveAoFileSystem(this)
            if this.liveAO && this.saveAnyTot && ~this.extProject
                this.turnOffLogFile();
            end
        end
        
        function updateSplitInd(this, Ind)
           this.splitInd = Ind; 
        end
        
    end
end

