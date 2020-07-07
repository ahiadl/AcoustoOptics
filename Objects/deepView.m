classdef deepView < handle
    % SCAN3D This object perform 2D scan and extract 3D acousto optics data. 
    % This object contain 2 main subobject: scan2D, which contains the AO #
    % and stages. 
    
    properties
        s2D
        stages      % if required
        graphics
        fileSystem
        sf
        owner
        
        res
        timeTable
        
        % Vars
        uVars
        s2DOwnerVars
        s2DVars
        time
        figsVars
        fileSystemVars
        curPos
        curPosIdx
        
        owned
        extStages
    end
    
    methods (Static)
        function uVars = uVarsCreate()
            uVars.ao         = acoustoOptics.uVarsCreate();
            uVars.s2D        = scan2DAO.uVarsCreate();
            uVars.figs       = dvGraphics.createUserVars();
            uVars.fileSystem = fileSystemDV.uVarsCreate();
           
            % Add here the vars needed for operation 
            
            % Grid
            uVars.grid.repeats = 0;
            
            uVars.grid.firstStart  = 0;
            uVars.grid.firstStride = 0;
            uVars.grid.firstEnd    = 0;

            uVars.grid.secondPos  = 0;
            
            uVars.grid.firstAxis  = 'X';
            uVars.grid.secondAxis = 'Y';
            
            % Stages
            uVars.stages.firstAxStageId  = 2;
            uVars.stages.secondAxStageId = 1;
            
            % Time To Sample
            uVars.time.start  = 2;
            uVars.time.stride = 2;
            uVars.time.end    = 20;
        end
    end
    
    methods
        function this = deepView(hS2D)
            this.sf       = statsFunctions("DV");
            this.graphics = dvGraphics();
            
            if (nargin > 0)
                this.s2D = hS2D;
            else
                this.s2D = scan2DAO();
            end 
            
            hS2DFS = this.s2D.fileSystem;
            this.fileSystem = fileSystemDV(this, hS2DFS);
        end

        % Manage Variables
        function setUserVars(this, uVars)
            fprintf("DV: Setting Variables to 2D-Scan.\n");
            this.uVars = uVars();
            
            fprintf("DV: 1. Calulate Time.\n")
            this.calcTime(uVars.time);
            
            fprintf("DV: 1. Setting FileSystem variables.\n");
            this.fileSystem.setUserVars(uVars.fileSystem);

            fprintf("DV: 2. Forming scan 2D variables.\n");
            s2DuVars = scan2DAO.uVarsCreate();
            
            s2DuVars.ao = uVars.ao;
            s2DuVars.ao.ao.skipParamsCheck = true;
            
            %General
            s2DuVars.grid.repeats = uVars.grid.repeats;
            s2DuVars.general.externalScan = true;
            
            %Grid
            s2DuVars.grid.firstStart  = uVars.grid.firstStart;
            s2DuVars.grid.firstStride = uVars.grid.firstStride;
            s2DuVars.grid.firstEnd    = uVars.grid.firstEnd;

            s2DuVars.grid.secondPos  = uVars.grid.secondPos;
            
            s2DuVars.grid.firstAxis  = uVars.grid.firstAxis;
            s2DuVars.grid.secondAxis = uVars.grid.secondAxis;
            
            %Graphics
            s2DuVars.figs = uVars.s2D.figs;
            s2DuVars.figs.depthIdx           = uVars.figs.depthIdx;
            s2DuVars.figs.firstAxType        = uVars.figs.firstAxType;
            s2DuVars.figs.depthAxType        = uVars.figs.depthAxType;
            s2DuVars.figs.intExt             = uVars.figs.intExt;
            s2DuVars.figs.hOwnerGraphics     = this.graphics;
            s2DuVars.figs.validOwnerGraphics = true;
            s2DuVars.figs.useExtClims        = false;
            s2DuVars.figs.reopenFigures      = uVars.figs.reopenFigures;
            
            s2DuVars.figs.fonts              = uVars.figs.fonts;
            
            %Filesystem
            s2DuVars.fileSystem.extProject   = true;
            s2DuVars.fileSystem.dontSaveVars = false;
            s2DuVars.fileSystem.saveResults  = uVars.fileSystem.saveResults;
            
            this.s2DOwnerVars = s2DuVars;
            this.s2D.setUserVars(s2DuVars)
            this.s2DVars = this.s2D.getVars();

            fprintf("DV: 5. Setting strings models and variables\n");
            strings{1} = "startDeepView";
            models{1}  = "Start Scan for (T) = (%.2f)";
            strings{2} = "timeDeepView";
            models{2}  = "T%.2f";
            this.sf.setStringModels(strings, models);
            
            fprintf("DV: 6. Setting figures variables.\n");
            figs = this.graphics.createOwnerVars();
            
            % Misc
            figs.depthIdx = uVars.figs.depthIdx;
            figs.reopenFigures = uVars.figs.reopenFigures;
            % Vectors
            figs.repeats    = uVars.grid.repeats;
            figs.firstAxis  = this.s2DVars.grid.firstVec;
            figs.depthAxis  = this.s2DVars.grid.depthVec;
            
            figs.timeFrames = this.time.timeVec;
            figs.numOfTimeFrames = this.time.numOfTimeFrames;
            
            % Labels
            figs.firstAxLabel = uVars.grid.firstAxis;
            
            % Handles
            figs.intExt      = uVars.figs.intExt;
            figs.validStruct = uVars.figs.validStruct;
            figs.extH        = uVars.figs.extH;

            this.graphics.setUserVars(figs);
            this.figsVars = figs;
        end
        
        function vars = getVars(this)
            vars.s2D           = this.s2DVars;
            vars.time          = this.time;
            vars.fileSystem    = this.fileSystemVars;
            vars.uVars         = this.uVars;
            vars.uVars.figs    = scan3DGraphics.createUserVars();
            vars.uVars.s2D.figs = scan2DGraphics.createOwnerVars();
            vars.uVars.ao.figs = AOGraphics.createOwnerVars();
        end
        
        % Results Functions
        
        function initResultsArrays3D(this)      
            firstAxLen     = this.s2DVars.grid.firstIdxLen;
            timeAxLen      = this.time.numOfTimeFrames;
            depthLen       = this.s2DVars.grid.depthIdxLen;
            
            repeats    = this.s2DVars.grid.repeats;
%             numOfQuant = this.s2DVars.ao.measVars.algo.samples.numOfQuant;
%             channels   = this.s2DVars.ao.measVars.algo.digitizer.channels;
            
%             this.res.phiChCmplx = zeros(firstAxLen, timeAxLen, depthLen, repeats, numOfQuant, channels);
%             this.res.phiCh      = zeros(firstAxLen, timeAxLen, depthLen, repeats, numOfQuant, channels);
%             this.res.phiQuant   = zeros(firstAxLen, timeAxLen, depthLen, repeats, numOfQuant);
            
            this.res.phiChCmplx = cell(timeAxLen);
            this.res.phiCh      = cell(timeAxLen);
            this.res.phiQuant   = cell(timeAxLen);
            
            this.res.phiStd     = zeros(firstAxLen, timeAxLen, depthLen, repeats);
            this.res.phi        = zeros(firstAxLen, timeAxLen, depthLen, repeats);
            this.res.phiAvg     = zeros(firstAxLen, timeAxLen, depthLen);
            this.res.phiAvgStd  = zeros(firstAxLen, timeAxLen, depthLen);
            
            this.curPosIdx = 0;
        end

        function putAOResToResultsArray(this, res2D)
            this.res.phi(        :, this.curPosIdx, :, :)       = res2D.phi;
            this.res.phiAvg(     :, this.curPosIdx, :)          = res2D.phiAvg;
            this.res.phiAvgStd(  :, this.curPosIdx, :)          = res2D.phiAvgStd;
        end
        
        % Config Functions
        function configureScan(this)
            this.timeTable.scan = struct();
            this.initResultsArrays3D();

            this.graphics.setGraphicsScanVars();

            this.curPos    = this.time.start;
            this.curPosIdx = 1;
            
            this.fileSystem.configFileSystem();
            this.fileSystem.saveVarsToDisk();
%             this.s2D.configureScan();
        end
        
        % Scan Functions
        function res = scanDeepView(this)
            for i=1:this.time.numOfTimeFrames
                this.updateCurPosAndIdx(i);
                this.sf.startScanTime("timeDeepView", 'singlePlane', this.curPos);
                this.sf.printStrModel("startDeepView", this.curPos, true);
                
                this.sf.startScanTime("timeDeepView", 'prepare', this.curPosIdx);
                this.graphics.notifyNewTimeFrame();
                this.fileSystem.updateFileSystem(this.curPos);
                this.sf.stopScanTime("timeDeepView", 'prepare', this.curPosIdx);
                
                this.s2DOwnerVars.ao.ao.timeToSample = this.time.timeVec(i);
                this.s2D.setUserVars(this.s2DOwnerVars);
                this.s2DVars = this.s2D.getVars();
                this.s2D.configureScan();
                this.s2DOwnerVars.ao.ao.skipParamsCheck = false;
                res2D = this.s2D.scan2D();

                this.sf.startScanTime("timeDeepView", 'copyPlane', this.curPos);
                this.putAOResToResultsArray(res2D);
                this.sf.stopScanTime("timeDeepView", 'copyPlane', this.curPos);

                this.sf.stopScanTime("timeDeepView", 'singlePlane', this.curPos);
            end
            this.fileSystem.saveResultsToDisk(this.res);
            this.fileSystem.closeFileSystem();
            res =  this.res;
            this.sf.printStr(sprintf("Done DV scan.\n"), true);
        end

        % Position Functions
        function updateCurPosAndIdx(this, timeIdx)
            % curScanIdx - always aligned to (R, 1st, 2nd)[#]
            % curScan    - always aligned to ( 1st, 2nd)[mm]
            this.curPosIdx = timeIdx;
            this.curPos    = this.time.timeVec(timeIdx);
            
            this.graphics.updateTimeAndInd(this.curPos, this.curPosIdx);
        end
        
        function pns = getPos(this)
            s2DPos = this.s2D.getPos();
            pns.curPos(1) = s2DPos.curPos(1);
            pns.curPos(2) = this.curPos;
            pns.curPosIdx(1) = s2DPos.curPosIdx(1);
            pns.curPosIdx(2) = s2DPos.curPosIdx(2);
            pns.curPosIdx(3) = this.curPosIdx;
        end
        
        function calcTime(this, timeVars)          
            this.time.start  = timeVars.start;
            this.time.stride = timeVars.stride;
            this.time.end    = timeVars.end;
            
            this.time.timeVec = timeVars.start : timeVars.stride : timeVars.end;
            this.time.numOfTimeFrames = length(this.time.timeVec);
            this.time.timeIdxVec = 1:this.time.numOfTimeFrames;
        end

        % Graphics Functions
        function gH = getGraphicsHandle(this)
            gH = this.graphics;
        end
        
        % Loaded data vars
        function setLoadedVars(this, uVars)
            % Deactivate FileSystem
            uVars.fileSystem.saveAny     = false;
            uVars.fileSystem.saveVars    = false;
            uVars.fileSystem.saveResults = false;
            uVars.fileSystem.saveFigs    = false;
            
            uVars.ao.fileSystem.saveRawData        = false;
            uVars.ao.fileSystem.saveNetSignal      = false;
            uVars.ao.fileSystem.saveDemultiplexed  = false;
            uVars.ao.fileSystem.saveReshapedSignal = false;
            uVars.ao.fileSystem.saveFFT            = false;
            uVars.ao.fileSystem.savePhiChCmplx     = false;
            uVars.ao.fileSystem.saveResults        = false;
            uVars.ao.fileSystem.saveFigs           = false;
            uVars.ao.fileSystem.saveVars           = false;
            uVars.ao.fileSystem.saveAny            = false;
            
            this.setUserVars(uVars);
        end
        
        function setLoaded3DData(this, data, vars)
            setLoadedVars(this, vars.uVars);
            this.res3D = data;
            this.graphics.setGraphicsScanVars();
            this.graphics.updateGraphicsConstruction();
            this.graphics.setLoadedData(data);
        end
    end
end

