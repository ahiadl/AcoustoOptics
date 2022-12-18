classdef scan2DAO < handle
    %SCAN3D This object is the merging of 3 separate objects: scan 3D and
    %scan 2D, this is for conviniece reasons that take effect in the manner
    %the stages are managed and graphics are managed. Both to of the scan
    %types have separate filesystem and operating mechanisms.
    
    properties
        ao
        stages
        graphics
        fileSystem
        sf
        owner
        
        res2D
        res1D
        timeTable
        
        % Vars
        uVars
        grid
        aoVars
        figsVars
        fileSystemVars
        generalVars
        
        curPos      %[first, second]              [mm]
        curPosIdx   %[r, firstIdx, secondIdx][#]
        owned
        extStages
    end
    
    methods (Static) 
        function uVars = createUserVars()
            uVars.ao         = acoustoOptics.createUserVars();
            uVars.figs       = scan2DGraphics.createUserVars();
            uVars.fileSystem = fileSystemS2D.uVarsFileSystemCreate();

            uVars.grid.repeats   = 1;
            
            uVars.grid.scanStart     = 0;
            uVars.grid.scanEnd       = 0;
            uVars.grid.scanStride    = 0;
            
            uVars.grid.secondPos = 0;
            uVars.grid.thirdPos = 0;
            
            uVars.grid.scanAxis   = 'Y';
            uVars.grid.secondAxis = 'X';
            uVars.grid.thirdAxis  = 'Z';
            uVars.grid.depthAxis  = 'Z';     
            
            uVars.stages.moveNonScanStages   = true;
            uVars.stages.stagesOwnedByParent = false;
            
            uVars.general.externalScan      = true;
        end
    end
    
    methods
        function this = scan2DAO(acoustoOpticHandle, stagesHandle)
            this.sf = statsFunctions("S2D");
            % Acousto optics initialization must be conducted before stages
            % connect since AO includes instrreset, what will damage stages
            % connection if called perliminary
            if nargin>0 && ~isempty(acoustoOpticHandle)
                this.ao = acoustoOpticHandle;
            else
                fprintf("S2D: Creating and connecting to Acousto Optics object...\n");
                this.ao = acoustoOptics();
                this.ao.init();
            end
            
            if nargin>1
                this.stages = stagesHandle;
                this.extStages = true;
            else
                fprintf("S2D: Creating and connecting to a stages object...\n");
                this.stages = stages('Zaber', 'COM3');
                this.extStages = false;
                this.stages.connect();
                fprintf("S2D: Done connecting to stages.\n");
            end
            this.graphics = scan2DGraphics();
            this.fileSystem = fileSystemS2D(this, this.ao.fileSystem);
        end
 
        % Manage Variables
        function setUserVars(this, uVars)
            fprintf("S2D: Setting Variables to 2D-Scan.\n");
            this.uVars = uVars;
            
            %-------------
            % General
            %-------------
            fprintf("S2D: 1. Setting general variables.\n"); 
            this.generalVars.externalScan          = uVars.general.externalScan;
            this.generalVars.moveNonScanStages     = uVars.stages.moveNonScanStages;
            this.generalVars.singleChannelAnalysis = uVars.ao.ao.analyzeSingleCh;
            
            %-------------
            % File System
            %-------------
            fprintf("S2D: 2. Setting S2D FileSystem variables.\n");
            this.fileSystem.setUserVars(uVars.fileSystem);
            
            %---------------
            % Acousto Optics
            %---------------
            fprintf("S2D: 3. Setting AcoustoOptics variables.\n");
            uVars.ao.fileSystem.extProject   = true;
            uVars.ao.fileSystem.saveResults  = uVars.fileSystem.saveResults;
            uVars.ao.fileSystem.dontSaveVars = true;
            
            uVars.ao.figs.intExt   = uVars.figs.intExt;
            uVars.ao.reopenFigures = uVars.figs.reopenFigures;
            
            this.ao.setVars(uVars.ao);
            this.aoVars = this.ao.getVars();            
            
            %---------------
            % Grid
            %---------------
            fprintf("S2D: 4. Setting Grid variables.\n");
            this.grid.repeats  = uVars.grid.repeats;
            
            this.grid.scanStart   = uVars.grid.scanStart;
            this.grid.scanStride  = uVars.grid.scanStride;
            this.grid.scanEnd     = uVars.grid.scanEnd;
            
            this.grid.secondPos   = uVars.grid.secondPos;
            this.grid.thirdPos    = uVars.grid.thirdPos;

            this.grid.scanAxis   = uVars.grid.scanAxis;
            this.grid.depthAxis  = uVars.grid.depthAxis;
            this.grid.secondAxis = uVars.grid.secondAxis;
            this.grid.thirdAxis  = uVars.grid.thirdAxis;
            
            this.calc2DGrid();

            %---------------
            % Stages
            %---------------
            fprintf("S2D: 5. Set Stages variables and Axes.\n");
            % only in case stages are owned by a parent object dont perform
            % calibration
            this.generalVars.stagesOwnedByParent = uVars.stages.stagesOwnedByParent;
            if ~this.generalVars.stagesOwnedByParent
                this.stages.assignStagesAxes([this.grid.scanAxis, this.grid.thirdAxis]);
            end

            fprintf("S2D: 5. Setting strings models and variables\n");
                       
            strings{1} = "start2DScan";
            models{1}  = sprintf("Start Scan for (%s,%s) = (%s, %s)", this.grid.scanAxis, this.grid.thirdAxis, "%.2f", "%.2f");
            strings{2} = "start1DScan";
            models{2}  = sprintf("Start Scan for (R,%s,%s) = (%s, %s, %s)", this.grid.scanAxis, this.grid.thirdAxis, "%d", "%.2f", "%.2f");
            strings{3} = "timeTable2D";
            models{3}  = sprintf("%s%s%s%s", this.grid.scanAxis, "%.2f", this.grid.thirdAxis, "%.2f");
%             models{3}  = sprintf("%s%s%s%s", this.grid.scanAxis, "%.2f");
            strings{4} = "timeTable1D";
            models{4}  = sprintf("R%s%s%s%s%s", "%d", this.grid.scanAxis, "%.2f", this.grid.thirdAxis, "%.2f");
%             models{4}  = sprintf("R%s%s%s", "%d", this.grid.scanAxis, "%.2f");
            this.sf.setStringModels(strings, models);
            
            %---------------
            % Figures
            %---------------
            fprintf("S2D: 6. Setting figures variables.\n");
            figs = this.graphics.createOwnerVars();

            figs.depthIdx      = uVars.figs.depthIdx;
            figs.useFrame      = uVars.ao.ao.useFrame;
            figs.repeats       = this.grid.repeats;
            figs.useExtClims   = uVars.figs.useExtClims;
            figs.reopenFigures = uVars.figs.reopenFigures;
            
            figs.channelsInRes  = this.aoVars.measVars.algo.general.analysisReps;
            figs.singleChannelAnalysis = this.generalVars.singleChannelAnalysis;
            figs.sepChIdx    = uVars.figs.sepChIdx;

            if this.generalVars.externalScan
                figs.hOwnerGraphics     = uVars.figs.hOwnerGraphics;
                figs.validOwnerGraphics = uVars.figs.validOwnerGraphics;
            else
                figs.validOwnerGraphics = false;
                figs.hOwnerGraphics = [];
            end
            
            figs.scanAxType   = uVars.figs.scanAxType;
            figs.scanAxisNorm = this.grid.scanVec;
            
            figs.depthAxType   = uVars.figs.depthAxType;
            figs.depthAxisNorm = this.aoVars.measVars.algo.len.depthVec;
            
            figs.secondPos     = this.grid.secondPos;
            figs.thirdPos      = this.grid.thirdPos;
                        
            figs.scanAxLabel   = uVars.grid.scanAxis;
            figs.secondAxLabel = uVars.grid.secondAxis;
            figs.thirdAxLabel  = uVars.grid.thirdAxis;
            figs.depthAxLabel  = uVars.grid.depthAxis;
            figs.mainPlane     = this.grid.mainPlane;
            
            figs.intExt        = uVars.figs.intExt;
            figs.validStruct   = uVars.figs.validStruct;
            figs.extH          = uVars.figs.extH;
            
            figs.fonts       = uVars.figs.fonts;

            this.graphics.setUserVars(figs);
            this.figsVars = figs;
        end
        
        function vars = getVars(this)
            vars.ao            = this.aoVars;
            vars.grid          = this.grid;
            vars.general       = this.generalVars;
            vars.fileSystem    = this.fileSystemVars;
            vars.uVars         = this.uVars;
            vars.uVars.figs    = scan2DGraphics.createOwnerVars();
            vars.uVars.ao.figs = AOGraphics.createUserVars();
        end
        
        function configureScan(this)
            this.timeTable.scan = struct(); 
            this.graphics.setGraphicsScanVars();
            
            this.curPos = [this.grid.scanStart, this.grid.thirdPos];
            this.curPosIdx = [1, 1, 1];
            
            %Config subObjects
            this.fileSystem.configFileSystem(this.grid.scanAxis);
            this.fileSystem.saveVarsToDisk();
            this.ao.configPeripherals();
        end
        
        % Results Functions
        function initResultsArrays2D(this)      
            resArr = this.get2DResultsArrayModel();
            
            this.res2D.phi        = resArr.phi;
            this.res2D.phiAvg     = resArr.phiAvg;
            this.res2D.phiAvgStd  = resArr.phiAvgStd;
            this.res2D.phiNorm    = resArr.phi;
            this.res2D.phiLog     = resArr.phi;
            
            this.res2D.rawPhi     = resArr.rawPhi;
            
            this.res2D.gMin  = inf;
            this.res2D.gMax  = 0;
            this.res2D.gSpan = 0;
            
            resArr = get1DResultsArrayModel(this);
            
            this.res1D.phi       = resArr.phi;
            this.res1D.phiAvg    = resArr.phi;
            this.res1D.phiAvgStd = resArr.phi;
            
            this.res1D.rawPhi    = resArr.rawPhi;
            
            this.curPosIdx = zeros(1,3);
        end
        
        function resArr = get1DResultsArrayModel(this)
           separateCh  = this.aoVars.measVars.algo.general.analysisReps;
           depthLen    = this.aoVars.measVars.algo.samples.numOfPos;
           depthLenRaw = this.aoVars.measVars.algo.samples.numOfPosAlgo;
           
           resArr.phi = zeros(depthLen, separateCh);
           resArr.rawPhi = zeros(depthLenRaw, separateCh);
        end
        
        function resArr = get2DResultsArrayModel(this)
            scanAxLen         = this.grid.scanIdxLen;
            depthLen          = this.aoVars.measVars.algo.samples.numOfPos;
            depthLenRaw       = this.aoVars.measVars.algo.samples.numOfPosAlgo;
            repeats           = this.grid.repeats;
            separateCh        = this.aoVars.measVars.algo.general.analysisReps;
            
            resArr.phi        = zeros(scanAxLen, 1, depthLen, repeats, separateCh);
            resArr.phiAvg     = zeros(scanAxLen, 1, depthLen, separateCh);
            resArr.phiAvgStd  = zeros(scanAxLen, 1, depthLen, separateCh);
            
            resArr.rawPhi     = zeros(scanAxLen, 1, depthLenRaw, repeats, separateCh);
        end
        
        function putAOResTo1DResultsArray(this, res)
            separateCh        = this.aoVars.measVars.algo.general.analysisReps;
            for i = 1:separateCh
                this.res1D.phi(:,i)    = res(i).phi';
                this.res1D.rawPhi(:,i) = res(i).rawPhi';
            end
        end
        
        function putAOResTo2DResultsArray(this)
            separateCh        = this.aoVars.measVars.algo.general.analysisReps;
            
            for i = 1:separateCh
                this.res2D.phi(   this.curPosIdx(2), 1, :, this.curPosIdx(1), i) = this.res1D.phi(:,i);
                this.res2D.rawPhi(this.curPosIdx(2), 1, :, this.curPosIdx(1), i) = this.res1D.rawPhi(:,i);

                %globals Min/Max/Span - for normalization
                this.res2D.gMin(i)  = min(min(min(min(this.res2D.phi(1:this.curPosIdx(2), :, :, 1: this.curPosIdx(1),i)))));
                this.res2D.gMax(i)  = max(max(max(max(this.res2D.phi(1:this.curPosIdx(2), :, :, 1: this.curPosIdx(1),i)))));
                this.res2D.gSpan(i) = this.res2D.gMax(i) - this.res2D.gMin(i);

                this.res2D.phiNorm(1:this.curPosIdx(2), :, :, 1: this.curPosIdx(1),i) =...
                                                            (this.res2D.phi(1:this.curPosIdx(2), :, :, 1: this.curPosIdx(1), i) - this.res2D.gMin(i))/this.res2D.gSpan(i);
                this.res2D.phiLog(1:this.curPosIdx(2), :, :, 1: this.curPosIdx(1),i) =...
                                                            log(this.res2D.phiNorm(1:this.curPosIdx(2), :, :, 1: this.curPosIdx(1),i));
            end
        end
        
        function averageRepeats1D(this)
            separateCh        = this.aoVars.measVars.algo.general.analysisReps;
            for i=1:separateCh
                this.res1D.phiAvg(:, i)    = mean(this.res2D.phi(this.curPosIdx(2), 1, :, 1:this.curPosIdx(1), i), 4);
                this.res1D.phiAvgStd(:, i) = std(this.res2D.phi(this.curPosIdx(2), 1, :, 1:this.curPosIdx(1), i), 0, 4);

                this.res2D.phiAvg(this.curPosIdx(2), 1, :, i)    = this.res1D.phiAvg(:, i);
                this.res2D.phiAvgStd(this.curPosIdx(2), 1, :, i) = this.res1D.phiAvgStd(:, i);
            end
        end
   
        function res = scan2D(this)
            this.initResultsArrays2D();
            this.graphics.resetGraphics();
            this.graphics.resetPlots();
            this.updateCurPosAndIdx( 1, 1);
            if this.generalVars.moveNonScanStages
                this.stages.moveAbsAx(this.grid.secondAxis, this.curPos(2))
                this.stages.moveAbsAx(this.grid.thirdAxis, this.grid.thirdPos)
            end
            for r = 1:this.grid.repeats
                this.updateCurPosAndIdx( r, 1);
                this.sf.startScanTime("timeTable2D", 'singleRep', [this.curPosIdx(1), this.curPos(2)]);
                this.sf.printStrModel("start2DScan", [this.curPosIdx(1), this.curPos(2)], true);

                for i=1:this.grid.scanIdxLen
                    this.updateCurPosAndIdx( r, i);
                    this.sf.startScanTime("timeTable1D", 'singlePos', [this.curPosIdx(1), this.curPos(1:2)]);
                    this.sf.printStrModel("start1DScan", [this.curPosIdx(1), this.curPos], true);
                    
                    this.sf.startScanTime("timeTable1D", 'prepare', [this.curPosIdx(1), this.curPos(1:2)]);
                    this.fileSystem.updateFileSystem([r, this.curPos(1)]);
                    this.sf.stopScanTime("timeTable1D", 'prepare', [this.curPosIdx(1), this.curPos(1:2)]);
                    
                    this.stages.moveAbsAx(this.grid.scanAxis, this.curPos(1))

                    this.sf.startScanTime("timeTable1D", 'netAcoustoOptics', [this.curPosIdx(1), this.curPos(1:2)]);
                    res = this.ao.runAcoustoOptics();
                    this.sf.stopScanTime("timeTable1D", 'netAcoustoOptics', [this.curPosIdx(1), this.curPos(1:2)]);

                    this.sf.startScanTime("timeTable1D", 'scan1DProcessingTime', [this.curPosIdx(1), this.curPos(1:2)]);
                    this.putAOResTo1DResultsArray(res);
                    this.putAOResTo2DResultsArray();
                    this.averageRepeats1D();
                    this.sf.stopScanTime("timeTable1D", 'scan1DProcessingTime', [this.curPosIdx(1), this.curPos(1:2)]);
                    
                    this.sf.startScanTime("timeTable1D", 'scanPlotting', [this.curPosIdx(1), this.curPos(1:2)]);
                    this.plotResults();
                    this.sf.stopScanTime("timeTable1D", 'scanPlotting', [this.curPosIdx(1), this.curPos(1:2)]);
                    
                    this.sf.stopScanTime("timeTable1D", 'singlePos', [this.curPosIdx(1), this.curPos(1:2)]);
                end
                this.sf.stopScanTime("timeTable2D", 'singleRep', [this.curPosIdx(1), this.curPos(2)]);
            end
            
            this.fileSystem.saveResultsToDisk(this.res2D);
            this.fileSystem.closeFileSystem();
            res = this.res2D;
            
            this.sf.printStr(sprintf("Done 2D scan.\n"), true);
        end
        
        % Position Functions
        function updateCurPosAndIdx(this, r, scanIdx)
            % curPosIdx - always aligned to (R, 1st, 2nd)[#]
            % curPos    - always aligned to ( 1st, 2nd)[mm]
            this.curPosIdx(1) = r;
            this.curPosIdx(2) = scanIdx;
            this.curPosIdx(3) = 1;
            this.curPos(1)    = this.grid.scanVec(scanIdx);
            if r==1 && scanIdx ==1
                this.curPos(2) = this.grid.secondPos;
                this.curPos(3) = this.grid.thirdPos;
            end
            this.graphics.updateCurPosAndIdx(this.getPos());
        end

        function pns = getPos(this)
            pns.curPos    = this.curPos;
            pns.curPosIdx = this.curPosIdx;
        end
        
        function updateThirdAxPos(this, thirdPos)
            this.curPos(3)     = thirdPos;
            this.grid.thirdPos = thirdPos;
        end
        
        function calc2DGrid(this)
            this.grid.scanVec     = this.grid.scanStart  : this.grid.scanStride  : this.grid.scanEnd;
            this.grid.scanLen     = abs(this.grid.scanStart - this.grid.scanEnd);
            this.grid.scanIdxLen  = length(this.grid.scanVec);
            this.grid.scanIdx     = 1:1:this.grid.scanIdxLen;        
            this.grid.scanCntr    = this.grid.scanVec - mean(this.grid.scanVec);
            this.grid.scanZero    = abs(this.grid.scanVec - this.grid.scanVec(1));
            
            this.grid.depthVec     = this.aoVars.measVars.algo.len.depthVec;
            dDepth                 = abs(this.grid.depthVec(2) - this.grid.depthVec(1));
            this.grid.depthLen     = this.grid.depthVec(end)+dDepth;
            this.grid.depthIdxLen  = length(this.grid.depthVec);
            this.grid.depthIdx     = 1:1:this.grid.depthIdxLen;
            this.grid.depthCntr    = this.grid.depthVec - mean(this.grid.depthVec);
            this.grid.depthZero    = abs(this.grid.depthVec - this.grid.depthVec(1));
            
            this.grid.thirdLen    = 0;
            this.grid.thirdIdxLen = 1;
            this.grid.thirdIdx    = 1;
            this.grid.thirdNorm   = 0;
            
            switch this.grid.scanAxis
                case 'X'
                    switch this.grid.thirdAxis
                        case 'Y'
                            this.grid.mainPlane = 'XZ';
                        case 'Z'
                            this.grid.mainPlane = 'XY';
                    end
                case 'Y'
                    switch this.grid.thirdAxis
                        case 'X'
                            this.grid.mainPlane = 'YZ';
                        case 'Z'
                            this.grid.mainPlane = 'XY';
                    end   
                case 'Z'
                    switch this.grid.thirdAxis
                        case 'X'
                            this.grid.mainPlane = 'YZ';
                        case 'Y'
                            this.grid.mainPlane = 'XZ';
                    end
            end
        end

        % Graphics Functions

        function gH = getGraphicsHandle(this)
            gH = this.graphics;
        end
        
        function plotAvgPlots(this)
            if this.uVars.figs.validStruct.curMainAxisAvg
                this.graphics.dispCurMainAxisAvg()
            end
            if this.uVars.figs.validStruct.curMainPlaneAvg
                this.graphics.dispCurMainPlaneAvg()
            end
        end
        
        function plotCurPlots(this)    
            if this.graphics.figs.validStruct.curMainAxis
                this.graphics.dispCurMainAxis()
            end
            if this.graphics.figs.validStruct.curMainAxisLog
                this.graphics.dispCurMainAxisLog()
            end
            if this.graphics.figs.validStruct.curMainPlane
                this.graphics.dispCurMainPlane()
            end  
            if this.graphics.figs.validStruct.curMainPlaneLog
                this.graphics.dispCurMainPlaneLog()
            end
        end
        
        function plotResults(this)
            this.graphics.set1DData(this.res1D);
            this.graphics.setLogData(this.res2D.phiLog); 
            this.graphics.setAvg1DData(this.res1D);
            this.plotCurPlots();
            this.plotAvgPlots();   
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
        
        function setLoaded2DData(this, data, vars)
            setLoadedVars(this, vars.uVars);
            this.res2D = data;
            this.graphics.setGraphicsScanVars();
            this.graphics.updateGraphicsConstruction();
            this.graphics.setLoadedData(data);
        end
    end
end

