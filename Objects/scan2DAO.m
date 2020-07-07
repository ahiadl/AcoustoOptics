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
        function uVars = uVarsCreate()
            uVars.ao         = acoustoOptics.uVarsCreate();
            uVars.figs       = scan2DGraphics.createUserVars();
            uVars.fileSystem = fileSystemS2D.uVarsFileSystemCreate();

            uVars.grid.repeats   = 1;
            
            uVars.grid.firstStart     = 0;
            uVars.grid.firstEnd       = 0;
            uVars.grid.firstStride    = 0;

            uVars.grid.secondPos = 0;
            
            uVars.grid.firstAxis  = 'Y';
            uVars.grid.secondAxis = 'X';
            
            uVars.stages.firstAxStageId  = 1;
            uVars.stages.secondAxStageId = 2;
            uVars.stages.moveSecondStage = true;
            uVars.stages.stagesOwnedByParent   = false;
            
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
                this.stages = stages('COM3');
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
            
            fprintf("S2D: 1. Setting general variables.\n"); 
            this.generalVars.externalScan = uVars.general.externalScan;
            this.generalVars.moveSecondStage = uVars.stages.moveSecondStage;
            
            fprintf("S2D: 2. Setting S2D FileSystem variables.\n");
            this.fileSystem.setUserVars(uVars.fileSystem);
            
            fprintf("S2D: 3. Setting AcoustoOptics variables.\n");
            uVars.ao.fileSystem.extProject   = true;
            uVars.ao.fileSystem.saveResults  = uVars.fileSystem.saveResults;
            uVars.ao.fileSystem.dontSaveVars = true;
            
            uVars.ao.figs.intExt   = uVars.figs.intExt;
            uVars.ao.reopenFigures = uVars.figs.reopenFigures;
            
            this.ao.setMeasVars(uVars.ao);
            this.aoVars = this.ao.getVars();
            
            fprintf("S2D: 4. Setting Grid variables.\n");
            this.grid.repeats  = uVars.grid.repeats;
            
            this.grid.firstStart   = uVars.grid.firstStart;
            this.grid.firstStride  = uVars.grid.firstStride;
            this.grid.firstEnd     = uVars.grid.firstEnd;
            
            this.grid.secondPos    = uVars.grid.secondPos;

            this.grid.firstAxis  = uVars.grid.firstAxis;
            this.grid.secondAxis = uVars.grid.secondAxis;

            this.calc2DGrid();

            fprintf("S2D: 5. Set Stages variables and Axes.\n");
            % only in case stages are owned by a parent object dont perform
            % calibration
            this.generalVars.stagesOwnedByParent = uVars.stages.stagesOwnedByParent;
            if ~this.generalVars.stagesOwnedByParent
                this.stages.assignStagesAxes([this.grid.firstAxis, this.grid.secondAxis],...
                                             [this.uVars.stages.firstAxStageId, this.uVars.stages.secondAxStageId]);
            end

            fprintf("S2D: 5. Setting strings models and variables\n");
                       
            strings{1} = "start2DScan";
            models{1}  = sprintf("Start Scan for (%s,%s) = (%s, %s)", this.grid.firstAxis, this.grid.secondAxis, "%.2f", "%.2f");
            strings{2} = "start1DScan";
            models{2}  = sprintf("Start Scan for (R,%s,%s) = (%s, %s, %s)", this.grid.firstAxis, this.grid.secondAxis, "%d", "%.2f", "%.2f");
            strings{3} = "timeTable2D";
            models{3}  = sprintf("%s%s%s%s", this.grid.firstAxis, "%.2f", this.grid.secondAxis, "%.2f");
            strings{4} = "timeTable1D";
            models{4}  = sprintf("R%s%s%s%s%s", "%d", this.grid.firstAxis, "%.2f", this.grid.secondAxis, "%.2f");
            this.sf.setStringModels(strings, models);
            
            fprintf("S2D: 6. Setting figures variables.\n");
            figs = this.graphics.createOwnerVars();

            figs.depthIdx      = uVars.figs.depthIdx;
            figs.useQuant      = uVars.ao.ao.useQuant;
            figs.repeats       = this.grid.repeats;
            figs.useExtClims   = uVars.figs.useExtClims;
            figs.reopenFigures = uVars.figs.reopenFigures;
            
            if this.generalVars.externalScan
                figs.hOwnerGraphics     = uVars.figs.hOwnerGraphics;
                figs.validOwnerGraphics = uVars.figs.validOwnerGraphics;
            else
                figs.validOwnerGraphics = false;
                figs.hOwnerGraphics = [];
            end
            
            figs.firstAxType     = uVars.figs.firstAxType;
            figs.firstAxisNorm = this.grid.firstVec;
            
            figs.depthAxType     = uVars.figs.depthAxType;
            figs.depthAxisNorm = this.aoVars.measVars.algo.len.zVecUSRes;
            
            figs.secondPos     = this.grid.secondPos;
                        
            figs.firstAxLabel = uVars.grid.firstAxis;
            figs.secondAxLabe = uVars.grid.secondAxis;
            figs.depthAxLabel = this.grid.depthAxis;
            figs.mainPlane    = this.grid.mainPlane;
            
            figs.intExt      = uVars.figs.intExt;
            figs.validStruct = uVars.figs.validStruct;
            figs.extH        = uVars.figs.extH;
            
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
            this.initResultsArrays2D();
            
            this.timeTable.scan = struct(); 
            this.graphics.setGraphicsScanVars();
            
            
            this.curPos = [this.grid.firstStart, this.grid.secondPos];
            this.curPosIdx = [1, 1, 1];
            
            %Config subObjects
            this.fileSystem.configFileSystem(this.grid.firstAxis);
            this.fileSystem.saveVarsToDisk();
            this.ao.configPeripherals();
        end
        
        % Results Functions
        function initResultsArrays2D(this)      
            resArr = this.get2DResultsArrayModel();
            
            this.res2D.phi        = resArr.phi;
            this.res2D.phiAvg     = resArr.phiAvg;
            this.res2D.phiAvgStd  = resArr.phiAvgStd;
            
            this.curPosIdx = zeros(1,3);
        end
        
        function resArr = get2DResultsArrayModel(this)
            firstAxLen    = this.grid.firstIdxLen;
            depthLen      = this.aoVars.measVars.algo.samples.numOfPos;
            
            repeats       = this.grid.repeats;
            
            resArr.phi        = zeros(firstAxLen, 1, depthLen, repeats);
            resArr.phiAvg     = zeros(firstAxLen, 1, depthLen);
            resArr.phiAvgStd  = zeros(firstAxLen, 1, depthLen);
        end
        
        function putAOResTo1DResultsArray(this, res)
            this.res1D.phi = permute(gather(res.phi),[1, 3, 2, 4]);
        end
        
        function putAOResTo2DResultsArray(this)
            this.res2D.phi(this.curPosIdx(2), 1, :, this.curPosIdx(1)) = this.res1D.phi;
        end
        
        function averageRepeats1D(this)
            this.res1D.phiAvg    = mean(this.res2D.phi(this.curPosIdx(2), 1, :, 1:this.curPosIdx(1)), 4);
            this.res1D.phiAvgStd = std(this.res2D.phi(this.curPosIdx(2), 1, :, 1:this.curPosIdx(1)), 0, 4);
            
            this.res2D.phiAvg(this.curPosIdx(2), 1, :)    = this.res1D.phiAvg;
            this.res2D.phiAvgStd(this.curPosIdx(2), 1, :) = this.res1D.phiAvgStd;
        end
   
        function res = scan2D(this)
            this.graphics.resetGraphics();
            if this.generalVars.moveSecondStage
                this.stages.moveStageAxisAbs(this.grid.secondAxis, this.curPos(2))
            end
            for r = 1:this.grid.repeats
                this.updateCurPosAndIdx( r, 1);
                this.sf.startScanTime("timeTable2D", 'singleRep', [this.curPosIdx(1), this.curPos(2)]);
                this.sf.printStrModel("start2DScan", [this.curPosIdx(1), this.curPos(2)], true);

                for i=1:this.grid.firstIdxLen
                    this.updateCurPosAndIdx( r, i);
                    this.sf.startScanTime("timeTable1D", 'singlePos', [this.curPosIdx(1), this.curPos]);
                    this.sf.printStrModel("start1DScan", [this.curPosIdx(1), this.curPos], true);
                    
                    this.sf.startScanTime("timeTable1D", 'prepare', [this.curPosIdx(1), this.curPos]);
                    this.fileSystem.updateFileSystem([r, this.curPos(1)]);
                    this.sf.stopScanTime("timeTable1D", 'prepare', [this.curPosIdx(1), this.curPos]);
                    
                    this.stages.moveStageAxisAbs(this.grid.firstAxis, this.curPos(1))

                    this.sf.startScanTime("timeTable1D", 'netAcoustoOptics', [this.curPosIdx(1), this.curPos]);
                    res = this.ao.runAcoustoOptics();
                    this.sf.stopScanTime("timeTable1D", 'netAcoustoOptics', [this.curPosIdx(1), this.curPos]);

                    this.sf.startScanTime("timeTable1D", 'scan1DProcessingTime', [this.curPosIdx(1), this.curPos]);
                    this.putAOResTo1DResultsArray(res);
                    this.putAOResTo2DResultsArray();
                    this.averageRepeats1D();
                    this.sf.stopScanTime("timeTable1D", 'scan1DProcessingTime', [this.curPosIdx(1), this.curPos]);
                    
                    this.sf.startScanTime("timeTable1D", 'scanPlotting', [this.curPosIdx(1), this.curPos]);
                    this.plotResults();
                    this.sf.stopScanTime("timeTable1D", 'scanPlotting', [this.curPosIdx(1), this.curPos]);
                    
                    this.sf.stopScanTime("timeTable1D", 'singlePos', [this.curPosIdx(1), this.curPos]);
                end
                this.sf.stopScanTime("timeTable2D", 'singleRep', [this.curPosIdx(1), this.curPos(2)]);
            end
            
            this.fileSystem.saveResultsToDisk(this.res2D);
            this.fileSystem.closeFileSystem();
            res = this.res2D;
            
            this.sf.printStr(sprintf("Done 2D scan.\n"), true);
        end
        
        % Position Functions
        function updateCurPosAndIdx(this, r, first)
            % curPosIdx - always aligned to (R, 1st, 2nd)[#]
            % curPos    - always aligned to ( 1st, 2nd)[mm]
            this.curPosIdx(1) = r;
            this.curPosIdx(2) = first;
            this.curPos(1)    = this.grid.firstVec(first);
            if r==1 && first ==1
                this.curPos(2) = this.grid.secondPos;
            end
            
            this.graphics.updateCurPosAndIdx(this.getPos());
        end

        function pns = getPos(this)
            pns.curPos    = this.curPos;
            pns.curPosIdx = this.curPosIdx;
        end
        
        function updateSecondAxPos(this, secondPos)
            this.curPos(2)      = secondPos;
            this.grid.secondPos = secondPos;
            switch this.grid.firstAxis
                case 'X'
                    switch this.grid.secondAxis
                        case 'Y'
                            this.grid.yPos = this.grid.secondPos;
                        case 'Z'
                            this.grid.zPos = this.grid.secondPos;
                    end
                case 'Y'
                    switch this.grid.secondAxis
                        case 'X'
                            this.grid.xPos = this.grid.secondPos;
                        case 'Z'
                            this.grid.zPos = this.grid.secondPos;
                    end   
                case 'Z'
                    switch this.grid.secondAxis
                        case 'X'
                            this.grid.xPos = this.grid.secondPos;
                        case 'Y'
                            this.grid.yPos = this.grid.secondPos;
                    end
            end
        end
        
        function calc2DGrid(this)
            this.grid.firstVec     = this.grid.firstStart  : this.grid.firstStride  : this.grid.firstEnd;
            this.grid.firstLen     = abs(this.grid.firstStart - this.grid.firstEnd);
            this.grid.firstIdxLen  = length(this.grid.firstVec);
            this.grid.firstIdx     = 1:1:this.grid.firstIdxLen;        
            this.grid.firstCntr    = this.grid.firstVec - mean(this.grid.firstVec);
            this.grid.firstZero    = abs(this.grid.firstVec - this.grid.firstVec(1));
            
            this.grid.depthVec     = this.aoVars.measVars.algo.len.zVecUSRes;
            dDepth                 = abs(this.grid.depthVec(2) - this.grid.depthVec(1));
            this.grid.depthLen     = this.grid.depthVec(end)+dDepth;
            this.grid.depthIdxLen  = length(this.grid.depthVec);
            this.grid.depthIdx     = 1:1:this.grid.depthIdxLen;
            this.grid.depthCntr    = this.grid.depthVec - mean(this.grid.depthVec);
            this.grid.depthZero    = abs(this.grid.depthVec - this.grid.depthVec(1));
            
            this.grid.secondLen    = 0;
            this.grid.secondIdxLen = 1;
            this.grid.secondIdx    = 1;
            this.grid.secondNorm   = 0;

            switch this.grid.firstAxis
                case 'X'
                    this.grid.xVec = this.grid.firstVec;
                    switch this.grid.secondAxis
                        case 'Y'
                            this.grid.yPos = this.grid.secondPos;
                            this.grid.zVec = this.grid.depthVec;
                            this.grid.depthAxis = 'Z';
                            this.grid.mainPlane = 'XZ';
                        case 'Z'
                            this.grid.zPos = this.grid.secondPos;
                            this.grid.yVec = this.grid.depthVec;
                            this.grid.depthAxis = 'Y';
                            this.grid.mainPlane = 'XY';
                    end
                case 'Y'
                    this.grid.yVec = this.grid.firstVec;
                    switch this.grid.secondAxis
                        case 'X'
                            this.grid.xPos = this.grid.secondPos;
                            this.grid.zVec = this.grid.depthVec;
                            this.grid.depthAxis = 'Z';
                            this.grid.mainPlane = 'YZ';
                        case 'Z'
                            this.grid.zPos = this.grid.secondPos;
                            this.grid.xVec = this.grid.depthVec;
                            this.grid.depthAxis = 'X';
                            this.grid.mainPlane = 'XY';
                    end   
                case 'Z'
                    this.grid.zVec = this.grid.firstVec;
                    switch this.grid.secondAxis
                        case 'X'
                            this.grid.xPos = this.grid.secondPos;
                            this.grid.yVec = this.grid.depthVec;
                            this.grid.depthAxis = 'Y';
                            this.grid.mainPlane = 'YZ';
                        case 'Y'
                            this.grid.yPos = this.grid.secondPos;
                            this.grid.xVec = this.grid.depthVec;
                            this.grid.depthAxis = 'X';
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
            if this.graphics.figs.validStruct.curMainPlane
                this.graphics.dispCurMainPlane()
            end  
        end
        
        function plotResults(this)
            this.graphics.set1DData(this.res1D);
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

