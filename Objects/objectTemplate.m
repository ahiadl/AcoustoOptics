classdef scan3DAO < handle
    % SCAN3D This object perform 2D scan and extract 3D acousto optics data. 
    % This object contain 2 main subobject: scan2D, which contains the AO #
    % and stages. 
    
    properties
        s2D
        stages
        graphics
        fileSystem
        sf
        owner
        
        res
        timeTable
        
        % Vars
        uVars
        s2DVars
        figsVars
        fileSystemVars
        
        owned
        extStages
    end
    
    methods (Static)
        function uVars = uVarsCreate()
            uVars.ao         = acoustoOptics.uVarsCreate();
            uVars.s2D        = scan2DAO.uVarsCreate();
            uVars.figs       = templateGraphics.createUserVars();
            uVars.fileSystem = templateFS.uVarsCreate();
           
            % add here the vars needed for operation 
            
        end
    end
    
    methods
        function this = scan3DAO(handles)
            % handles should have 2 fields:
            % a. s2D - handle to initialized s2D object with initialized
            % stages and AO
            % b. stages - handle to the stages object
            % stages can be given solely, but if s2D handle is given its
            % stages handle must be given as well
            this.sf       = statsFunctions("S3D");
            this.graphics = scan3DGraphics();
            
            if (nargin > 0) && isfield(handles, 'stages')
                this.stages = handles.stages;
                this.extStages = true;
            else
                this.stages = stages('COM3');
                this.extStages = false; 
            end
            
            if (nargin > 0) && isfield(handles, 's2D')
                this.s2D = handles.s2d;
            else
                this.s2D = scan2DAO([], this.stages);
            end 
            
            % this must be here so AO instrreset will not ruin stages
            % connection
            if ~this.extStages
                fprintf("S3D: Creating a stages object...\n");
                this.stages.connect();
                fprintf("S3D: Done connecting to stages.\n");
            end
            hS3DFS = this.s2D.fileSystem;
            this.fileSystem = fileSystemS3D(this, hS3DFS);
        end

        % Manage Variables
        function setUserVars(this, uVars)
            fprintf("S3D: Setting Variables to 2D-Scan.\n");
            this.uVars = uVars();
            
            fprintf("S3D: 4. Set Stages variables and Axes.\n");
            % Only in case of internal stages (e.g. running from script and
            % not from GUI).
            % This describes in general which stage belongs to what axis.
            if ~this.extStages
                this.stages.assignStagesAxes([uVars.grid.firstAxis,        uVars.grid.secondAxis],...
                                             [uVars.stages.firstAxStageId, uVars.stages.secondAxStageId]);
            end
            
            fprintf("S3D: 2. Setting FileSystem variables.\n");
            this.fileSystem.setUserVars(uVars.fileSystem);

            fprintf("S3D: 1. Forming scan 2D variables.\n");
            s2DuVars = scan2DAO.uVarsCreate();
            
            s2DuVars.ao = uVars.ao;
            
            %General
            s2DuVars.grid.repeats = uVars.grid.repeats;
            s2DuVars.general.externalScan = true;
            
            %Grid
            s2DuVars.grid.firstStart  = uVars.grid.firstStart;
            s2DuVars.grid.firstStride = uVars.grid.firstStride;
            s2DuVars.grid.firstEnd    = uVars.grid.firstEnd;

            s2DuVars.grid.secondPos  = uVars.grid.secondStart;
            
            s2DuVars.grid.firstAxis  = uVars.grid.firstAxis;
            s2DuVars.grid.secondAxis = uVars.grid.secondAxis;
            
            %Graphics
            s2DuVars.figs = uVars.s2D.figs;

            s2DuVars.figs.intExt             = uVars.figs.intExt;
            s2DuVars.figs.hOwnerGraphics     = this.graphics;
            s2DuVars.figs.validOwnerGraphics = true;
            s2DuVars.figs.useExtClims        = uVars.figs.normColorsToPlane;
            
            s2DuVars.figs.fonts              = uVars.figs.fonts;
            
            %Filesystem
            s2DuVars.fileSystem.extProject   = true;
            s2DuVars.fileSystem.dontSaveVars = true;
            s2DuVars.fileSystem.saveResults  = uVars.fileSystem.saveResults;
            
            this.s2D.setUserVars(s2DuVars)
            this.s2DVars = this.s2D.getVars();
            
            fprintf("S3D: 3. Setting Grid variables.\n");
            this.grid.repeats      = uVars.grid.repeats;
            
            this.grid.firstStart   = this.s2DVars.grid.firstStart;
            this.grid.firstStride  = this.s2DVars.grid.firstStride;
            this.grid.firstEnd     = this.s2DVars.grid.firstEnd;
            
            this.grid.secondStart  = uVars.grid.secondStart;
            this.grid.secondStride = uVars.grid.secondStride;
            this.grid.secondEnd    = uVars.grid.secondEnd;
            
            this.grid.firstAxis  = uVars.grid.firstAxis;
            this.grid.secondAxis = uVars.grid.secondAxis;

            this.calc3DGrid();

            fprintf("S3D: 5. Setting strings models and variables\n");
            strings{1} = "start3DScan";
            models{1}  = sprintf("Start Scan for (%s) = (%s)", this.grid.secondAxis, "%.2f");
            strings{2} = "timeTable3D";
            models{2}  = sprintf("%s%s", this.grid.secondAxis, "%.2f");

            this.sf.setStringModels(strings, models);
            
            fprintf("S3D: 6. Setting figures variables.\n");
            figs = this.graphics.createOwnerVars();
            
            % Misc
            figs.normColorsToPlane = uVars.figs.normColorsToPlane;
            
            % Vectors
            figs.repeats    = uVars.grid.repeats;
            figs.firstAxis  = this.grid.firstVec;
            figs.secondAxis = this.grid.secondVec;
            figs.depthAxis  = this.grid.depthVec;
            
            % Labels
            figs.firstAxLabel = uVars.grid.firstAxis;
            figs.secondAxLabe = uVars.grid.secondAxis;
            figs.depthAxLabel = this.grid.depthAxis;
            figs.mainPlane    = this.grid.mainPlane;
            
            % Handles
            figs.intExt      = uVars.figs.intExt;
            figs.validStruct = uVars.figs.validStruct;
            figs.extH        = uVars.figs.extH;

            this.graphics.setUserVars(figs);
            this.figsVars = figs;
        end
        
        function vars = getVars(this)
            vars.s2D           = this.s2DVars;
            vars.grid          = this.grid;
            vars.general       = this.generalVars;
            vars.fileSystem    = this.fileSystemVars;
            vars.uVars         = this.uVars;
            vars.uVars.figs    = scan3DGraphics.createUserVars();
            vars.uVars.s2D.figs = scan2DGraphics.createOwnerVars();
            vars.uVars.ao.figs = AOGraphics.createOwnerVars();
        end
        
        % Results Functions
        
        function initResultsArrays3D(this)      
            firstAxLen     = this.grid.firstIdxLen;
            secondAxLen    = this.grid.secondIdxLen;
            depthLen       = this.grid.depthIdxLen;
            
            repeats    = this.grid.repeats;
            numOfQuant = this.s2DVars.ao.measVars.algo.samples.numOfQuant;
            channels   = this.s2DVars.ao.measVars.algo.digitizer.channels;
            
            this.res.phiChCmplx = zeros(firstAxLen, secondAxLen, depthLen, repeats, numOfQuant, channels);
            this.res.phiCh      = zeros(firstAxLen, secondAxLen, depthLen, repeats, numOfQuant, channels);
            this.res.phiQuant   = zeros(firstAxLen, secondAxLen, depthLen, repeats, numOfQuant);
            this.res.phiStd     = zeros(firstAxLen, secondAxLen, depthLen, repeats);
            this.res.phi        = zeros(firstAxLen, secondAxLen, depthLen, repeats);
            this.res.phiAvg     = zeros(firstAxLen, secondAxLen, depthLen);
            this.res.phiAvgStd  = zeros(firstAxLen, secondAxLen, depthLen);
            
            this.curPosIdx = zeros(1,3);
        end

        function putAOResTo3DResultsArray(this, res2D)
            this.res.phiChCmplx( :, this.curPosIdx, :, :, :, :) = res2D.phiChCmplx;
            this.res.phiCh(      :, this.curPosIdx, :, :, :, :) = res2D.phiCh;
            this.res.phiQuant(   :, this.curPosIdx, :, :, :)    = res2D.phiQuant;
            this.res.phi(        :, this.curPosIdx, :, :)       = res2D.phi;
            this.res.phiStd(     :, this.curPosIdx, :, :)       = res2D.phiStd;
            this.res.phiAvg(     :, this.curPosIdx, :)          = res2D.phi;
            this.res.phiAvgStd(  :, this.curPosIdx, :)          = res2D.phiStd;
        end
        
        % Config Functions
        function configureScan(this)
            this.timeTable.scan = struct();
            this.initResultsArrays3D();

            this.graphics.setGraphicsScanVars();
            this.graphics.updateGraphicsConstruction();
            
            this.curPos    = this.grid.secondStart;
            this.curPosIdx = 1;
            
            this.fileSystem.configFileSystem(this.grid.secondAxis);
            this.fileSystem.saveVarsToDisk();
            this.s2D.configureScan();
        end
        
        % Scan Functions
        function res = scan3D(this)
            for i=1:this.grid.secondIdxLen
                this.updateCurPosAndIdx(i);
                this.sf.startScanTime("timeTable3D", 'singlePlane', this.curPos);
                this.sf.printStrModel("start3DScan", this.curPos, true);
                
                this.sf.startScanTime("timeTable3D", 'prepare', this.curPosIdx);
                this.fileSystem.updateFileSystem(this.curPos);
                this.s2D.updateSecondAxPos(this.curPos);
                this.sf.stopScanTime("timeTable3D", 'prepare', this.curPosIdx);

                this.stages.moveStageAxisAbs(this.grid.secondAxis, this.curPos)
                res2D = this.s2D.scan2D();

                this.sf.startScanTime("timeTable3D", 'copyPlane', this.curPos);
                this.putAOResTo3DResultsArray(res2D);
                this.sf.stopScanTime("timeTable3D", 'copyPlane', this.curPos);

                this.sf.stopScanTime("timeTable3D", 'singlePlane', this.curPos);
            end
            this.fileSystem.saveResultsToDisk(this.res);
            this.fileSystem.closeFileSystem();
            res = this.res;
            this.sf.printStr(sprintf("Done 3D scan.\n"), true);
        end

        % Position Functions
        function updateCurPosAndIdx(this, second)
            % curScanIdx - always aligned to (R, 1st, 2nd)[#]
            % curScan    - always aligned to ( 1st, 2nd)[mm]
            this.curPosIdx = second;
            this.curPos    = this.grid.secondVec(second);
            
            this.graphics.updateS3DCurPosAndIdx(this.curPosIdx, this.curPos);
        end
        
        function pns = getPos(this)
            s2DPos = this.s2D.getPos();
            pns.curPos(1) = s2DPos.curPos(1);
            pns.curPos(2) = this.curPos;
            pns.curPosIdx(1) = s2DPos.curPosIdx(1);
            pns.curPosIdx(2) = s2DPos.curPosIdx(2);
            pns.curPosIdx(3) = this.curPosIdx;
        end
        
        function calc3DGrid(this)          
            this.grid.firstVec     = this.s2DVars.grid.firstVec;
            this.grid.firstLen     = this.s2DVars.grid.firstLen;
            this.grid.firstIdxLen  = this.s2DVars.grid.firstIdxLen;
            this.grid.firstIdx     = this.s2DVars.grid.firstIdx;
            this.grid.firstNorm    = this.s2DVars.grid.firstNorm;
            
            this.grid.depthVec     = this.s2DVars.grid.depthVec;
            this.grid.depthLen     = this.s2DVars.grid.depthLen;
            this.grid.depthIdxLen  = this.s2DVars.grid.depthIdxLen;
            this.grid.depthIdx     = this.s2DVars.grid.depthIdx;
            this.grid.depthNorm    = this.s2DVars.grid.depthNorm;
            
            this.grid.secondVec    = this.grid.secondStart : this.grid.secondStride : this.grid.secondEnd;
            this.grid.secondLen    = abs(this.grid.secondStart - this.grid.secondEnd);
            this.grid.secondIdxLen = length(this.grid.secondVec);
            this.grid.secondIdx    = 1:1:this.grid.secondIdxLen;
            this.grid.secondNorm   = this.grid.secondVec - min(this.grid.secondVec);

            switch this.grid.firstAxis
                case 'X'
                    this.grid.xVec = this.grid.firstVec;
                    switch this.grid.secondAxis
                        case 'Y'
                            this.grid.yVec = this.grid.secondVec;
                            this.grid.zVec = this.grid.depthVec;
                            this.grid.depthAxis = 'Z';
                            this.grid.mainPlane = 'XZ';
                        case 'Z'
                            this.grid.zVec = this.grid.secondVec;
                            this.grid.yVec = this.grid.depthVec;
                            this.grid.depthAxis = 'Y';
                            this.grid.mainPlane = 'XY';
                    end
                case 'Y'
                    this.grid.yVec = this.grid.firstVec;
                    switch this.grid.secondAxis
                        case 'X'
                            this.grid.xVec = this.grid.secondVec;
                            this.grid.zVec = this.grid.depthVec;
                            this.grid.depthAxis = 'Z';
                            this.grid.mainPlane = 'YZ';
                        case 'Z'
                            this.grid.zVec = this.grid.secondVec;
                            this.grid.xVec = this.grid.depthVec;
                            this.grid.depthAxis = 'X';
                            this.grid.mainPlane = 'XY';
                    end   
                case 'Z'
                    this.grid.zVec = this.grid.firstVec;
                    switch this.grid.secondAxis
                        case 'X'
                            this.grid.xVec = this.grid.secondVec;
                            this.grid.yVec = this.grid.depthVec;
                            this.grid.depthAxis = 'Y';
                            this.grid.mainPlane = 'YZ';
                        case 'Y'
                            this.grid.yVec = this.grid.secondVec;
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

