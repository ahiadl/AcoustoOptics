classdef scan3DAO < handle
    % SCAN3D This object performs 2D scan and extract 3D acousto optics data. 
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
        grid
        s2DVars
        figsVars
        fileSystemVars
        generalVars
        
        curPos      %[third]              [mm]
        curPosIdx   %[thirdIdx][#]
        owned
        extStages
    end
    
    methods (Static)
        function uVars = uVarsCreate()
            uVars.ao         = acoustoOptics.createUserVars();
            uVars.s2D        = scan2DAO.createUserVars();
            uVars.figs       = scan3DGraphics.createUserVars();
            uVars.fileSystem = fileSystemS3D.uVarsCreate();
           
            uVars.grid.repeats = 1;
            
            uVars.grid.scan1Start  = 0;
            uVars.grid.scan1Stride = 0;
            uVars.grid.scan1End    = 0;

            uVars.grid.scan2Start  = 0;
            uVars.grid.scan2Stride = 0;
            uVars.grid.scan2End    = 0;
            
            uVars.grid.thirdPos = 0;
            
            uVars.grid.scan1Axis  = 'Y';
            uVars.grid.scan2Axis  = 'Z';
            uVars.grid.thirdAxis  = 'X';
            uVars.grid.depthAxis  = 'X';
            
%             uVars.stages.scan1AxStageId = 1;
%             uVars.stages.scan2AxStageId = 2;
            uVars.stages.stagesOwnedByParent = false;
            uVars.stages.dontMoveNonScaning  = true;
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
                this.s2D = handles.s2D;
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
            
            %--------------------------------------------------------------
            fprintf("S3D: 1. Set Stages variables and Axes.\n");
            % Only in case of internal stages (e.g. running from script and
            % not from GUI).
            % This describes in general which stage belongs to what axis.
            this.generalVars.dontMoveNonScaning = uVars.stages.dontMoveNonScaning;
            this.generalVars.stagesOwnedByParent = uVars.stages.stagesOwnedByParent;
            if ~this.generalVars.stagesOwnedByParent
                this.stages.assignStagesAxes([uVars.grid.scan1Axis, uVars.grid.scan2Axis]);
            end
           
            %--------------------------------------------------------------
            fprintf("S3D: 2. Setting FileSystem variables.\n");
            this.fileSystem.setUserVars(uVars.fileSystem);
            
            %--------------------------------------------------------------
            fprintf("S3D: 3. Forming scan 2D variables.\n");
            s2DuVars = scan2DAO.createUserVars();
            
            s2DuVars.ao = uVars.ao;
            
            %General
            s2DuVars.grid.repeats = uVars.grid.repeats;
            s2DuVars.general.externalScan = true;
            
            %Grid
            s2DuVars.grid.scanStart  = uVars.grid.scan1Start;
            s2DuVars.grid.scanStride = uVars.grid.scan1Stride;
            s2DuVars.grid.scanEnd    = uVars.grid.scan1End;
            
            s2DuVars.grid.secondPos = uVars.grid.scan2Start;
            s2DuVars.grid.thirdPos  = uVars.grid.thirdPos;
            
            s2DuVars.grid.scanAxis   = uVars.grid.scan1Axis;
            s2DuVars.grid.secondAxis = uVars.grid.scan2Axis;
            s2DuVars.grid.thirdAxis  = uVars.grid.thirdAxis;
            s2DuVars.grid.depthAxis  = uVars.grid.depthAxis;
           
            %Stages
            s2DuVars.stages.moveNonScanStages   = false;
            s2DuVars.stages.stagesOwnedByParent = true;
            
            %Graphics
            s2DuVars.figs = uVars.s2D.figs;
            
            s2DuVars.figs.depthIdx = uVars.figs.depthIdx;
            
            s2DuVars.figs.intExt             = uVars.figs.intExt;
            s2DuVars.figs.hOwnerGraphics     = this.graphics;
            s2DuVars.figs.validOwnerGraphics = true;
            s2DuVars.figs.useExtClims        = uVars.figs.normColorsToPlane;
            
            s2DuVars.figs.scanAxType    = uVars.figs.scan1AxType;
            s2DuVars.figs.depthAxType   = uVars.figs.depthAxType;
            s2DuVars.figs.reopenFigures = uVars.figs.reopenFigures;
            
            s2DuVars.figs.fonts         = uVars.figs.fonts;
            
            %Filesystem
            s2DuVars.fileSystem.extProject   = true;
            s2DuVars.fileSystem.dontSaveVars = true;
            s2DuVars.fileSystem.saveResults  = uVars.fileSystem.saveResults;
            
            this.s2D.setUserVars(s2DuVars)
            this.s2DVars = this.s2D.getVars();
                        
            %--------------------------------------------------------------
            fprintf("S3D: 4. Setting Grid variables.\n");
            this.grid.repeats     = uVars.grid.repeats;
            
            this.grid.scan1Start   = uVars.grid.scan1Start;
            this.grid.scan1Stride  = uVars.grid.scan1Stride;
            this.grid.scan1End     = uVars.grid.scan1End;
            
            this.grid.scan2Start  = uVars.grid.scan2Start;
            this.grid.scan2Stride = uVars.grid.scan2Stride;
            this.grid.scan2End    = uVars.grid.scan2End;
            
            this.grid.thirdPos    = uVars.grid.thirdPos;
            
            this.grid.scan1Axis   = uVars.grid.scan1Axis;
            this.grid.scan2Axis   = uVars.grid.scan2Axis;
            this.grid.thirdAxis   = uVars.grid.thirdAxis;
            this.grid.depthAxis   = uVars.grid.depthAxis;
            
            this.calc3DGrid();
            
            %--------------------------------------------------------------
            fprintf("S3D: 5. Setting strings models and variables\n");
            strings{1} = "start3DScan";
            models{1}  = sprintf("Start Scan for (%s) = (%s)", this.grid.scan2Axis, "%.2f");
            strings{2} = "timeTable3D";
            models{2}  = sprintf("%s%s", this.grid.scan2Axis, "%.2f");

            this.sf.setStringModels(strings, models);
            
            %--------------------------------------------------------------
            fprintf("S3D: 6. Setting figures variables.\n");
            figs = this.graphics.createOwnerVars();
            
            % Misc
            figs.normColorsToPlane = uVars.figs.normColorsToPlane;
            figs.reopenFigures     = uVars.figs.reopenFigures;
            
            % Vectors
            figs.repeats    = uVars.grid.repeats;
            
            figs.scan1AxType    = uVars.figs.scan1AxType;
            figs.scan1AxisNorm  = this.grid.scan1Vec;
            figs.scan2AxType    = uVars.figs.scan2AxType;
            figs.scan2AxisNorm  = this.grid.scan2Vec;
            figs.depthAxType    = uVars.figs.depthAxType;
            figs.depthAxisNorm  = this.s2DVars.ao.measVars.algo.len.depthVec;

            % Labels
            figs.scan1AxLabel = uVars.grid.scan1Axis;
            figs.scan2AxLabel = uVars.grid.scan2Axis;
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
            vars.s2D            = this.s2DVars;
            vars.grid           = this.grid;
            vars.general        = this.generalVars;
            vars.fileSystem     = this.fileSystemVars;
            vars.uVars          = this.uVars;
            vars.uVars.figs     = scan3DGraphics.createUserVars();
            vars.uVars.s2D.figs = scan2DGraphics.createOwnerVars();
            vars.uVars.ao.figs  = AOGraphics.createOwnerVars();
        end
        
        % Results Functions
        
        function initResultsArrays3D(this)      
            scan1AxLen    = this.grid.scan1IdxLen;
            scan2AxLen    = this.grid.scan2IdxLen;
            depthLen      = this.grid.depthIdxLen;
            
            repeats    = this.grid.repeats;
            
            this.res.phi        = zeros(scan1AxLen, scan2AxLen, depthLen, repeats);
            this.res.phiAvg     = zeros(scan1AxLen, scan2AxLen, depthLen);
            this.res.phiAvgStd  = zeros(scan1AxLen, scan2AxLen, depthLen);
            
            this.curPosIdx = zeros(1,3);
        end

        function put2DResTo3DResultsArray(this, res2D)
            this.res.phi(        :, this.curPosIdx, :, :)       = res2D.phi;
            this.res.phiAvg(     :, this.curPosIdx, :)          = res2D.phiAvg;
            this.res.phiAvgStd(  :, this.curPosIdx, :)          = res2D.phiAvgStd;
        end
        
        % Config Functions
        function configureScan(this)
            this.timeTable.scan = struct();
            this.initResultsArrays3D();

            this.graphics.setGraphicsScanVars();
            this.graphics.updateGraphicsConstruction();
            
            this.curPos    = this.grid.scan2Start;
            this.curPosIdx = 1;
            
            this.fileSystem.configFileSystem(this.grid.scan2Axis);
            this.fileSystem.saveVarsToDisk();
            this.s2D.configureScan();
        end
        
        % Scan Functions
        function res = scan3D(this)
            if ~this.generalVars.dontMoveNonScaning
                this.stages.moveAbsAx(this.grid.thirdAxis, this.grid.thirdPos)
            end
            for i=1:this.grid.scan2IdxLen
                this.updateCurPosAndIdx(i);
                this.sf.startScanTime("timeTable3D", 'singlePlane', this.curPos);
                this.sf.printStrModel("start3DScan", this.curPos, true);
                
                this.sf.startScanTime("timeTable3D", 'prepare', this.curPosIdx);
                this.fileSystem.updateFileSystem(this.curPos);
                this.s2D.updateThirdAxPos(this.curPos);
                this.sf.stopScanTime("timeTable3D", 'prepare', this.curPosIdx);

                this.stages.moveAbsAx(this.grid.scan2Axis, this.curPos)
                res2D = this.s2D.scan2D();

                this.sf.startScanTime("timeTable3D", 'copyPlane', this.curPos);
                this.put2DResTo3DResultsArray(res2D);
                this.sf.stopScanTime("timeTable3D", 'copyPlane', this.curPos);

                this.sf.stopScanTime("timeTable3D", 'singlePlane', this.curPos);
            end
            this.fileSystem.saveResultsToDisk(this.res);
            this.fileSystem.closeFileSystem();
            res = this.res;
            this.sf.printStr(sprintf("Done 3D scan.\n"), true);
        end

        % Position Functions
        function updateCurPosAndIdx(this, scan2)
            % curScanIdx - always aligned to (R, 1st, 2nd)[#]
            % curScan    - always aligned to ( 1st, 2nd)[mm]
            this.curPosIdx = scan2;
            this.curPos    = this.grid.scan2Vec(scan2);
            
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
            this.grid.scan1Vec     = this.s2DVars.grid.scanVec;
            this.grid.scan1Len     = this.s2DVars.grid.scanLen;
            this.grid.scan1IdxLen  = this.s2DVars.grid.scanIdxLen;
            this.grid.scan1Idx     = this.s2DVars.grid.scanIdx;
            this.grid.scan1Cntr    = this.s2DVars.grid.scanCntr;
            this.grid.scan1Zero    = this.s2DVars.grid.scanZero;
            
            this.grid.depthVec     = this.s2DVars.grid.depthVec;
            this.grid.depthLen     = this.s2DVars.grid.depthLen;
            this.grid.depthIdxLen  = this.s2DVars.grid.depthIdxLen;
            this.grid.depthIdx     = this.s2DVars.grid.depthIdx;
            this.grid.depthCntr    = this.s2DVars.grid.depthCntr;
            this.grid.depthZero    = this.s2DVars.grid.depthZero;
            
            this.grid.scan2Vec    = this.grid.scan2Start : this.grid.scan2Stride : this.grid.scan2End;
            this.grid.scan2Len    = abs(this.grid.scan2Start - this.grid.scan2End);
            this.grid.scan2IdxLen = length(this.grid.scan2Vec);
            this.grid.scan2Idx    = 1:1:this.grid.scan2IdxLen;
            this.grid.scan2Cntr   = this.grid.scan2Vec - min(this.grid.scan2Vec);
            this.grid.scan2Zero   = abs(this.grid.scan2Vec - this.grid.scan2Vec(1));

            switch this.grid.scan1Axis
                case 'X'
                    this.grid.xVec = this.grid.scan1Vec;
                    switch this.grid.scan2Axis
                        case 'Y'
                            this.grid.yVec = this.grid.scan2Vec;
                            this.grid.zVec = this.grid.depthVec;
                            this.grid.depthAxis = 'Z';
                            this.grid.mainPlane = 'XZ';
                        case 'Z'
                            this.grid.zVec = this.grid.scan2Vec;
                            this.grid.yVec = this.grid.depthVec;
                            this.grid.depthAxis = 'Y';
                            this.grid.mainPlane = 'XY';
                    end
                case 'Y'
                    this.grid.yVec = this.grid.scan1Vec;
                    switch this.grid.scan2Axis
                        case 'X'
                            this.grid.xVec = this.grid.scan2Vec;
                            this.grid.zVec = this.grid.depthVec;
                            this.grid.depthAxis = 'Z';
                            this.grid.mainPlane = 'YZ';
                        case 'Z'
                            this.grid.zVec = this.grid.scan2Vec;
                            this.grid.xVec = this.grid.depthVec;
                            this.grid.depthAxis = 'X';
                            this.grid.mainPlane = 'XY';
                    end   
                case 'Z'
                    this.grid.zVec = this.grid.scan1Vec;
                    switch this.grid.scan2Axis
                        case 'X'
                            this.grid.xVec = this.grid.scan2Vec;
                            this.grid.yVec = this.grid.depthVec;
                            this.grid.depthAxis = 'Y';
                            this.grid.mainPlane = 'YZ';
                        case 'Y'
                            this.grid.yVec = this.grid.scan2Vec;
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
        
        function setNavVars(this, ax, idx, rep)
            navVars.navAx  = ax;
            navVars.navIdx = idx;
            navVars.navRep = rep; 
            
            this.graphics.setNavVars(navVars);
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

