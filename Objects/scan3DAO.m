classdef scan3DAO < handle
    %SCAN2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ao
        stages
        graphics
        fileSystem
        sf
        owner
        
        res3D
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
        
%         strings
        curPos      %[first, second]              [mm]
        curPosIdx   %[r, firstIdx, secondIdx][#]
        owned
        extStages
    end
    
    methods (Static)
        function figs = uVarsFiguresCreate()
            figs.zIdx        = 0;
            figs.intExt      = 'int';
            
            figsNames = scan3DGraphics.getGraphicsNames();
            
            for i=1:length(figsNames)
                figs.validStruct.(figsNames{i}) = false;
            end
            
            for i=1:length(figsNames)
                figs.extH.(figsNames{i}) =  Graphics.createHandlesStruct();
            end
            
            figs.fonts.type       = [];
            figs.fonts.titleSize  = 18;
            figs.fonts.labelsSize = 18;
            figs.fonts.axisSize   = 18;
        end
        
        function uVars = uVarsCreate()
            uVars.ao   = acoustoOptics.uVarsCreate();
            uVars.figs = scan3DAO.uVarsFiguresCreate();
            uVars.fileSystem = fileSystemS3D.uVarsCreate();
           
            uVars.grid.startFirst     = 0;
            uVars.grid.startSecond    = 0;
            uVars.grid.endFirst       = 0;
            uVars.grid.endSecond      = 0;
            uVars.grid.strideFirst    = 0;
            uVars.grid.strideSecond   = 0;
            
            uVars.grid.secondPos = 0;
            
            uVars.grid.firstAxis  = 'Y';
            uVars.grid.secondAxis = 'X';
                        
            uVars.general.repeats   = 1;
            uVars.general.useQuant  = true;
            uVars.general.scan2D    = 0;
        end
    end
    
    methods
        function this = scan3DAO(acoustoOpticHandle, stagesHandle)
            this.sf = statsFunctions("S3D");
            
            if nargin>0
                this.ao = acoustoOpticHandle;
            else
                this.ao = acoustoOptics();
            end
            
            if nargin>1
                this.stages = stagesHandle;
                this.extStages = true;
            else
                this.stages = stages();
                this.extStages = false;
                this.stages.connect();
            end
            
            this.graphics = scan3DGraphics();
            this.fileSystem = fileSystemS3D();
        end

        function configure3DScan(this)
            this.fileSystem.turnLogFileOn(); 
            this.timeTable.scan = struct();
            this.initResultsArrays3D();
            this.initResultsArrays2D();
            
            
            this.ao.configPeripherals();
            this.graphics.setGraphicsScanVars();
            this.graphics.updateGraphicsConstruction();
        end
        
        function configure2DScan(this)
            %Config scan3D Object
            this.initResultsArrays2D();
            this.timeTable.scan = struct();
            
            %Config subObjects
            this.fileSystem.configFileSystem();
            this.graphics.setGraphicsScanVars();
            this.graphics.updateGraphicsConstruction();
            
            %Config AO
            this.updateAOFileSystem();
            this.ao.configPeripherals();
        end
        
        function updateAOFileSystem(this)
            if ~uVars.fileSystem.extProject
                uVars.ao.fileSystem.varsFileNameModel  = "AOVars.mat";
                uVars.ao.fileSystem.extProjPath        = this.fileSystemVars.projPath;
                uVars.ao.fileSystem.extProjResultsPath = this.fileSystemVars.rawDataPath;
                uVars.ao.fileSystem.extProjFigsPath    = this.fileSystemVars.figsPath;
                
                filenameAO = sprintf("AO-R%s%s%s%s%s.mat", "%d", this.grid.firstAxis, "%.2f", this.grid.secondAxis,"%.2f");
            else
                
                uVars.ao.fileSystem.varsFileNameModel  = this.fileSystemVars.aoVarsFileNameModel;
                uVars.ao.fileSystem.extProjPath        = this.fileSystemVars.extProjPath;
                uVars.ao.fileSystem.extProjResultsPath = this.fileSystemVars.rawDataPath;
                uVars.ao.fileSystem.extProjFigsPath    = this.fileSystemVars.figsPath;
                
                filenameAO = uVars.fileSystem.fileNameModelAO;
            end
            
            this.ao.fileSystem.setVarsFilenameVariables({});
            this.ao.fileSystem.setDataFilenameModel(filenameAO);
        end
        
        % Manage Variables
        function setUserVars(this, uVars)
            fprintf("S3D: Setting Variables to 2D-Scan.\n");
            this.uVars = uVars();
            
            fprintf("S3D: 1. Setting general variables.\n");
            this.generalVars.useQuant = uVars.general.useQuant;
            this.generalVars.repeats  = uVars.general.repeats;
            this.generalVars.scan2D   = uVars.general.scan2D;
            
            fprintf("S3D: 2. Setting FileSystem variables.\n");
            % Here, just notice wether or not this is an external project;
            
            uVars.fileSystem.saveAny = uVars.fileSystem.saveVars || uVars.fileSystem.saveResults || uVars.fileSystem.saveFigs;
            
            this.fileSystem.setUserVars(uVars.fileSystem);
            this.fileSystemVars = this.fileSystem.configFileSystem();
            
            % Config AO File System
            % this is just to notify AO that this is an external project.
            % The AO object must be updated later with the folders paths 
            % made for this project.
            
            uVars.ao.fileSystem.extProject = true;
            uVars.ao.fileSystem.saveVars   = false;
            
            if ~uVars.fileSystem.extProject
                filename = sprintf("2DScan-(%s-%s).mat", this.grid.secondAxis , "%.2f");
                this.fileSystem.set2DDataFilenameModel(filename)
            end

            fprintf("S3D: 3. Setting AcoustoOptivs variables.\n");
            this.ao.setMeasVars(uVars.ao);
            this.aoVars = this.ao.getAOVars();
            
            fprintf("S3D: 4. Setting Grid variables.\n");
            this.grid.startFirst   = uVars.grid.startFirst;
            this.grid.strideFirst  = uVars.grid.strideFirst;
            this.grid.endFirst     = uVars.grid.endFirst;
            this.grid.startSecond  = uVars.grid.startSecond;
            this.grid.strideSecond = uVars.grid.strideSecond;
            this.grid.endSecond    = uVars.grid.endSecond;
            
            this.grid.firstAxis  = uVars.grid.firstAxis;
            this.grid.secondAxis = uVars.grid.secondAxis;

            this.grid.secondPos = uVars.grid.secondPos;
            
            this.calc3DGrid();

            fprintf("S3D: 5. Set Stages variables and Axes.\n");
            % Only in case of internal stages (e.g. running from script and
            % not from GUI).
            % This describes in general which stage belongs to what axis.
            if ~this.extStages
                if this.generalVars.scan2D
                    this.stages.assignStagesAxes(this.grid.firstAxis, this.uVars.stages.firstAxStageId);
                else
                    this.stages.assignStagesAxes(this.grid.firstAxis, this.uVars.stages.firstAxStageId,...
                                                 this.grid.secondAxis, this.uVars.stages.secondAxStageId);
                end
            end
            
            fprintf("S3D: 5. Setting strings models and variables\n");
                       
            strings{1} = "start3DScan";
            models{1}  = sprintf("Start Scan for (%s) = (%s)", this.grid.secondAxis, "%.2f");
            strings{2} = "start2DScan";
            models{2}  = sprintf("Start Scan for (%s,%s) = (%s, %s)", this.grid.firstAxis, this.grid.secondAxis, "%.2f", "%.2f");
            strings{3} = "start1DScan";
            models{3}  = sprintf("Start Scan for (R,%s,%s) = (%s, %s, %s)", this.grid.firstAxis, this.grid.secondAxis, "%d", "%.2f", "%.2f");
            strings{4} = "timeTable3D";
            models{4}  = sprintf("%s%s", this.grid.secondAxis, "%.2f");
            strings{5} = "timeTable2D";
            models{5}  = sprintf("%s%s%s%s", this.grid.firstAxis, "%.2f", this.grid.secondAxis, "%.2f");
            strings{6} = "timeTable1D";
            models{6}  = sprintf("R%s%s%s%s%s", "%d", this.grid.firstAxis, "%.2f", this.grid.secondAxis, "%.2f");
            this.sf.setStringModels(strings, models);
            
            fprintf("S3D: 6. Setting figures variables.\n");
            figs = this.graphics.createGraphicsUserVars();

            figs.depthIdx = uVars.figs.depthIdx;

            figs.useQuant   = uVars.general.useQuant;
            figs.repeats    = uVars.general.repeats;
            figs.firstAxis  = this.grid.firstVec;
            figs.secondAxis = this.grid.secondVec;
            figs.depthAxis  = this.aoVars.measVars.algo.len.zVecUSRes;

            figs.firstAxLabel = uVars.grid.firstAxis;
            figs.secondAxLabe = uVars.grid.secondAxis;
            figs.depthAxLabel = this.grid.depthAxis;
            figs.mainPlane    = this.grid.mainPlane;
            
            figs.intExt      = uVars.figs.intExt;
            figs.validStruct = uVars.figs.validStruct;
            figs.extH        = uVars.figs.extH;
            figs.fonts       = uVars.figs.fonts;
            
            figs.normColorsToPlane = uVars.figs.normColorsToPlane;
            
            this.graphics.setUserVars(figs);
            this.figsVars = figs;
        end
        
        function vars = getScanVars(this)
            vars.ao            = this.aoVars;
            vars.grid          = this.grid;
            vars.general       = this.generalVars;
            vars.fileSystem    = this.fileSystemVars;
            vars.uVars         = this.uVars;
            vars.uVars.figs    = scan3DGraphics.createGraphicsUserVars();
            vars.uVars.ao.figs = AOGraphics.createGraphicsUserVars();
        end
        
        
        % Results Functions
        
        function initResultsArrays3D(this)      
            firstAxLen     = this.grid.firstIdxLen;
            secondAxLen    = this.grid.secondIdxLen;
            depthLen       = this.aoVars.measVars.algo.samples.numOfPos;
            
            repeats    = this.generalVars.repeats;
            numOfQuant = this.aoVars.measVars.algo.samples.numOfQuant;
            channels   = this.aoVars.measVars.algo.digitizer.channels;
            
            this.res3D.phiChCmplx = zeros(firstAxLen, secondAxLen, depthLen, repeats, numOfQuant, channels);
            this.res3D.phiCh      = zeros(firstAxLen, secondAxLen, depthLen, repeats, numOfQuant, channels);
            this.res3D.phiQuant   = zeros(firstAxLen, secondAxLen, depthLen, repeats, numOfQuant);
            this.res3D.phiStd     = zeros(firstAxLen, secondAxLen, depthLen, repeats);
            this.res3D.phi        = zeros(firstAxLen, secondAxLen, depthLen, repeats);
            this.res3D.phiAvg     = zeros(firstAxLen, secondAxLen, depthLen);
            this.res3D.phiAvgStd  = zeros(firstAxLen, secondAxLen, depthLen);
            
            this.curPosIdx = zeros(1,3);
        end
        
        function initResultsArrays2D(this)      
            resArr = this.get2DResultsArrayModel();
            
            this.res2D.phiChCmplx = resArr.phiChCmplx;
            this.res2D.phiCh      = resArr.phiCh;
            this.res2D.phiQuant   = resArr.phiQuant;
            this.res2D.phiStd     = resArr.phiStd;
            this.res2D.phi        = resArr.phi;
            this.res2D.phiAvg     = resArr.phiAvg;
            this.res2D.phiAvgStd  = resArr.phiAvgStd;
            
            this.curPosIdx = zeros(1,3);
        end
        
        function resArr = get2DResultsArrayModel(this)
            firstAxLen    = this.grid.firstIdxLen;
            depthLen      = this.aoVars.measVars.algo.samples.numOfPos;
            
            repeats    = this.generalVars.repeats;
            numOfQuant = this.aoVars.measVars.algo.samples.numOfQuant;
            channels   = this.aoVars.measVars.algo.digitizer.channels;
            
            resArr.phiChCmplx = zeros(firstAxLen, 1, depthLen, repeats, numOfQuant, channels);
            resArr.phiCh      = zeros(firstAxLen, 1, depthLen, repeats, numOfQuant, channels);
            resArr.phiQuant   = zeros(firstAxLen, 1, depthLen, repeats, numOfQuant);
            resArr.phiStd     = zeros(firstAxLen, 1, depthLen, repeats);
            resArr.phi        = zeros(firstAxLen, 1, depthLen, repeats);
            resArr.phiAvg     = zeros(firstAxLen, 1, depthLen);
            resArr.phiAvgStd  = zeros(firstAxLen, 1, depthLen);
        end
        
        function putAOResTo1DResultsArray(this, res)
            this.res1D.phiChCmplx = permute(gather(res.phiChCmplx), [4, 5, 3, 6, 1, 2]);
            this.res1D.phiCh      = permute(gather(res.phiCh),      [4, 5, 3, 6, 1, 2]);
            this.res1D.phiQuant   = permute(gather(res.phiQuant),   [3, 4, 2, 5, 1]);
            this.res1D.phi        = permute(gather(res.phi),        [1, 3, 2, 4]);
            this.res1D.phiStd     = permute(gather(res.phiStd),     [1, 3, 2, 4]);
        end
        
        function putAOResTo2DResultsArray(this)
            this.res2D.phiChCmplx( this.curPosIdx(2), 1, :, this.curPosIdx(1), :, :) = this.res1D.phiChCmplx;
            this.res2D.phiCh(      this.curPosIdx(2), 1, :, this.curPosIdx(1), :, :) = this.res1D.phiCh;
            this.res2D.phiQuant(   this.curPosIdx(2), 1, :, this.curPosIdx(1), :)    = this.res1D.phiQuant;
            this.res2D.phi(        this.curPosIdx(2), 1, :, this.curPosIdx(1))       = this.res1D.phi;
            this.res2D.phiStd(     this.curPosIdx(2), 1, :, this.curPosIdx(1))       = this.res1D.phiStd;
        end
        
        function putAOResTo3DResultsArray(this)
            this.res3D.phiChCmplx( :, this.curPosIdx(3), :, :, :, :) = this.res2D.phiChCmplx;
            this.res3D.phiCh(      :, this.curPosIdx(3), :, :, :, :) = this.res2D.phiCh;
            this.res3D.phiQuant(   :, this.curPosIdx(3), :, :, :)    = this.res2D.phiQuant;
            this.res3D.phi(        :, this.curPosIdx(3), :, :)       = this.res2D.phi;
            this.res3D.phiStd(     :, this.curPosIdx(3), :, :)       = this.res2D.phiStd;
            this.res3D.phiAvg(     :, this.curPosIdx(3), :)          = this.res2D.phi;
            this.res3D.phiAvgStd(  :, this.curPosIdx(3), :)          = this.res2D.phiStd;
        end

        function averageRepeats1D(this)
            this.res1D.phiAvg    = mean(this.res2D.phi(this.curPosIdx(2), 1, :, 1:this.curPosIdx(1)), 4);
            this.res1D.phiAvgStd = std(this.res2D.phi(this.curPosIdx(2), 1, :, 1:this.curPosIdx(1)), 0, 4);
            
            this.res2D.phiAvg(this.curPosIdx(2), 1, :)    = this.res1D.phiAvg;
            this.res2D.phiAvgStd(this.curPosIdx(2), 1, :) = this.res1D.phiAvgStd;
        end
        
        % Scan Functions
        function res = scan3D(this)
            this.configure3DScan();
            if this.fileSystemVars.saveVars
                this.fileSystem.saveVarsToDisk(this.getScanVars(), "");
                this.ao.saveVarsToDisk("");
            end
            
            for i=1:this.grid.secondIdxLen
                this.updateCurPosAndIdx( [], [], i);
                this.sf.startScanTime("timeTable3D", 'singlePlane', this.curPos(2));
                this.sf.printStrModel("start3DScan", this.curPos, true);

                this.graphics.updateCurPosAndIdx(this.getPos());    
                this.stages.moveStageAxisAbs(this.grid.secondAxis, this.curPos(2))
                this.scan2D();

                this.sf.startScanTime("timeTable3D", 'copyPlane', this.curPos(2));
                this.putAOResTo3DResultsArray();
                this.sf.stopScanTime("timeTable3D", 'copyPlane', this.curPos(2));

                this.sf.stopScanTime("timeTable3D", 'singlePlane', this.curPos(2));
            end

            if this.fileSystemVars.saveResults
                this.fileSystem.save3DResultsToDisk(this.res3D);
            end
            
            res = this.res3D;
            this.sf.printStr(sprintf("Done 3D scan.\n"), true);
            this.fileSystem.turnLogFileOff();
        end
        
        function res = scan2D(this)
            this.fileSystem.set2DDataFilenameVariables({this.curPos(2)});
            for r = 1:this.generalVars.repeats
                this.updateCurPosAndIdx( r, [], []);
                this.sf.startScanTime("timeTable2D", 'singleRep', [this.curPosIdx(1), this.curPos(2)]);
                this.sf.printStrModel("start2DScan", [this.curPosIdx(1), this.curPos(2)], true);

                for i=1:this.grid.firstIdxLen
                    this.updateCurPosAndIdx( [], i, []);
                    this.sf.startScanTime("timeTable1D", 'singlePos', [this.curPosIdx(1), this.curPos]);
                    this.sf.printStrModel("start1DScan", [this.curPosIdx(1), this.curPos], true);
                    
                    this.sf.startScanTime("timeTable1D", 'prepare', [this.curPosIdx(1), this.curPos]);
                    this.graphics.updateCurPosAndIdx(this.getPos()); 
                    this.ao.fileSystem.setDataFilenameVariables({r, this.curPos(1), this.curPos(2)}); 
                    this.sf.stopScanTime("timeTable1D", 'prepare', [this.curPosIdx(1), this.curPos]);
                    
                    this.stages.moveStageAxisAbs(this.grid.firstAxis, this.curPos(1))

                    this.sf.startScanTime("timeTable1D", 'netAcoustoOptics', [this.curPosIdx(1), this.curPos]);
                    res = this.ao.measureAndAnalyse();
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
            
            if this.fileSystemVars.saveResults
                this.fileSystem.save2DResultsToDisk(this.res2D);
            end
            
            res = this.res2D;
            this.sf.printStr(sprintf("Done 1D scan.\n"), true);
        end
        
        % Position Functions
        function updateCurPosAndIdx(this, r, first, second)
            % curScanIdx - always aligned to (R, 1st, 2nd)[#]
            % curScan    - always aligned to ( 1st, 2nd)[mm]
            if ~isempty(r)
                this.curPosIdx(1) = r;
            end
            if~isempty(first)
                this.curPosIdx(2) = first;
                this.curPos(1)    = this.grid.firstVec(first);
            end
            if ~isempty(second)
                this.curPosIdx(3) = second;
                this.curPos(2)    = this.grid.secondVec(second);
            end
        end
        
        function pns = getPos(this)
%             pns.curPos     = this.curPos;
            pns.curPos    = this.curPos;
%             pns.curPosIdx  = this.curPosIdx;
            pns.curPosIdx = this.curPosIdx;
        end
        
        function calc3DGrid(this)
            this.grid.firstVec     = this.grid.startFirst  : this.grid.strideFirst  : this.grid.endFirst;
            this.grid.firstLen     = abs(this.grid.startFirst - this.grid.endFirst);
            this.grid.firstIdxLen  = length(this.grid.firstVec);
            this.grid.firstIdx     = 1:1:this.grid.firstIdxLen;        
            this.grid.firstNorm    = this.grid.firstVec - mean(this.grid.firstVec);
            this.grid.depthVec     = this.aoVars.measVars.algo.len.zVecUSRes;
            
            if this.generalVars.scan2D
                this.grid.secondVec    = this.grid.secondPos;
                this.grid.secondLen    = 0;
                this.grid.secondIdxLen = 1;
                this.grid.secondIdx    = 1;
                this.grid.secondNorm   = 0;
            else
                this.grid.secondVec = this.grid.startSecond : this.grid.strideSecond : this.grid.endSecond;
                this.grid.secondLen    = abs(this.grid.startSecond - this.grid.endSecond);
                this.grid.secondIdxLen = length(this.grid.secondVec);
                this.grid.secondIdx    = 1:1:this.grid.secondIdxLen;
                this.grid.secondNorm   = this.grid.secondVec - min(this.grid.secondVec);
            end
            
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
            if this.graphics.figs.validStruct.nav
                this.graphics.dispNavPlane()
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
        
        function setLoaded3DData(this, data, vars)
            setLoadedVars(this, vars.uVars);
            this.res3D = data;
            this.graphics.setGraphicsScanVars();
            this.graphics.updateGraphicsConstruction();
            this.graphics.setLoadedData(data);
        end
    end
end

