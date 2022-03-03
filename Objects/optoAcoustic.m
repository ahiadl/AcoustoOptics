classdef optoAcoustic < handle
    % SCAN3D This object perform 2D scan and extract 3D acousto optics data. 
    % This object contain 2 main subobject: scan2D, which contains the AO #
    % and stages. 
    
    properties
        daq
        bp
        stages
        graphics
        fileSystem
        sf
        owner
        
        res
        res3D
        resLive
        timeTable
        
        % Vars
        uVars
        extVars
        measVars
        grid
        
        curPos
        curPosIdx
        
        owned
        runLive
    end
    
    methods (Static)
        function uVars = uVarsCreate()
            uVars.activeChannels  = 1:256; %Mask
            uVars.numOfSamples    = 2030;
            uVars.numOfAvg        = 1;
            uVars.numOfFrames     = 1;
            uVars.delay           = 0;

            uVars.c          = 1480;
            uVars.fs         = 40e6;
            uVars.geometry   = 'Circular';
            uVars.imageWidth = 30e-3; %[m]
            uVars.BPmode     = 3;

            uVars.firstAxis  = 'X';
            uVars.secondAxis = 'Y';
            uVars.scanAxis   = 'Z';
            
            uVars.scanStart  = 190;
            uVars.scanStride = 1;
            uVars.scanEnd    = 200;
            
            uVars.figs         = OAGraphics.createUserVars();
            uVars.fileSystem   = fileSystemOA.uVarsCreate();
            
        end
    end
    
    methods
        function this = optoAcoustic(stagesH)
            if (nargin > 0) && ischar(stagesH)
                this.stages = stages(stagesH);
            elseif (nargin > 0 )
                this.stages = stagesH;
            else
                this.stages = [];
            end
            
            this.daq        = DAQ();
            this.bp         = BackProjection();
            this.graphics   = OAGraphics();
            this.fileSystem = fileSystemOA(this);
            this.sf         = statsFunctions("OA");
 
        end

        function setVars(this, uVars)
            this.uVars = uVars;
            
            % DAQ
            uVarsDAQ = this.daq.uVarsCreate();
            
            uVarsDaq.delay          = 0; %Currently not supported
            uVarsDaq.activeChannels = uVars.activeChannels;
            uVarsDaq.numOfSamples   = uVars.numOfSamples;
            uVarsDaq.numOfAvg       = uVars.numOfAvg;
            uVarsDaq.numOfFrames    = uVars.numOfFrames;
            
            this.daq.setUserVars(uVarsDaq);
            
            this.extVars.daq = uVarsDaq;
            this.measVars.daq = this.daq.getVars();
            
            % Backprojection
            uVarsBP.c          = uVars.c;
            uVarsBP.fs         = uVars.fs;
            uVarsBP.geometry   = uVars.geometry;
            uVarsBP.imageWidth = uVars.imageWidth; %[m]
            uVarsBP.mode       = uVars.BPmode;

            this.extVars.bp = uVarsDaq;
            this.bp.setVars(uVarsBP);
            this.measVars.bp = this.bp.getVars();
            
            % Grid
            this.grid.scanStart  = uVars.scanStart;
            this.grid.scanStride = uVars.scanStride;
            this.grid.scanEnd    = uVars.scanEnd;
            
            this.grid.firstAxis  = uVars.firstAxis;
            this.grid.secondAxis = uVars.secondAxis;
            this.grid.scanAxis   = uVars.scanAxis;
            
            this.calcGrid();
            this.measVars.grid = this.grid;
            
            % Graphics
            uVarsGraphics = this.graphics.createOwnerVars();
            
            uVarsGraphics.intExt         = uVars.figs.intExt;
            uVarsGraphics.validStruct    = uVars.figs.validStruct;
            uVarsGraphics.extH           = uVars.figs.extH;
            uVarsGraphics.reopenFigures  = uVars.figs.reopenFigures;
            
            uVarsGraphics.numOfChannels  = length(uVars.activeChannels);
            uVarsGraphics.numOfSamples   = uVars.numOfSamples;
            uVarsGraphics.reconSize      = this.measVars.bp.recon.n;
            uVarsGraphics.imageWidth     = this.measVars.bp.recon.imageWidth;
            
            uVarsGraphics.tVec           = this.measVars.bp.timing.t;
            uVarsGraphics.reconAxis      = this.measVars.bp.recon.reconAxis;
            
            uVarsGraphics.firstAxLabel   = this.grid.firstAxis;
            uVarsGraphics.secondAxLabel  = this.grid.secondAxis;
            uVarsGraphics.scanAxLabel    = this.grid.scanAxis;
            
            uVarsGraphics.firstAxType    = uVars.figs.firstAxType;
            uVarsGraphics.secondAxType   = uVars.figs.secondAxType;
            
            uVarsGraphics.grid = this.grid;
%             uVarsGraphics.displayAxAsIdx = uVars.figs.displayAxAsIdx;
            
            uVarsGraphics.useExtClims   = uVars.figs.useExtClims;
            uVarsGraphics.scanAxCoor    = uVars.figs.scanAxCoor;
            
            this.extVars.figs = uVarsGraphics;
            this.graphics.setUserVars(uVarsGraphics);
            
            % Filesystem
            uVarsFS = this.fileSystem.uVarsCreate();
            
            uVarsFS.saveResults  = uVars.fileSystem.saveResults;
%             uVarsFS.saveFigs     = uVars.fileSystem.saveFigs;
            uVarsFS.dontSaveVars = uVars.fileSystem.dontSaveVars;
            
            uVarsFS.projName     = uVars.fileSystem.projName;
            uVarsFS.dirPath      = uVars.fileSystem.dirPath;
            uVarsFS.resDirName   = uVars.fileSystem.resDirName;
            
            uVarsFS.extProject   = uVars.fileSystem.extProject; 
            
            this.extVars.fs = uVarsFS;
            this.fileSystem.setUserVars(uVarsFS);
            
            fprintf("OA: 5. Setting strings models and variables\n");
                       
            strings{1} = "start3DScan";
            models{1}  = sprintf("Start 3D OA Scan");
            strings{2} = "start2DScan";
            models{2}  = sprintf("Start Scan for (%s) = (%s)", this.grid.scanAxis, "%.2f");
            strings{3} = "timeTable3D";
            models{3}  = sprintf("%s", this.grid.scanAxis);
            strings{4} = "timeTable2D";
            models{4}  = sprintf("%s%s", this.grid.scanAxis, "%.2f");
            this.sf.setStringModels(strings, models);
            
        end
        
        function vars = getVars(this)
            vars.measVars = this.measVars;
            vars.extVars  = this.extVars;
            vars.uVars    = this.uVars;
        end
        
        function configPeripherals(this)
           % config fileSystem
           this.fileSystem.configFileSystem();
           this.daq.configAndConnect();
           this.graphics.setGraphicsDynamicVars();
           this.initResultsArray();
        end
        
        function configScan(this)
           this.configPeripherals();
           this.initResultsArray();
        end
        
        function res = acqAndBP(this)
            this.res.sigMat = this.daq.acquire();
%             tmp = load("./sigMatSample.mat");
%             this.res.sigMat = tmp.A;
            this.res.recon = this.bp.calcBP(this.res.sigMat);
            
            res = this.res;
        end
        
        function res = runOA(this)
            this.fileSystem.saveVarsToDisk();
            res = acqAndBP(this);
            this.plotAll();
            this.fileSystem.saveResultsToDisk(this.res);
        end
        
        function res = runLiveOA(this)
%             this.fileSystem.saveVarsToDisk();
            this.runLive = true;
            i=1;
            while (this.runLive)
                this.resLive{i} = acqAndBP(this);
                this.plotAll();
                i=i+1;
            end
            this.fileSystem.saveResultsToDisk(this.res);
            res = this.resLive;
            
        end
        
        function updateCurPosAndIdx(this, idx)
           this.curPosIdx = idx;
           this.curPos = this.grid.scanVec(idx);
        end
        
        function initResultsArray(this)
            this.res.recon = zeros(this.grid.firstIdxLen, this.grid.secondIdxLen);
            this.res.sigMat = zeros( this.measVars.daq.numOfSamples, this.measVars.daq.numOfActiveChannels);
            
            this.res3D.recon = zeros(this.grid.firstIdxLen, this.grid.secondIdxLen, this.grid.scanIdxLen);
            this.res3D.sigMat = zeros(this.measVars.daq.numOfSamples, this.measVars.daq.numOfActiveChannels, this.grid.scanIdxLen);
        end
        
        function putCurResToScanRes(this, curRes, idx)
            this.res3D.recon(:,:, idx) = curRes.recon;
            this.res3D.sigMat(:,:, idx) = curRes.sigMat;
        end
        
        function res = scan(this)
%             this.graphics.resetGraphics();
            this.updateCurPosAndIdx(1);
%             this.fileSystem.saveVarsToDisk();
            this.sf.startScanTime("timeTable3D", 'CompleteScan', []);
            this.sf.printStrModel("start3DScan", [], true);
            for i=1:this.grid.scanIdxLen
                this.updateCurPosAndIdx(i);
                this.sf.startScanTime("timeTable2D", 'singleImage', [this.curPos]);
                this.sf.printStrModel("start2DScan", [this.curPos], true);

                this.sf.startScanTime("timeTable2D", 'prepare', [this.curPos]);
%                 this.fileSystem.updateFileSystem([this.curPos(1)]);
                this.sf.stopScanTime("timeTable2D", 'prepare', [this.curPos]);

                this.stages.moveAbsAx(this.grid.scanAxis, this.curPos(1))

                this.sf.startScanTime("timeTable2D", 'acqAndBP', [this.curPos]);
                curRes = this.acqAndBP();
                this.putCurResToScanRes(curRes, i);
                this.sf.stopScanTime("timeTable2D", 'acqAndBP', [this.curPos]);

                this.sf.startScanTime("timeTable2D", 'scanPlotting', [this.curPos]);
                this.plotScan();
                this.sf.stopScanTime("timeTable2D", 'scanPlotting', [this.curPos]);

                this.sf.stopScanTime("timeTable2D", 'singleImage', [this.curPos]);
            end
            
            this.sf.stopScanTime("timeTable3D", 'CompleteScan', []);

%             this.fileSystem.saveResultsToDisk(this.res3D);
%             this.fileSystem.closeFileSystem();
            res = this.res3D;
            
            this.sf.printStr(sprintf("Done 3D scan.\n"), true);
        end
        
        function plotAll(this)
%             this.graphics.updateCurPos(this.curPos, this.curPosIdx);
            this.graphics.setData(this.res);
            this.graphics.dispSinogram();
            this.graphics.dispRecon();
        end
        
        function plotScan(this)
            this.graphics.updateCurPos(this.curPos, this.curPosIdx);
            this.graphics.setDataScan(this.res3D);
            this.graphics.dispSinogram();
            this.graphics.dispRecon();
        end
        
        function calcGrid(this)
            this.grid.scanVec     = this.grid.scanStart  : this.grid.scanStride  : this.grid.scanEnd;
            this.grid.scanLen     = abs(this.grid.scanStart - this.grid.scanEnd);
            this.grid.scanIdxLen  = length(this.grid.scanVec);
            this.grid.scanIdx     = 1:1:this.grid.scanIdxLen;        
            this.grid.scanCntr    = this.grid.scanVec - mean(this.grid.scanVec);
            this.grid.scanZero    = abs(this.grid.scanVec - this.grid.scanVec(1));
            this.grid.scanNormal  = this.grid.scanVec;
            
            this.grid.firstVec     = this.measVars.bp.recon.reconAxis;
            this.grid.dFirst       = this.measVars.bp.recon.dxdy;
            this.grid.firstLen     = this.measVars.bp.recon.reconAxis(end)+this.grid.dFirst;
            this.grid.firstIdxLen  = length(this.grid.firstVec);
            this.grid.firstIdx     = 1:1:this.grid.firstIdxLen;
            this.grid.firstCntr    = this.grid.firstVec - mean(this.grid.firstVec);
            this.grid.firstZero    = abs(this.grid.firstVec - this.grid.firstVec(1));
            this.grid.firstNormal  = this.grid.firstVec - min(this.grid.firstVec);
            
            this.grid.secondVec     = this.measVars.bp.recon.reconAxis;
            this.grid.dSecond       = this.measVars.bp.recon.dxdy;
            this.grid.secondLen     = this.measVars.bp.recon.reconAxis(end)+this.grid.dSecond;
            this.grid.secondIdxLen  = length(this.grid.secondVec);
            this.grid.secondIdx     = 1:1:this.grid.secondIdxLen;
            this.grid.secondCntr    = this.grid.secondVec - mean(this.grid.secondVec);
            this.grid.secondZero    = abs(this.grid.secondVec - this.grid.secondVec(1));
            this.grid.secondNormal  = this.grid.secondVec - min(this.grid.secondVec);
        end
        
    end
end

