classdef scan2D < scanObj
    %SCAN2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods (Static)
        function figs = uVarsFiguresCreate()
            figs.zIdx        = 0;
            figs.intExt      = 'int';
            
            figsNames = scan2DGraphics.getGraphicsNames();
            
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
            uVars.figs = scan2D.uVarsFiguresCreate();
            uVars.fileSystem = fileSystemS2D.uVarsCreate();
            
%             uVars.fileSystem.scanName        = [];
%             uVars.fileSystem.resDirPath      = [];
%             uVars.fileSystem.saveFullData    = false;
%             uVars.fileSystem.saveReducedData = false;
%             uVars.fileSystem.saveFigs        = false;
%             uVars.fileSystem.savePhiChCmplx  = false;
           
            uVars.grid.startX    = 0;
            uVars.grid.startY    = 0;
            uVars.grid.endX      = 0;
            uVars.grid.endY      = 0;
            uVars.grid.strideX   = 0;
            uVars.grid.strideY   = 0;
            
            uVars.grid.firstAxis = 'Y';
            
            uVars.general.repeats   = 1;
            uVars.general.useQuant  = true;
        end
        
    end
    
    methods
        function this = scan2D(acoustoOpticHandle, owner)
            this@scanObj(acoustoOpticHandle, owner);
            
            this.strings.scan = "Done Scan for (R,X,Y) = (%d, %.2f, %.2f)";
            this.strings.timeTable = "R%dX%.2fY%.2f";
            
            this.graphics = scan2DGraphics();
            this.fileSystem = fileSystemS2D();
        end
        
        function initResultsArrays(this)      
            xLen       = this.grid.xIdxLen;
            yLen       = this.grid.yIdxLen;
            zLen       = this.aoVars.measVars.algo.samples.numOfPos;
            repeats    = this.generalVars.repeats;
            numOfQuant = this.aoVars.measVars.algo.samples.numOfQuant;
            channels   = this.aoVars.measVars.algo.digitizer.channels;
            
            this.results.phiChCmplx = zeros(xLen, yLen, zLen, repeats, numOfQuant, channels);
            this.results.phiCh      = zeros(xLen, yLen, zLen, repeats, numOfQuant, channels);
            this.results.phiQuant   = zeros(xLen, yLen, zLen, repeats, numOfQuant);
            this.results.phiStd     = zeros(xLen, yLen, zLen, repeats);
            this.results.phi        = zeros(xLen, yLen, zLen, repeats);
            this.results.phiAvg     = zeros(xLen, yLen, zLen);
            this.results.phiAvgStd  = zeros(xLen, yLen, zLen);
            
            this.curScanIdx = zeros(1,3);
        end
        
        function configureScan(this)
            this.fileSystem.turnLogFileOn();
            this.timeTable.scan = struct();
            this.initResultsArrays();
            this.ao.configPeripherals();
            this.graphics.setGraphicsScanVars();
            this.graphics.updateGraphicsConstruction();
        end

        % Manage Variables
        function setUserVars(this, uVars)
            fprintf("S2D: Setting Variables to 2D-Scan.\n");
            this.uVars = uVars();
            
            fprintf("S2D: 1. Setting FileSystem variables.\n");
            
            uVars.fileSystem.saveAny = uVars.fileSystem.saveVars || uVars.fileSystem.saveResults || uVars.fileSystem.saveFigs;
            
            this.fileSystem.setUserVars(uVars.fileSystem);
            this.fileSystemVars = this.fileSystem.configFileSystem();
            
            uVars.ao.fileSystem.saveVars = false; %don't save vars each measurement made
            
            uVars.ao.fileSystem.dataFileNameModel      = "R%dX%.2fY%.2f.mat";
            uVars.ao.fileSystem.varsFileNameModel      = "AOVars.mat";
            
            uVars.ao.fileSystem.extProject         = true;
            uVars.ao.fileSystem.extProjPath        = this.fileSystemVars.projPath;
            uVars.ao.fileSystem.extProjResultsPath = this.fileSystemVars.rawDataPath;
            uVars.ao.fileSystem.extProjFigsPath    = this.fileSystemVars.figsPath;
            
            fprintf("S2D: 3. Setting AcoustoOptivs variables.\n");
            this.ao.setMeasVars(uVars.ao);
            this.ao.fileSystem.setVarsFilenameVariables({});
            
            this.aoVars = this.ao.getAOVars(); 
            
            fprintf("S2D: 4. Setting scan variables.\n");
            this.grid.startX  = uVars.grid.startX;
            this.grid.strideX = uVars.grid.strideX;
            this.grid.endX    = uVars.grid.endX;
            this.grid.startY  = uVars.grid.startY;
            this.grid.strideY = uVars.grid.strideY;
            this.grid.endY    = uVars.grid.endY;
            
            this.grid.firstAxis = uVars.grid.firstAxis;
            
            this.calcGrid();
            
            fprintf("S2D: 5. Setting general variables.\n");
            this.generalVars.useQuant = uVars.general.useQuant;
            this.generalVars.repeats  = uVars.general.repeats;
            
            fprintf("S2D: 6. Setting figures variables.\n");
            figs = this.graphics.createGraphicsUserVars();

            figs.zIdx = uVars.figs.zIdx;

            figs.useQuant = uVars.general.useQuant;
            figs.repeats  = uVars.general.repeats;
            figs.xAxis    = this.grid.xVec ;
            figs.yAxis    = this.grid.yVec;
            figs.zAxis    = this.aoVars.measVars.algo.len.zVecUSRes;
            
            figs.scanFirstAx = uVars.grid.firstAxis;
            
            figs.intExt      = uVars.figs.intExt;
            figs.validStruct = uVars.figs.validStruct;
            figs.extH        = uVars.figs.extH;
            figs.fonts       = uVars.figs.fonts;
            
            this.graphics.setUserVars(figs);
            this.figsVars = figs;
            
            fprintf("S2D: 7. Setting filesystem variables.\n");
%             this.fileSystem.seterVars();
        end
        
        function vars = getScanVars(this)
            vars.ao         = this.aoVars;
            vars.grid       = this.grid;
            vars.general    = this.generalVars;
            vars.fileSystem = this.fileSystemVars;
            vars.uVars      = this.uVars;
            vars.uVars.figs    = scan2DGraphics.createGraphicsUserVars();
            vars.uVars.ao.figs = AOGraphics.createGraphicsUserVars();
        end
        
        % Results Functions
        function putAOResToResultsArray(this, res)
            this.results.phiChCmplx( this.curPosIdx(2), this.curPosIdx(3), :, this.curPosIdx(1), :, :) = permute(gather(res.phiChCmplx), [4, 5, 3, 6, 1, 2]);
            this.results.phiCh(      this.curPosIdx(2), this.curPosIdx(3), :, this.curPosIdx(1), :, :) = permute(gather(res.phiCh),      [4, 5, 3, 6, 1, 2]);
            this.results.phiQuant(   this.curPosIdx(2), this.curPosIdx(3), :, this.curPosIdx(1), :)    = permute(gather(res.phiQuant),   [3, 4, 2, 5, 1]);
            this.results.phi(        this.curPosIdx(2), this.curPosIdx(3), :, this.curPosIdx(1))       = permute(gather(res.phi),        [1, 3, 2, 4]);
            this.results.phiStd(     this.curPosIdx(2), this.curPosIdx(3), :, this.curPosIdx(1))       = permute(gather(res.phiStd),     [1, 3, 2, 4]);
            this.results.phiAvg(     this.curPosIdx(2), this.curPosIdx(3), :)       = permute(gather(res.phi),        [1, 3, 2]);
            this.results.phiAvgStd(  this.curPosIdx(2), this.curPosIdx(3), :)       = permute(gather(res.phiStd),     [1, 3, 2]);
        end

        function averageRepeats(this)
            if strcmp(this.grid.firstAxis, 'Y')
                this.results.phiAvg    = mean(this.results.phi(:,:,:,1:this.curScanIdx(1)), 4);
                this.results.phiAvgStd = std(this.results.phi(:,:,:,1:this.curScanIdx(1)), 0, 4);
            elseif strcmp(this.grid.firstAxis, 'X')
                this.results.phiAvg    = mean(this.results.phi(:,:,:,1:this.curScanIdx(1)), 4);
                this.results.phiAvgStd = std(this.results.phi(:,:,:,1:this.curScanIdx(1)), 0, 4);
            end 
        end
        
        % Scan Functions
        function res = startScan(this)
            this.configureScan();
            if this.fileSystemVars.saveVars
                this.fileSystem.saveVarsToDisk(this.getScanVars(), "");
                this.ao.saveVarsToDisk("");
            end
            
            for r = 1:this.generalVars.repeats
                this.curScanIdx(1) = r;
                for i=1:this.grid.secondIdxLen
                    
                    this.updateCurPosAndIdx( r, 1, i);
                    this.graphics.updateCurPosAndIdx(this.getPosAndScan()); 
                    
                    this.plotAvgPlots()
                    
                    for j=1:this.grid.firstIdxLen
                        
                        this.updateCurPosAndIdx( r, j, i);
                        this.graphics.updateCurPosAndIdx(this.getPosAndScan());
                        this.ao.fileSystem.setDataFilenameVariables({r, this.curPos(1), this.curPos(2)}); 
                        
                        this.printStr(sprintf("S2D: Scaninng on Position: (X,Y) = (%.2f, %.2f)", this.curPos(1), this.curPos(2)), true);
                        this.startScanTime('singlePos');
                        
                        this.startScanTime('netAcoustoOptics');
                        res = this.ao.moveMeasureAndAnalyse(this.curPos);
                        this.stopScanTime('netAcoustoOptics');
                        
                        this.startScanTime('copyTime');
                        this.putAOResToResultsArray(res);
                        this.stopScanTime('copyTime');
                        
                        this.graphics.setData (this.results.phi, this.results.phiStd,...
                                               this.results.phiAvg, this.results.phiAvgStd);
                        this.plotCurPlots();
                        
                        this.stopScanTime('singlePos');
                    end
                end
                this.curScanIdx(2:3) = 0; %for time table

                this.averageRepeats();
                this.plotAvgPlots();
                
                % NOTICE:It doesn't have to be here but it's here so graphics could
                % create a real time mean illustration
            end
            if this.fileSystemVars.saveResults
                this.fileSystem.saveResultsToDisk(this.results);
            end
            res = this.results;
            this.printStr(sprintf("S2D: Done 2D scan.\n"), true);
            this.fileSystem.turnLogFileOff();
        end
        
        % Position Functions
        function updateCurPosAndIdx(this, r, first, second)
            % curScanIdx - always aligned to (R, 1st, 2nd)[#]
            % curScan    - always aligned to ( 1st, 2nd)[mm]
            this.curScanIdx  = [r, first, second];
            this.curScan     = [this.grid.firstVec(first), this.grid.secondVec(second)];
            
            % curPosIdx - always aligned to (R, xIdx, yIdx)[#]
            % curPos    - always aligned to (X, Y] [mm]
            if strcmp(this.grid.firstAxis, 'Y')
                this.curPosIdx = [r, second, first];
                this.curPos    = [this.grid.secondVec(this.curPosIdx(2)), ...
                                  this.grid.firstVec(this.curPosIdx(3))];
            elseif strcmp(this.grid.firstAxis, 'X')
                this.curPosIdx = [r,first, second];
                this.curPos    = [this.grid.firstVec(this.curPosIdx(2)), ...  
                                  this.grid.secondVec(this.curPosIdx(3))];
            end
        end
        
        function pns = getPosAndScan(this)
            pns.curPos     = this.curPos;
            pns.curScan    = this.curScan;
            pns.curPosIdx  = this.curPosIdx;
            pns.curScanIdx = this.curScanIdx;
        end
        
        function calcGrid(this)
            this.grid.xVec    = this.grid.startX:this.grid.strideX:this.grid.endX;
            this.grid.yVec    = this.grid.startY:this.grid.strideY:this.grid.endY;
            
            this.grid.xLen    = abs(this.grid.startX - this.grid.endX);
            this.grid.yLen    = abs(this.grid.startY - this.grid.endY);
            this.grid.xIdxLen = length(this.grid.xVec);
            this.grid.yIdxLen = length(this.grid.yVec);
            this.grid.xIdx    = 1:1:this.grid.xIdxLen;        
            this.grid.yIdx    = 1:1:this.grid.yIdxLen;
            
            switch this.grid.firstAxis
                case 'Y'
                    this.grid.startFirst   = this.grid.startY;
                    this.grid.strideFirst  = this.grid.strideY;
                    this.grid.endFirst     = this.grid.endY;
                    
                    this.grid.startSecond  = this.grid.startX;
                    this.grid.strideSecond = this.grid.strideX;
                    this.grid.endSecond    = this.grid.endX;
                    
                    this.grid.firstVec     = this.grid.yVec;
                    this.grid.secondVec    = this.grid.xVec;
                    
                    this.grid.firstLen     = this.grid.yLen;
                    this.grid.secondLen    = this.grid.xLen;
                    this.grid.firstIdxLen  = this.grid.yIdxLen;
                    this.grid.secondIdxLen = this.grid.xIdxLen;
                    this.grid.firstIdx     = this.grid.yIdx;
                    this.grid.secondIdx    = this.grid.xIdx;
                    
                case 'X'
                    this.grid.startFirst   = this.grid.startX;
                    this.grid.strideFirst  = this.grid.strideX;
                    this.grid.endFirst     = this.grid.endX;
                    
                    this.grid.startSecond  = this.grid.startY;
                    this.grid.strideSecond = this.grid.strideY;
                    this.grid.endSecond    = this.grid.endY;
                    
                    this.grid.firstVec     = this.grid.xVec;
                    this.grid.secondVec    = this.grid.yVec;
                    
                    this.grid.firstLen     = this.grid.xLen;
                    this.grid.secondLen    = this.grid.yLen;
                    this.grid.firstIdxLen  = this.grid.xIdxLen;
                    this.grid.secondIdxLen = this.grid.yIdxLen;
                    this.grid.firstIdx     = this.grid.xIdx;
                    this.grid.secondIdx    = this.grid.yIdx;
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
        
        function setLoadedData(this, data, vars)
            setLoadedVars(this, vars.uVars);
            this.results = data;
            this.graphics.setGraphicsScanVars();
            this.graphics.updateGraphicsConstruction();
            this.graphics.setData(data.phi, data.phiStd, data.phiAvg, data.phiAvgStd);
        end
    end
end

