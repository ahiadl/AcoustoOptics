classdef scan2D < handle
    %SCAN2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ao
        graphics
        fileSystem
        sf
        owner
        
        results
        timeTable
        
        % Vars
        uVars
        grid
        aoVars
        figsVars
        fileSystemVars
        generalVars
        
        strings
        curPos      %[X,Y]          [mm]
        curPosIdx   %[r, xIdx, yIdx][#]
        
        owned
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
           
            uVars.grid.startFirst     = 0;
            uVars.grid.startSecond    = 0;
            uVars.grid.endFirst       = 0;
            uVars.grid.endSecond      = 0;
            uVars.grid.strideFirst    = 0;
            uVars.grid.strideSecond   = 0;
            
            uVars.grid.firstAxis  = 'Y';
            uVars.grid.secondAxis = 'X';
%             uVars.grid.depthAxis  = 'Z';
            
            uVars.general.repeats   = 1;
            uVars.general.useQuant  = true;
        end
        
    end
    
    methods
        function this = scan2D(acoustoOpticHandle)
            strings.scan      = "Done Scan for (R,X,Y) = (%d, %.2f, %.2f)";
            strings.timeTable = "R%dX%.2fY%.2f";
            
            this.sf = statsFunctions("S2D", strings);
            
            if ~isempty(acoustoOpticHandle)
                this.ao = acoustoOpticHandle;
            else
                this.ao = acoustoOptics();
            end

            this.graphics = scan2DGraphics();
            this.fileSystem = fileSystemS2D();
        end
        
        function initResultsArrays(this)      
            firstAxLen     = this.grid.firstIdxLen;
            secondAxLen    = this.grid.secondIdxLen;
            depthLen       = this.aoVars.measVars.algo.samples.numOfPos;
            
            repeats    = this.generalVars.repeats;
            numOfQuant = this.aoVars.measVars.algo.samples.numOfQuant;
            channels   = this.aoVars.measVars.algo.digitizer.channels;
            
            this.results.phiChCmplx = zeros(firstAxLen, secondAxLen, depthLen, repeats, numOfQuant, channels);
            this.results.phiCh      = zeros(firstAxLen, secondAxLen, depthLen, repeats, numOfQuant, channels);
            this.results.phiQuant   = zeros(firstAxLen, secondAxLen, depthLen, repeats, numOfQuant);
            this.results.phiStd     = zeros(firstAxLen, secondAxLen, depthLen, repeats);
            this.results.phi        = zeros(firstAxLen, secondAxLen, depthLen, repeats);
            this.results.phiAvg     = zeros(firstAxLen, secondAxLen, depthLen);
            this.results.phiAvgStd  = zeros(firstAxLen, secondAxLen, depthLen);
            
            this.curPosIdx = zeros(1,3);
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
            
            uVars.ao.fileSystem.dataFileNameModel      = "AO-R%dX%.2fY%.2f.mat";
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
            this.grid.startFirst   = uVars.grid.startFirst;
            this.grid.strideFirst  = uVars.grid.strideFirst;
            this.grid.endFirst     = uVars.grid.endFirst;
            this.grid.startSecond  = uVars.grid.startSecond;
            this.grid.strideSecond = uVars.grid.strideSecond;
            this.grid.endSecond    = uVars.grid.endSecond;
            
            this.grid.firstAxis  = uVars.grid.firstAxis;
            this.grid.secondAxis = uVars.grid.secondAxis;

            this.calcGrid();
            
            fprintf("S2D: 5. Setting general variables.\n");
            this.generalVars.useQuant = uVars.general.useQuant;
            this.generalVars.repeats  = uVars.general.repeats;
            
            fprintf("S2D: 6. Setting figures variables.\n");
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
                this.results.phiAvg    = mean(this.results.phi(:,:,:,1:this.curPosIdx(1)), 4);
                this.results.phiAvgStd = std(this.results.phi(:,:,:,1:this.curPosIdx(1)), 0, 4);
            elseif strcmp(this.grid.firstAxis, 'X')
                this.results.phiAvg    = mean(this.results.phi(:,:,:,1:this.curPosIdx(1)), 4);
                this.results.phiAvgStd = std(this.results.phi(:,:,:,1:this.curPosIdx(1)), 0, 4);
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
                this.curPosIdx(1) = r;
                for i=1:this.grid.secondIdxLen
                    
                    this.updateCurPosAndIdx( r, 1, i);
                    this.graphics.updateCurPosAndIdx(this.getPos()); 
                    
                    this.plotAvgPlots()
                    
                    for j=1:this.grid.firstIdxLen
                        
                        this.updateCurPosAndIdx( r, j, i);
                        this.graphics.updateCurPosAndIdx(this.getPos());
                        this.ao.fileSystem.setDataFilenameVariables({r, this.curPos(1), this.curPos(2)}); 
                        
                        this.sf.printStr(sprintf("Scaninng on Position: (X,Y) = (%.2f, %.2f)", this.curPos(1), this.curPos(2)), true);
                        this.sf.startScanTime('singlePos', this.curPos);
                        
                        this.sf.startScanTime('netAcoustoOptics', this.curPos);
                        res = this.ao.moveMeasureAndAnalyse(this.curPos);
                        this.sf.stopScanTime('netAcoustoOptics', this.curPos);
                        
                        this.sf.startScanTime('copyTime', this.curPos);
                        this.putAOResToResultsArray(res);
                        this.sf.stopScanTime('copyTime', this.curPos);
                        
                        this.graphics.setData (this.results.phi, this.results.phiStd,...
                                               this.results.phiAvg, this.results.phiAvgStd);
                        this.plotCurPlots();
                        
                        this.sf.stopScanTime('singlePos', this.curPos);
                    end
                end
                this.curPosIdx(2:3) = 0; %for time table

                this.averageRepeats();
                this.plotAvgPlots();
                
                % NOTICE:It doesn't have to be here but it's here so graphics could
                % create a real time mean illustration
            end
            if this.fileSystemVars.saveResults
                this.fileSystem.saveResultsToDisk(this.results);
            end
            res = this.results;
            this.sf.printStr(sprintf("Done 2D scan.\n"), true);
            this.fileSystem.turnLogFileOff();
        end
        
        % Position Functions
        function updateCurPosAndIdx(this, r, first, second)
            % curScanIdx - always aligned to (R, 1st, 2nd)[#]
            % curScan    - always aligned to ( 1st, 2nd)[mm]
            this.curPosIdx  = [r, first, second];
            this.curPos     = [this.grid.firstVec(first), this.grid.secondVec(second)];
        end
        
        function pns = getPos(this)
%             pns.curPos     = this.curPos;
            pns.curPos    = this.curPos;
%             pns.curPosIdx  = this.curPosIdx;
            pns.curPosIdx = this.curPosIdx;
        end
        
        function calcGrid(this)
            this.grid.firstVec  = this.grid.startFirst  : this.grid.strideFirst  : this.grid.endFirst;
            this.grid.secondVec = this.grid.startSecond : this.grid.strideSecond : this.grid.endSecond;
            this.grid.depthVec = this.aoVars.measVars.algo.len.zVecUSRes;
            
            this.grid.firstLen     = abs(this.grid.startFirst - this.grid.endFirst);
            this.grid.secondLen    = abs(this.grid.startSecond - this.grid.endSecond);
            this.grid.firstIdxLen  = length(this.grid.firstVec);
            this.grid.secondIdxLen = length(this.grid.secondVec);
            this.grid.firstIdx     = 1:1:this.grid.firstIdxLen;        
            this.grid.secondIdx    = 1:1:this.grid.secondIdxLen;
            
            this.grid.firstNorm    = this.grid.firstVec - mean(this.grid.firstVec);
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

