classdef scanAO < handle
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
        
        res
%         res1D
        timeTable
        
        % Vars
        uVars
        grid
        aoVars
        figsVars
        fileSystemVars
        generalVars
        
        curPos   % [first, second] [mm]
        curIdx   % [firstIdx, secondIdx, r] [#]
        timing
        owned
        extStages
    end
    
    methods (Static) 
        function uVars = createUserVars()
            uVars.ao         = acoustoOptics.createUserVars();
            uVars.figs       = AOSGraphics.createUserVars();
            uVars.fileSystem = fileSystemAOS.uVarsFileSystemCreate();

            uVars.grid.scan1Start     = 0;
            uVars.grid.scan1End       = 0;
            uVars.grid.scan1Stride    = 0;
            
            uVars.grid.scan2Start     = 0;
            uVars.grid.scan2End       = 0;
            uVars.grid.scan2Stride    = 0;

            uVars.grid.thirdPos = 0;
            
            uVars.grid.scan1Lable  = 'Y';
            uVars.grid.scan2Label  = 'X';
            uVars.grid.thirdLabel  = 'Z';
            uVars.grid.depthLabel  = 'Z';     
            
            uVars.stages.moveNonScanStages   = true;
            uVars.stages.stagesOwnedByParent = false;
            
            uVars.general.externalScan      = true;
        end
    end
    
    methods
        function this = scanAO(acoustoOpticHandle, stagesHandle)
            this.sf = statsFunctions("AOS");
            % Acousto optics initialization must be conducted before stages
            % connect since AO includes instrreset, what will damage stages
            % connection if called perliminary
            if nargin>0 && ~isempty(acoustoOpticHandle)
                this.ao = acoustoOpticHandle;
            else
                fprintf("AOS: Creating and connecting to Acousto Optics object...\n");
                this.ao = acoustoOptics();
                this.ao.init();
            end
            
            if nargin>1
                this.stages = stagesHandle;
                this.extStages = true;
            else
                fprintf("AOS: Creating and connecting to a stages object...\n");
                this.stages = stages('Zaber', 'COM3');
                this.extStages = false;
                this.stages.connect();
                fprintf("AOS: Done connecting to stages.\n");
            end
            this.graphics = AOSGraphics();
            this.fileSystem = fileSystemAOS(this, this.ao.fileSystem);
        end
 
        % Manage Variables
        function setUserVars(this, uVars)
            fprintf("AOS: Setting Variables to 2D-Scan.\n");
            this.uVars = uVars;
            
            %-------------
            % General
            %-------------
            fprintf("AOS: 1. Setting general variables.\n"); 
            this.generalVars.externalScan          = uVars.general.externalScan;
            this.generalVars.moveNonScanStages     = uVars.stages.moveNonScanStages;
            this.generalVars.singleChannelAnalysis = uVars.ao.ao.analyzeSingleCh;
            
            this.generalVars.tg.send = uVars.general.tgSend;
            this.generalVars.tg.chatID = uVars.general.chatID;

            %-------------
            % File System
            %-------------
            fprintf("AOS: 2. Setting AOS FileSystem variables.\n");
            this.fileSystem.setUserVars(uVars.fileSystem);
            
            %---------------
            % Acousto Optics
            %---------------
            fprintf("AOS: 3. Setting AcoustoOptics variables.\n");
            uVars.ao.fileSystem.extProject   = true;
            uVars.ao.fileSystem.saveResults  = uVars.fileSystem.saveResults;
            uVars.ao.fileSystem.dontSaveVars = true;
            
            uVars.ao.figs.intExt   = uVars.figs.intExt;
            uVars.ao.reopenFigures = uVars.figs.reopenFigures;
            
            uVars.ao.ao.uploadToTelegram = false;

            this.ao.setVars(uVars.ao);
            this.aoVars = this.ao.getVars();

            %---------------
            % Grid
            %---------------
            fprintf("AOS: 4. Setting Grid variables.\n");
            
            this.grid.scan1Start  = uVars.grid.scan1Start;
            this.grid.scan1Stride = uVars.grid.scan1Stride;
            this.grid.scan1End    = uVars.grid.scan1End;
             
            this.grid.scan2Start  = uVars.grid.scan2Start;
            this.grid.scan2Stride = uVars.grid.scan2Stride;
            this.grid.scan2End    = uVars.grid.scan2End;

            this.grid.thirdPos = uVars.grid.thirdPos;

            this.grid.scan1Label = uVars.grid.scan1Label;
            this.grid.scan2Label = uVars.grid.scan2Label;
            this.grid.thirdLabel = uVars.grid.thirdLabel;
            this.grid.depthLabel = uVars.grid.depthLabel;
            
            this.calcGrid();

            %---------------
            % Stages
            %---------------
            fprintf("AOS: 5. Set Stages variables and Axes.\n");
            % only in case stages are owned by a parent object dont perform
            % calibration
            this.generalVars.stagesOwnedByParent = uVars.stages.stagesOwnedByParent;
            if ~this.generalVars.stagesOwnedByParent
                this.stages.assignStagesAxes([this.grid.scan1Label, this.grid.thirdLabel]);
            end

            fprintf("AOS: 5. Setting strings models and variables\n");
                       
            strings{1} = "scan1Done";
            models{1}  = sprintf("AO on pos (%s,%s,%s) = (%s, %s, %s)", this.grid.scan1Label, this.grid.scan2Label, this.grid.thirdLabel, "%.2f", "%.2f", "%.2f");
            strings{2} = "scan2Start";
            models{2}  = sprintf("Startin New Line: %s = %s", this.grid.scan2Label, "%.2f");
            this.sf.setStringModels(strings, models);
            
            %---------------
            % Figures
            %---------------
            fprintf("AOS: 6. Setting figures variables.\n");
            figs = this.graphics.createOwnerVars();

            figs.channelsInRes         = this.aoVars.measVars.algo.general.analysisReps;
            figs.singleChannelAnalysis = this.generalVars.singleChannelAnalysis;
            figs.sepChIdx              = uVars.figs.sepChIdx;
            figs.normDispPlaneColors   = uVars.figs.normDispPlaneColors;
            figs.reopenFigures         = uVars.figs.reopenFigures;

            figs.scan1AxType  = uVars.figs.scan1AxType;
            figs.scan1VecNorm = this.grid.scan1Vec;

            figs.scan2AxType  = uVars.figs.scan2AxType;
            figs.scan2VecNorm = this.grid.scan2Vec;           

            figs.depthAxType  = uVars.figs.depthAxType;
            figs.depthVecNorm = this.aoVars.measVars.algo.len.depthVec;
            
            figs.thirdPos = this.grid.thirdPos;
                        
            figs.scan1Label  = uVars.grid.scan1Label;
            figs.scan2Label  = uVars.grid.scan2Label;
            figs.thirdLabel  = uVars.grid.thirdLabel;
            figs.depthLabel  = uVars.grid.depthLabel;
            figs.mainPlane     = this.grid.mainPlane;
            
            figs.varsNavLine  = uVars.figs.varsNavLine;
            figs.varsNavPlane = uVars.figs.varsNavPlane;
            
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
            
            this.curPos = [this.grid.scan1Start, this.grid.scan2Start, this.grid.thirdPos];
            this.curIdx = [1, 1];
            
            %Config subObjects
            this.fileSystem.configFileSystem(this.grid.scan1Label, this.grid.scan2Label);
            this.fileSystem.saveVarsToDisk();
            this.ao.configPeripherals();
        end
        
        % Results Functions
        function initResultsArrays(this)      
            resArr = this.getResultsArrayModel();
            
            this.res.phi        = resArr.phi;
            this.res.phiNorm    = resArr.phi;
            this.res.phiLog     = resArr.phi;
           
            this.res.rawPhi     = resArr.rawPhi;
           
            separateCh        = this.aoVars.measVars.algo.general.analysisReps;

            this.res.globalMin  = inf*ones(1,separateCh);
            this.res.globalMax  = zeros(1,separateCh);
            this.res.globalSpan = zeros(1,separateCh);

            this.res.globalMinLog  = zeros(1,separateCh);

            this.curIdx = zeros(1,2);
            this.timing.progress = 0;
            this.timing.averageTime = 0;
        end

        function resArr = getResultsArrayModel(this)
            scan1Len   = this.grid.scan1IdxLen;
            scan2Len   = this.grid.scan2IdxLen;
            depthLen   = this.aoVars.measVars.algo.samples.numOfPos;
            separateCh = this.aoVars.measVars.algo.general.analysisReps;
           
            resArr.phi        = zeros(scan1Len, scan2Len, depthLen, separateCh);
            resArr.phiAvg     = zeros(scan1Len, scan2Len, depthLen, separateCh);
            resArr.phiAvgStd  = zeros(scan1Len, scan2Len, depthLen, separateCh);
            
            resArr.rawPhi     = zeros(scan1Len, scan2Len, depthLen, separateCh);
            resArr.phiLog     = zeros(scan1Len, scan2Len, depthLen, separateCh);
            resArr.phiNorm    = zeros(scan1Len, scan2Len, depthLen, separateCh);
        end

        function putAOToResArray(this, res)
            separateCh        = this.aoVars.measVars.algo.general.analysisReps;
            
            for i = 1:separateCh
                phi = res(i).phi;
                this.res.phi(   this.curIdx(1), this.curIdx(2), :, i) = phi;
                this.res.rawPhi(this.curIdx(1), this.curIdx(2), :, i) = res(i).rawPhi;
                
                phi = this.res.phi(:,:,:,i);
                % Global normalization and elimination of zeros
                globalMin  = min(this.res.globalMin(i), min(phi(:)));
                globalMax  = max(this.res.globalMax(i), max(phi(:)));
                globalSpan = globalMax - globalMin;

                this.res.globalMin(i)  = globalMin;
                this.res.globalMax(i)  = globalMax;
                this.res.globalSpan(i) = globalSpan;

                phiNorm = (phi - globalMin) ./ globalSpan;
                minVals = mink(unique(phiNorm(:)), 2);
                phiNorm(phiNorm == 0) =  minVals(2);
                this.res.phiNorm(:,:,:,i) = phiNorm;

                % GLobal Log scale
                this.res.phiLog(:,:,:,i) = log(this.res.phiNorm());
                this.res.globalMinLog(i) =  log(minVals(2));
            end
        end
   
        function res = scan(this)
            this.initResultsArrays();
            this.graphics.resetGraphics();
            this.graphics.resetPlots();
            this.resetTimeTable();

            this.updateCurPosAndIdx( 1, 1);
            if this.generalVars.moveNonScanStages
                this.stages.moveAbsAx(this.grid.thirdLabel, this.grid.thirdPos)
            end
            
            i=1;

            for j=1:this.grid.scan2IdxLen
                this.updateCurPosAndIdx(i, j);
                this.stages.moveAbsAx(this.grid.scan2Label, this.curPos(2))
                
                curScan1Start = this.grid.scan1StartVec(j);
                curScan1Step  = this.grid.scan1StrideVec(j);
                curScan1End   = this.grid.scan1EndVec(j);

                this.sf.printStrModel("scan2Start", this.curPos(2), true);
                fprintf("-----------------------\n");
                for i = curScan1Start : curScan1Step : curScan1End
                    T = tic;
                    this.timing.progress = this.timing.progress + 1;
                    this.updateCurPosAndIdx(i, j);
                    this.sf.startScanTime("ttScan1", 'singlePos', [i, j]);
                    
                    
                    this.sf.startScanTime("ttScan1", 'prepare', [i, j]);
                    this.fileSystem.updateFileSystem(this.curPos(1:2));
                    this.sf.stopScanTime("ttScan1", 'prepare', [i, j]);
                    
                    this.stages.moveAbsAx(this.grid.scan1Label, this.curPos(1))
    
                    this.sf.startScanTime("ttScan1", 'AO', [i, j]);
                    res = this.ao.runAcoustoOptics();
                    this.sf.stopScanTime("ttScan1", 'AO', [i, j]);
    
                    this.sf.startScanTime("ttScan1", 'processing', [i, j]);
                    this.putAOToResArray(res);
                    this.sf.stopScanTime("ttScan1", 'processing', [i, j]);
                    
                    this.sf.startScanTime("ttScan1", 'plotting', [i, j]);
                    this.plotResults();
                    this.sf.stopScanTime("ttScan1", 'plotting', [i, j]);
                    
                    this.sf.stopScanTime("ttScan1", 'singlePos', [i, j]);

                    this.sf.printStrModel("scan1Done", this.curPos, true);
                    T = toc(T);
                    this.printProgress(T);
                end
                this.sendFigure();
            end
           
            this.fileSystem.saveResultsToDisk(this.res);
            this.fileSystem.closeFileSystem();
            res = this.res;
            
            this.sf.printStr(sprintf("Done AO Scan.\n"), true);
        end
        
        % Position Functions
        function updateCurPosAndIdx(this, scan1Idx, scan2Idx)
            % curPosIdx - always aligned to (R, 1st, 2nd)[#]
            % curPos    - always aligned to ( 1st, 2nd)[mm]
            this.curIdx(1) = scan1Idx;
            this.curIdx(2) = scan2Idx;

            this.curPos(1) = this.grid.scan1Vec(scan1Idx);
            this.curPos(2) = this.grid.scan2Vec(scan2Idx);
            
            this.graphics.updateCurPosAndIdx(this.getPos());
        end

        function pns = getPos(this)
            pns.curPos = this.curPos;
            pns.curIdx = this.curIdx;
        end

        function calcGrid(this)
            if this.grid.scan1Stride == 0
                this.grid.scan1Vec = this.grid.scan1Start;
                this.grid.scan1Len     = 0;
            else
                this.grid.scan1Vec = this.grid.scan1Start  : this.grid.scan1Stride  : this.grid.scan1End;
                this.grid.scan1Len = abs(this.grid.scan1Start - this.grid.scan1End);
            end
            this.grid.scan1IdxLen  = length(this.grid.scan1Vec);
            this.grid.scan1Idx     = 1:1:this.grid.scan1IdxLen;        
            this.grid.scan1Cntr    = this.grid.scan1Vec - mean(this.grid.scan1Vec);
            this.grid.scan1Zero    = abs(this.grid.scan1Vec - this.grid.scan1Vec(1));
            
            if this.grid.scan2Stride == 0
                this.grid.scan2Vec = this.grid.scan2Start;
                this.grid.scan2Len = 0;
            else
                this.grid.scan2Vec = this.grid.scan2Start  : this.grid.scan2Stride  : this.grid.scan2End;
                this.grid.scan2Len = abs(this.grid.scan2Start - this.grid.scan2End);
            end
            this.grid.scan2IdxLen  = length(this.grid.scan2Vec);
            this.grid.scan2Idx     = 1:1:this.grid.scan2IdxLen;        
            this.grid.scan2Cntr    = this.grid.scan2Vec - mean(this.grid.scan2Vec);
            this.grid.scan2Zero    = abs(this.grid.scan2Vec - this.grid.scan2Vec(1));

            this.grid.depthVec     = this.aoVars.measVars.algo.len.depthVec;
            dDepth                 = abs(this.grid.depthVec(2) - this.grid.depthVec(1));
            this.grid.depthLen     = this.grid.depthVec(end)+dDepth;
            this.grid.depthIdxLen  = length(this.grid.depthVec);
            this.grid.depthIdx     = 1:1:this.grid.depthIdxLen;
            this.grid.depthCntr    = this.grid.depthVec - mean(this.grid.depthVec);
            this.grid.depthZero    = abs(this.grid.depthVec - this.grid.depthVec(1));
            
            this.grid.totalScanPts = this.grid.scan2IdxLen * this.grid.scan1IdxLen;

            % Snake scan
            this.grid.scan1StartVec  = ones(1, this.grid.scan2IdxLen);
            this.grid.scan1StrideVec = ones(1, this.grid.scan2IdxLen);
            this.grid.scan1EndVec    = ones(1, this.grid.scan2IdxLen)*this.grid.scan1IdxLen;

            this.grid.scan1StartVec(2:2:end)  = this.grid.scan1IdxLen;
            this.grid.scan1StrideVec(2:2:end) = -1;
            this.grid.scan1EndVec(2:2:end)    = 1;

            switch this.grid.scan1Label
                case 'X'
                    switch this.grid.scan2Label
                        case 'Y'
                            this.grid.mainPlane = 'XZ';
                        case 'Z'
                            this.grid.mainPlane = 'XY';
                    end
                case 'Y'
                    switch this.grid.scan2Label
                        case 'X'
                            this.grid.mainPlane = 'YZ';
                        case 'Z'
                            this.grid.mainPlane = 'XY';
                    end   
                case 'Z'
                    switch this.grid.scan2Label
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
        
        function plotResults(this)
            this.graphics.setData(this.res)
            this.graphics.dispAll();
        end
        
        function resetTimeTable(this)
            this.sf.timeTable.ttScan1.singlePos  = zeros(this.grid.scan1IdxLen, this.grid.scan2IdxLen);
            this.sf.timeTable.ttScan1.prepare    = zeros(this.grid.scan1IdxLen, this.grid.scan2IdxLen);
            this.sf.timeTable.ttScan1.AO         = zeros(this.grid.scan1IdxLen, this.grid.scan2IdxLen);
            this.sf.timeTable.ttScan1.processing = zeros(this.grid.scan1IdxLen, this.grid.scan2IdxLen);
            this.sf.timeTable.ttScan1.plotting   = zeros(this.grid.scan1IdxLen, this.grid.scan2IdxLen);

            this.sf.timeTable.tic.ttScan1.singlePos  = zeros(this.grid.scan1IdxLen, this.grid.scan2IdxLen, 'uint64');
            this.sf.timeTable.tic.ttScan1.prepare    = zeros(this.grid.scan1IdxLen, this.grid.scan2IdxLen, 'uint64');
            this.sf.timeTable.tic.ttScan1.AO         = zeros(this.grid.scan1IdxLen, this.grid.scan2IdxLen, 'uint64');
            this.sf.timeTable.tic.ttScan1.processing = zeros(this.grid.scan1IdxLen, this.grid.scan2IdxLen, 'uint64');
            this.sf.timeTable.tic.ttScan1.plotting   = zeros(this.grid.scan1IdxLen, this.grid.scan2IdxLen, 'uint64');

            this.sf.timeTable.ttScan2.scan1Plane = zeros(1, this.grid.scan2IdxLen, 'uint64');
        end
        
        function printProgress(this, T)
           totalNumOfLines = this.grid.totalScanPts;
           curNumOfLines   = this.timing.progress;
           percent = curNumOfLines/totalNumOfLines * 100;

           this.timing.averageTime = (this.timing.averageTime *(curNumOfLines-1) + T )/curNumOfLines;
           
           totalTime = totalNumOfLines * this.timing.averageTime / 60;
           timeToNow = curNumOfLines * this.timing.averageTime / 60;

           str = sprintf("Done line %d/%d (%.2f [%%%%]). Progress: %.2f/%.2f mins.\n",...
                            curNumOfLines, totalNumOfLines, percent, timeToNow, totalTime);
           fprintf(str)
           this.generalVars.tg.str = str;
        end

        function sendFigure(this)
             if this.generalVars.tg.send 
                tgprintf(this.generalVars.tg.chatID, this.generalVars.tg.str);
                tgprint(this.generalVars.tg.chatID, this.graphics.figs.navPlane.handles.cur.ax, 'photo');
                tgprint(this.generalVars.tg.chatID, this.graphics.figs.navPlaneLog.handles.cur.ax, 'photo');
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
            this.res = data;
            this.graphics.setGraphicsScanVars();
            this.graphics.updateGraphicsConstruction();
            this.graphics.setLoadedData(data);
        end
    end
end