classdef pulseEchoUS < handle
    % Pulse Echo Ultrasound Object: performs an ultrasound imaging by
    % projecting an ultrasound pulse and measure its reflections
    % From each measurement a line reconstruction is made. 
    % By performing 1D scan a 2D image is received. 
    % In order to get best results, a pulse intensity calibration is
    % needed.
    
    properties
        scope
        stages
        graphics
        fileSystem
        sf
        owner
        
        res1D
        res2D
        timeTable
        
        % Vars
        grid
        uVars
        extVars
        measVars
        
        curPos      %[scan, third]      [mm]
        curPosIdx   %[scanIdx, thirdIdx][#]
        
        owned
        extStages
    end
    
    methods (Static) 
        function uVars = uVarsCreate()
            % Recons
            uVars.distFromCenter = [];
            uVars.imageWidth     = [];
            uVars.c              = [];
            uVars.usFreq         = [];
            
            % Scope
            uVars.signalChannel = [];
            uVars.trigCh        = [];
            uVars.avg           = [];
            uVars.usRepRate     = [];
            uVars.pulseIntVal   = [];
            
            % Grid & Stages
            uVars.scanStart  = [];
            uVars.scanStride = [];
            uVars.scanEnd    = [];
            
            uVars.thirdPos   = [];

            uVars.scanAxis   = 'X';  
            uVars.depthAxis  = 'Y';
            uVars.thirdAxis  = 'Z';
            
            uVars.stagesOwnedByParent = false;
            
            % Graphics
            uVars.figs = PEUSGraphics.createUserVars();
            
            % FileSystem
            uVars.fileSystem = fileSystemPEUS.uVarsCreate();
        end
    end
    
    methods
        function this = pulseEchoUS(scopeHandle, stagesHandle)
            this.sf = statsFunctions("PEUS");
            % Acousto optics initialization must be conducted before stages
            % connect since AO includes instrreset, what will damage stages
            % connection if called perliminary
            if nargin>0 && ~isempty(scopeHandle)
                this.scope = scopeHandle;
            else
                fprintf("PEUS: Creating and connecting to Scope...\n");
                Scope_visa_address = 'USB0::0x1AB1::0x04CE::DS1ZA172115085::0::INSTR';
                this.scope = Rigol_DS1074z('ni',Scope_visa_address);
            end
            
            if nargin>1
                this.stages = stagesHandle;
                this.extStages = true;
            else
                fprintf("PEUS: Creating and connecting to Stages...\n");
                this.stages = stages('COM3');
                this.extStages = false;
                this.stages.connect();
                fprintf("PEUS: Done connecting to stages.\n");
            end
            this.graphics = PEUSGraphics();
            this.fileSystem = fileSystemPEUS(this);
        end
 
        % Manage Variables
        function setVars(this, uVars)
            this.uVars = uVars;

            % Reconstruction parameters
            reconVars.distFromCenter = uVars.distFromCenter;
            reconVars.imageWidth     = uVars.imageWidth;
            reconVars.c              = uVars.c;
            reconVars.usFreq         = uVars.usFreq;
            reconVars.pulseIntVal    = uVars.pulseIntVal;
            
            reconVars.delay  = 2*(reconVars.distFromCenter/reconVars.c);
            reconVars.tScale = ((reconVars.imageWidth)/reconVars.c)/6;

            % Scope parameters
            scopeVars.ch.id = uVars.signalChannel;
            scopeVars.ch.coup   = 'AC';
            scopeVars.ch.invert = 'OFF';
            scopeVars.ch.offset = 0;
            scopeVars.ch.scale  = 0.005; %TODO: Logic to set the scale - take it from optimizeResolution in boardTesting project
                      
            scopeVars.trigger.mode      = 'EDGE';
            scopeVars.trigger.trigCh    = uVars.trigCh;
            scopeVars.trigger.edge      = 'POSitive';
            scopeVars.trigger.trigLevel = 1; %TODO: should be configured according to pulser trigger level/
            
            scopeVars.timebase.mode   = 'MAIN';
            scopeVars.timebase.scale  = reconVars.tScale;
            scopeVars.timebase.offset = reconVars.delay;
            
            scopeVars.waveform.chToSample    = uVars.signalChannel;
            scopeVars.waveform.WavPointsMode = 'NORMal';
            scopeVars.waveform.Format        = 'BYTE';
            
            factor = 2;
            scopeVars.acq.mode           = 'AVER';
            scopeVars.acq.numOfAvg       = uVars.avg;
            scopeVars.acq.waitAfterClear = (1/uVars.usRepRate)*uVars.avg*factor;
            
            this.scope.setGlobalConf(scopeVars);
            this.scope.setChannelVars(scopeVars.ch);
            this.configScope(scopeVars.waveform.chToSample); % TODO: think how to move this function to configure function
            
            this.extVars.scope  = scopeVars;
            this.measVars.scope = this.scope.getVars();
            
            reconVars.tVec     = this.measVars.scope.timebase.dispTimeVec;
            reconVars.fs       = 1./this.measVars.scope.timebase.dt;
            reconVars.depthVec = (1/2)*(reconVars.c*(reconVars.tVec));

            this.extVars.recon = reconVars;            
            
            % Grid & Stages parameters
            
            this.grid.scanStart   = uVars.scanStart;
            this.grid.scanStride  = uVars.scanStride;
            this.grid.scanEnd     = uVars.scanEnd;
            
            this.grid.thirdPos    = uVars.thirdPos;

            this.grid.scanAxis   = uVars.scanAxis;  %Label
            this.grid.depthAxis  = uVars.depthAxis; %Label

            this.calc2DGrid();
            
            stagesVars.moveSecondStage     = uVars.moveSecondStage;
            stagesVars.stagesOwnedByParent = uVars.stagesOwnedByParent;
            if ~stagesVars.stagesOwnedByParent
                this.stages.assignStagesAxes([this.grid.scanAxis, this.grid.thirdAxis]);
            end
            
            this.extVars.stages = stagesVars;
            this.measVars.stages = stagesVars;
            
            this.measVars.recon = reconVars;
            
            % Stats Functions
            
            strings{1} = "start2DScan";
            models{1}  = sprintf("Start Scan for (%s) = (%s)",  this.grid.thirdAxis, "%.2f");
            strings{2} = "start1DScan";
            models{2}  = sprintf("Start Scan for (%s,%s) = (%s, %s)", this.grid.scanAxis, this.grid.thirdAxis, "%.2f", "%.2f");
            strings{3} = "timeTable2D";
            models{3}  = sprintf("%s%s", this.grid.thirdAxis, "%.2f");
            strings{4} = "timeTable1D";
            models{4}  = sprintf("%s%s", this.grid.scanAxis, "%.2f");
            this.sf.setStringModels(strings, models);
            
            % Graphics parameters
            figsVars = PEUSGraphics.createOwnerVars();
            
            figsVars.intExt        = uVars.figs.intExt;
            figsVars.useExtClims   = uVars.figs.useExtClims;
            figsVars.reopenFigures = uVars.figs.reopenFigures;
            
            figsVars.numOfSamples  = 1200; %TODO: get this parameter from the scope
            figsVars.imageWidth    = uVars.imageWidth;
            
            figsVars.tVec          = this.measVars.recon.tVec;
            
            figsVars.depthAxis     = this.grid.depthVec;
            figsVars.depthAxisCntr = this.grid.depthCntr;
            
            figsVars.scanAxis      = this.grid.scanVec;
            figsVars.scanAxisCntr  = this.grid.scanCntr;
            
            figsVars.depthAxLabel  = this.grid.depthAxis;
            figsVars.scanAxLabel   = this.grid.scanAxis;
            figsVars.thirdAxLabel  = this.grid.thirdAxis;

            figsVars.depthAxType   = uVars.figs.depthAxType;
            figsVars.scanAxType    = uVars.figs.scanAxType;
            
            figsVars.thirdAxCoor   = uVars.thirdPos;
            
            figsVars.validStruct = uVars.figs.validStruct;
            figsVars.extH        = uVars.figs.extH;
            figsVars.fonts       = uVars.figs.fonts;
            
            this.graphics.setUserVars(figsVars);
            
            this.extVars.figs = figsVars;
            this.measVars.figs = figsVars;
            
            % FileSystem parameters
            this.fileSystem.setUserVars(uVars.fileSystem); 
        end
        
        function vars = getVars(this)
            vars.uVars    = this.uVars;
            vars.extVars  = this.extVars;
            vars.measVars = this.measVars;
        end
        
        function configureScan(this)
            this.initResultsArrays2D();
            
            this.timeTable.scan = struct(); 
            this.graphics.setGraphicsScanVars();

            this.curPos = [this.grid.scanStart, this.grid.thirdPos];
            this.curPosIdx = [1, 1];
            
            %Config subObjects
            this.fileSystem.configFileSystem();
            this.fileSystem.saveVarsToDisk();
            
        end
        
        function configScope(this, chToSample)
            this.scope.ConfChannelFromVars(chToSample);
            this.scope.ConfAcq();
            this.scope.ConfTrigger();
            this.scope.ConfTimebase();
            this.scope.ConfWaveform();
        end
        
        % Results Functions
        function initResultsArrays2D(this)   
            this.res1D.rawData = zeros(1, this.grid.depthIdxLen);
            this.res1D.recon   = zeros(1, this.grid.depthIdxLen);
            this.res2D.rawData = zeros(this.grid.scanIdxLen, this.grid.depthIdxLen);
            this.res2D.recon   = zeros(this.grid.scanIdxLen, this.grid.depthIdxLen); 
            
            this.curPosIdx     = 0;
        end

        function putResToArrays(this, res)
            this.res1D.rawData = res.rawData;
            this.res1D.recon   = res.recon;
            this.res2D.rawData(this.curPosIdx(1), :) = res.rawData;
            this.res2D.recon(this.curPosIdx(1), :) = res.recon;
        end
        
        % Scan Function
        function res = scan(this)
            this.graphics.resetGraphics();
            this.graphics.resetPlots();
%             if this.measVars.stages.moveSecondStage
%                 this.stages.moveAbsAx(this.grid.thirdAxis, this.curPos(2))
%             end
            
            this.sf.startScanTime("timeTable2D", 'complete2DScan', this.curPos(2));
            this.sf.printStrModel("start2DScan", this.curPos(2), true);
                
            for i=1:this.grid.scanIdxLen
                this.updateCurPosAndIdx(i);
                this.sf.startScanTime("timeTable1D", 'singlePos', this.curPos(1));
                this.sf.printStrModel("start1DScan", this.curPos, true);

%                 this.fileSystem.updateFileSystem([this.curPos(1)]);

%                 this.stages.moveAbsAx(this.grid.scanAxis, this.curPos(1))

                this.sf.startScanTime("timeTable1D", 'acquire', this.curPos(1));
                resScope = this.scope.acquire();
                this.sf.stopScanTime("timeTable1D", 'acquire', this.curPos(1));

                this.sf.startScanTime("timeTable1D", 'reconstruction', this.curPos(1));
                res.rawData = resScope.data;
                res.recon = this.reconstructLine(res);
                this.putResToArrays(res);
                this.sf.stopScanTime("timeTable1D", 'reconstruction', this.curPos(1));

                this.sf.startScanTime("timeTable1D", 'plotting', this.curPos(1));
                this.plotResults();
                this.sf.stopScanTime("timeTable1D", 'plotting', this.curPos(1));

                this.sf.stopScanTime("timeTable1D", 'singlePos',  this.curPos(1));
            end

            this.fileSystem.saveResultsToDisk(this.res2D);
            this.fileSystem.closeFileSystem();
            res = this.res2D;
            
            this.sf.printStr(sprintf("Done 2D scan.\n"), true);
        end
        
        % Reconstruction Functions
        function recon = reconstructLine(this, res)
            avg = mean(res.rawData);
            dataNorm = res.rawData - avg;
            [envUp, envDown] = envelope(dataNorm);
            envTot = envUp + abs(envDown);
            Pz     = this.measVars.recon.pulseIntVal - cumsum(envTot);
            recon  = envTot./Pz;
        end
        
        % Position Functions
        function updateCurPosAndIdx(this, idx)
            % curPosIdx - always aligned to (R, 1st, 2nd)[#]
            % curPos    - always aligned to ( 1st, 2nd)[mm]
            this.curPosIdx(1) = idx;
            this.curPos(1)    = this.grid.scanVec(idx);
            if idx == 1
                this.curPos(2) = this.grid.thirdPos;
            end
            
%             this.graphics.updateCurPosAndIdx(this.getPos());
        end

        function pns = getPos(this)
            pns.curPos    = this.curPos;
            pns.curPosIdx = this.curPosIdx;
        end
        
        function updateSecondAxPos(this, thirdPos)
            this.curPos(2)     = thirdPos;
            this.grid.thirdPos = thirdPos;
        end
        
        function calc2DGrid(this)
            this.grid.scanVec     = this.grid.scanStart  : this.grid.scanStride  : this.grid.scanEnd;
            this.grid.scanLen     = abs(this.grid.scanStart - this.grid.scanEnd);
            this.grid.scanIdxLen  = length(this.grid.scanVec);
            this.grid.scanIdx     = 1:1:this.grid.scanIdxLen;        
            this.grid.scanCntr    = this.grid.scanVec - mean(this.grid.scanVec);
            this.grid.scanZero    = abs(this.grid.scanVec - this.grid.scanVec(1));
            
            this.grid.depthVec     = this.extVars.recon.depthVec;
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
                    switch this.grid.depthAxis
                        case 'Y'
                            this.grid.thirdAxis = 'Z';
                        case 'Z'
                            this.grid.thirdAxis = 'Y';
                    end
                case 'Y'
                    switch this.grid.depthAxis
                        case 'X'
                            this.grid.thirdAxis = 'Z';
                        case 'Z'
                            this.grid.thirdAxis = 'X';
                    end   
                case 'Z'
                    switch this.grid.depthAxis
                        case 'X'
                            this.grid.thirdAxis = 'Y';
                        case 'Y'
                            this.grid.thirdAxis = 'X';
                    end
            end 
        end

        % Graphics Functions
        function gH = getGraphicsHandle(this)
            gH = this.graphics;
        end

        function plotResults(this)
            this.graphics.setData(this.res2D);
            this.graphics.plotRawData();
            this.graphics.plotRecon();
        end
        
%         % Loaded data vars 
%         function setLoadedVars(this, uVars)
%             %TODO
%         end
%         
%         function setLoaded2DData(this, data, vars)
%             %TODO
%         end
    end
end

