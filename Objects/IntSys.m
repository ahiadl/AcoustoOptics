classdef IntSys < handle
    %OASCAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        aoS2D
        oa
        stages
        peUS
        
        figs
        fs
        sf
        
        vars
        grid
        
        curPos
        curPosIdx
        
        extStages
        
        res
        
        timeTable
    end
    
    methods (Static)
        function uVars = uVarsCreate()
            uVars.c              = 1500;
            uVars.imageWidth     = 30e-3;
            uVars.distFromCenter = 9.2e-2;
            
            uVars.aoS2D = scan2DAO.createUserVars();
            uVars.oa    = optoAcoustic.uVarsCreate();
            uVars.peus  = pulseEchoUS.uVarsCreate();

            uVars.grid.aoScanStart   = 0;
            uVars.grid.aoScanStride  = 0;
            uVars.grid.aoScanEnd     = 0;

            uVars.grid.oaScanStart  = 0;
            uVars.grid.oaScanStride = 0;
            uVars.grid.oaScanEnd    = 0;
            
            uVars.grid.aoScanAxis  = 'X';
            uVars.grid.aoDepthAxis = 'Y';
            uVars.grid.oaScanAxis  = 'Z';
            
            uVars.figs.intExt = 'int';
            uVars.figs.validStruct.corrected = true;
            uVars.figs.reopenFigures = false;
            
            uVars.figs.aoScanAxType  = 'Normal';
            uVars.figs.aoDepthAxType = 'Normal';
            uVars.figs.oaScanAxType  = 'Normal';
            
            uVars.figs.fonts = [];
            
            uVars.fileSystem.saveResults = false;
        end
    end
    
    methods
        function this = IntSys(handles)
         % handles should have 2 fields:
            % a. s2D - handle to initialized s2D object with initialized
            % stages and AO
            % b. stages - handle to the stages object
            % stages can be given solely, but if s2D handle is given its
            % stages handle must be given as well

            if (nargin > 0) && isfield(handles, 'stages')
                this.stages = handles.stages;
                this.extStages = true;
            else
                this.stages = stagesNew('PI'); %TODO: fix with new stages
                this.extStages = false; 
            end
            
            if (nargin > 0) && isfield(handles, 'aoS2D')
                %External AOS2D
                this.aoS2D = handles.aoS2D; %assuming external S2D will be given with already initialized stages object
            else
                %Internal AOS2D
                this.aoS2D = scan2DAO([], this.stages);
            end 
            
            if (nargin > 0) && isfield(handles, 'oa')
                this.oa = handles.oa;
            else
                this.oa = optoAcoustic([], this.stages);
            end 
            
            if (nargin > 0) && isfield(handles, 'peUS')
                this.peUS = handles.peUS;
            else
                this.peUS = pulseEchoUS(scope, this.stages);
            end 
            
            % this must be here so AO instrreset will not ruin stages
            % connection
            if ~this.extStages
                fprintf("S3D: Creating a stages object...\n");
                this.stages.connect();
                fprintf("S3D: Done connecting to stages.\n");
            end
            
            this.curPos = [0,0];
            this.sf = statsFunctions("INTSYS");
%             this.figs = IntSysGraphics();
%             hS3DFS = this.s2D.fileSystem;
%             this.fs = fileSystemIntSys(this, hS3DFS);
        end
        
        function setUserVars(this, uVars)
            %----------------------------
            % Acousto Optics S2D vars ---
            %----------------------------
            aoS2DuVars = scan2DAO.uVarsCreate();
            
            %General
            aoS2DuVars.grid.repeats         = 1; %disabled
            aoS2DuVars.general.externalScan = true;
            
            %Grid
            aoS2DuVars.grid.scanStart  = uVars.grid.aoScanStart; %this should be constant definition
            aoS2DuVars.grid.scanStride = uVars.grid.aoScanStride;
            aoS2DuVars.grid.scanEnd    = uVars.grid.aoScanEnd;   %this should be constant definition

            aoS2DuVars.grid.thirdPos   = uVars.grid.oaScanStart;
            
            aoS2DuVars.grid.scanAxis   = uVars.grid.aoScanAxis;
            aoS2DuVars.grid.depthAxis  = uVars.grid.aoDepthAxis;
            aoS2DuVars.grid.thirdAxis  = uVars.grid.oaScanAxis;
            
            %Stages
            aoS2DuVars.stages.moveSecondStage     = false;
            aoS2DuVars.stages.stagesOwnedByParent = true;
            
            % AO
            aoS2DuVars.ao.ao.c = uVars.c;
            
            aoS2DuVars.ao.ao.fSin             = uVars.aoS2D.ao.ao.fSin;              
            aoS2DuVars.ao.ao.fTrain           = uVars.aoS2D.ao.ao.fTrain;
            aoS2DuVars.ao.ao.cycInPulse       = uVars.aoS2D.ao.ao.cycInPulse;
            aoS2DuVars.ao.ao.channels         = uVars.aoS2D.ao.ao.channels;

            aoS2DuVars.ao.ao.fExtClk          = uVars.aoS2D.ao.ao.fExtClk; 
            aoS2DuVars.ao.ao.extClkDcyc       = uVars.aoS2D.ao.ao.extClkDcyc;
            aoS2DuVars.ao.ao.fSclk            = uVars.aoS2D.ao.ao.fSclk;

            aoS2DuVars.ao.ao.timeToSample     = uVars.aoS2D.ao.ao.timeToSample;
            aoS2DuVars.ao.ao.useQuant         = uVars.aoS2D.ao.ao.useQuant;
            aoS2DuVars.ao.ao.quantTime        = uVars.aoS2D.ao.ao.quantTime;

            aoS2DuVars.ao.ao.phantomDepth     = uVars.aoS2D.ao.ao.phantomDepth;
            aoS2DuVars.ao.ao.distFromPhantom  = uVars.aoS2D.ao.ao.distFromPhantom; %TODO: should be replaced with the value from reconstruction

            aoS2DuVars.ao.ao.envDC            = uVars.aoS2D.ao.ao.envDC;
            aoS2DuVars.ao.ao.envUS            = uVars.aoS2D.ao.ao.envUS;
            aoS2DuVars.ao.ao.envHar           = uVars.aoS2D.ao.ao.envHar;

            aoS2DuVars.ao.ao.fastAnalysis     = uVars.aoS2D.ao.ao.fastAnalysis;
            aoS2DuVars.ao.ao.useHadamard      = uVars.aoS2D.ao.ao.useHadamard;
            aoS2DuVars.ao.ao.useGPU           = uVars.aoS2D.ao.ao.useGPU;

            aoS2DuVars.ao.ao.useVirtualData   = uVars.aoS2D.ao.ao.useVirtualData;
            aoS2DuVars.ao.ao.limitByN         = uVars.aoS2D.ao.ao.limitByN;
            aoS2DuVars.ao.ao.N                = uVars.aoS2D.ao.ao.N;
            aoS2DuVars.ao.ao.dispTimeTable    = uVars.aoS2D.ao.ao.dispTimeTable;

            aoS2DuVars.ao.ao.IOPort           = uVars.aoS2D.ao.ao.IOPort;
            aoS2DuVars.ao.ao.IOLine           = uVars.aoS2D.ao.ao.IOLine;

            %Graphics
            
            % s2D
            aoS2DuVars.figs.depthIdx           = uVars.aoS2D.figs.depthIdx;
            
            aoS2DuVars.figs.intExt             = uVars.figs.intExt;
            aoS2DuVars.figs.hOwnerGraphics     = this.figs;
            aoS2DuVars.figs.validOwnerGraphics = false;
            aoS2DuVars.figs.useExtClims        = false; %TODO: support
            
            aoS2DuVars.figs.scanAxType         = uVars.figs.aoScanAxType;
            aoS2DuVars.figs.depthAxType        = uVars.figs.aoDepthAxType;
            
            aoS2DuVars.figs.reopenFigures      = uVars.figs.reopenFigures;
            
%             aoS2DuVars.figs.fonts              = uVars.figs.fonts;
            
            aoS2DuVars.figs.validStruct        = uVars.aoS2D.figs.validStruct;
            aoS2DuVars.figs.extH               = uVars.aoS2D.figs.extH;
            
            % AO
            aoS2DuVars.ao.figs.validStruct     = uVars.aoS2D.ao.figs.validStruct;
            aoS2DuVars.ao.figs.extH            = uVars.aoS2D.ao.figs.extH;
            
            aoS2DuVars.ao.figs.depthIdx        = uVars.aoS2D.ao.figs.depthIdx;
            aoS2DuVars.ao.figs.ch              = uVars.aoS2D.ao.figs.ch;
            aoS2DuVars.ao.figs.quant           = uVars.aoS2D.ao.figs.quant;

            aoS2DuVars.ao.figs.phiAxType       = uVars.figs.aoDepthAxType;
            aoS2DuVars.ao.figs.displayFullFFT  = uVars.aoS2D.ao.figs.displayFullFFT;
            aoS2DuVars.ao.figs.FFTenv          = uVars.aoS2D.ao.figs.FFTenv;
                        
            %Filesystem
            aoS2DuVars.fileSystem.extProject   = true;
            aoS2DuVars.fileSystem.dontSaveVars = true;
            aoS2DuVars.fileSystem.saveResults  = uVars.fileSystem.saveResults;
            
            aoS2DuVars.ao.fileSystem.saveRawData        = uVars.aoS2D.ao.fileSystem.saveRawData;
            aoS2DuVars.ao.fileSystem.saveNetSignal      = uVars.aoS2D.ao.fileSystem.saveNetSignal;
            aoS2DuVars.ao.fileSystem.saveDemultiplexed  = uVars.aoS2D.ao.fileSystem.saveDemultiplexed;
            aoS2DuVars.ao.fileSystem.saveReshapedSignal = uVars.aoS2D.ao.fileSystem.saveReshapedSignal;
            aoS2DuVars.ao.fileSystem.saveFFT            = uVars.aoS2D.ao.fileSystem.saveFFT;
            
            % Set Vars to S2D Acousto Optics
            this.aoS2D.setUserVars(aoS2DuVars)
            this.vars.aoS2D = this.aoS2D.getVars();
            
            %---------------------------
            %--- Opto Acoustics Vars ---
            %---------------------------
            OAuVars = this.oa.uVarsCreate();

            OAuVars.activeChannels            = uVars.oa.activeChannels; %Mask
            OAuVars.numOfSamples              = uVars.oa.numOfSamples;
            OAuVars.numOfAvg                  = uVars.oa.numOfAvg;
            OAuVars.numOfFrames               = uVars.oa.numOfFrames;

            OAuVars.c                         = uVars.c;
            OAuVars.fs                        = uVars.oa.fs;
            OAuVars.geometry                  = uVars.oa.geometry;
            OAuVars.imageWidth                = uVars.imageWidth; %[m]
            OAuVars.BPmode                    = uVars.oa.BPmode;
            
            OAuVars.figs.intExt               = uVars.figs.intExt;
            OAuVars.figs.reopenFigures        = uVars.figs.reopenFigures;
            OAuVars.figs.useExtClims          = false; %TODO: support
            
            OAuVars.figs.firstAxLabel         = uVars.grid.aoScanAxis;
            OAuVars.figs.secondAxLabel        = uVars.grid.aoDepthAxis;
            OAuVars.figs.scanAxLabel          = uVars.grid.oaScanAxis;
            
            OAuVars.figs.firstAxType          = uVars.figs.aoScanAxType;
            OAuVars.figs.secondAxType         = uVars.figs.aoDepthAxType;

            OAuVars.figs.scanAxCoor           = uVars.grid.oaScanStart;

            OAuVars.figs.validStruct          = uVars.oa.figs.validStruct;
            OAuVars.figs.extH                 = uVars.oa.figs.extH;

            OAuVars.fileSystem.saveResults    = uVars.fileSystem.saveResults;
%             OAuVars.fileSystem.saveFigs       = false; %TODO: not supported 

            OAuVars.fileSystem.extProject        = true;
            OAuVars.fileSystem.dontSaveVars      = true;
%             OAuVars.fileSystem.stackAllSubObjRes = false;
%             
%             OAuVars.fileSystem.useExtVarsPath    = false;
%             OAuVars.fileSystem.extVarsPath       = [];

            this.oa.setVars(OAuVars);
            this.vars.OA = this.oa.getVars();
            
            %------------------------------
            %--- Pulse-Echo Ultrasound Vars
            %------------------------------
            PEUSuVars = pulseEchoUS.uVarsCreate();

            PEUSuVars.distFromCenter = uVars.distFromCenter;
            PEUSuVars.imageWidth     = uVars.imageWidth;
            PEUSuVars.c              = uVars.c;
            PEUSuVars.usFreq         = uVars.peus.usFreq;

            % Scope
            PEUSuVars.signalChannel = uVars.peus.signalChannel;
            PEUSuVars.trigCh        = uVars.peus.trigCh;
            PEUSuVars.avg           = uVars.peus.avg;
            PEUSuVars.usRepRate     = uVars.peus.usRepRate;
            PEUSuVars.pulseIntVal   = uVars.peus.pulseIntVal;

            % Grid & Stages
            PEUSuVars.scanStart  = uVars.grid.aoScanStart;
            PEUSuVars.scanStride = uVars.grid.aoScanStride;
            PEUSuVars.scanEnd    = uVars.grid.aoScanEnd;

            PEUSuVars.thirdPos   = uVars.grid.oaScanStart;

            PEUSuVars.scanAxis   = uVars.grid.aoScanAxis;  
            PEUSuVars.depthAxis  = uVars.grid.aoDepthAxis;

            PEUSuVars.moveSecondStage     = false;
            PEUSuVars.stagesOwnedByParent = true;

            % Graphics
            PEUSuVars.figs.intExt        = uVars.figs.intExt;
            PEUSuVars.figs.useExtClims   = false; %TODO: support
            PEUSuVars.figs.reopenFigures = uVars.figs.reopenFigures; 

            PEUSuVars.figs.scanAxType    = uVars.figs.aoScanAxType;
            PEUSuVars.figs.depthAxType   = uVars.figs.aoDepthAxType;
            
            PEUSuVars.figs.validStruct   = uVars.peus.figs.validStruct;
            PEUSuVars.figs.extH          = uVars.peus.figs.extH;

            % FileSystem
            PEUSuVars.fileSystem.saveResults  = uVars.fileSystem.saveResults;

            PEUSuVars.fileSystem.extProject        = true;
            PEUSuVars.fileSystem.dontSaveVars      = true;

            this.peUS.setUserVars(PEUSuVars);
            this.vars.peUS = this.peUS.getVars();
            
            %-------------------------
            %--- Scan Grid Calculation
            %-------------------------
            this.grid.aoScanAxis  = uVars.grid.aoScanAxis;
            this.grid.oaScanAxis  = uVars.grid.oaScanAxis;
            this.grid.aoDepthAxis = uVars.grid.aoDepthAxis;
            
            this.grid.aoScanStart   = uVars.grid.aoScanStart;
            this.grid.aoScanStride  = uVars.grid.aoScanStride;
            this.grid.aoScanEnd     = uVars.grid.aoScanEnd;
            
            this.grid.oaScanStart  = uVars.grid.oaScanStart;
            this.grid.oaScanStride = uVars.grid.oaScanStride;
            this.grid.oaScanEnd    = uVars.grid.oaScanEnd;
            
            this.calcSysGrid();
            
            %--------------------------
            %--- Stats Function
            %--------------------------
            fprintf("INTSYS: 5. Setting strings models and variables\n");
            strings{1} = "start3DScan";
            models{1}  = sprintf("Start Scan for (%s) = (%s)", this.grid.oaScanAxis, "%.2f");
            strings{2} = "timeTable3D";
            models{2}  = sprintf("%s%s", this.grid.oaScanAxis, "%.2f");

            this.sf.setStringModels(strings, models);
        end
        
        function vars = getVars(this)
            vars = this.vars;
        end
        
        function configureScan(this)
            this.timeTable.scan = struct();
            this.initResultsArrays();

%             this.figs.setGraphicsScanVars();
%             this.figs.updateGraphicsConstruction();
            
            this.curPos    = this.grid.aoScanStart;
            this.curPosIdx = 1;
            
%             this.fs.configFileSystem(this.grid.secondAxis);
%             this.fs.saveVarsToDisk();
            this.aoS2D.configureScan();
            this.oa.configPeripherals();
            this.peUS.configureScan();
        end

        function updateCurPosAndIdx(this, oaScanIdx)
            this.curPosIdx = oaScanIdx;
            this.curPos(1) = this.grid.oaScanVec(oaScanIdx);
            
%             this.updateIntSysCurPosAndIdx(this.curPosIdx, this.curPos);
        end
        
        function initResultsArrays(this)
            aoScanAxLen    = this.grid.aoScanIdxLen;
            oaScanAxLen    = this.grid.oaScanIdxLen;
            aoDepthLen     = this.grid.aoDepthIdxLen;
            usDepthLen     = length(this.vars.peUS.measVars.recon.depthVec);
            
            oaChannels     = this.vars.OA.measVars.daq.numOfActiveChannels;
            oaNumOfSamples = this.vars.OA.measVars.daq.numOfSamples;
            reconSize      = this.vars.OA.measVars.bp.recon.n;
            
            this.res.OA.sigMat = zeros(oaScanAxLen, oaNumOfSamples, oaChannels);
            this.res.OA.recond = zeros(oaScanAxLen, reconSize, reconSize);
            
            this.res.AO.phi    = zeros(oaScanAxLen, aoScanAxLen, aoDepthLen);
            
            this.res.US.rawData = zeros(oaScanAxLen, aoScanAxLen, usDepthLen);
            this.res.US.recon   = zeros(oaScanAxLen, aoScanAxLen, usDepthLen);
        end
        
        function putResultsToArrays(this, resOA, resAO, resUS)
            this.res.OA.sigMat(this.curPosIdx, :, :)  = permute(resOA.sigMat, [3,1,2]); %permute might be needed
            this.res.OA.recon(this.curPosIdx, :, :)   = permute(resOA.recon, [3,1,2]); %permute might be needed
            
            this.res.AO.phi( this.curPosIdx, :, :)    = permute(resAO.phi, [2,1,3]);

            this.res.US.rawData(this.curPosIdx, :, :) = permute(resUS.rawData, [3,1,2]); %permute might be needed
            this.res.US.recon(this.curPosIdx, :, :)   = permute(resUS.recon, [3,1,2]); %permute might be needed
        end
        
        function res = scan(this)
            
            for i=1:this.grid.oaScanIdxLen
                this.updateCurPosAndIdx(i);
%                 this.sf.startScanTime("timeTable3D", 'singlePlane', this.curPos);
                this.sf.printStrModel("start3DScan", this.curPos, true);
                
%                 this.sf.startScanTime("timeTable3D", 'prepare', this.curPosIdx);
%                 this.fileSystem.updateFileSystem(this.curPos);
%                 this.s2D.updateSecondAxPos(this.curPos);
%                 this.sf.stopScanTime("timeTable3D", 'prepare', this.curPosIdx);
                
                this.stages.moveAbsAx(this.grid.oaScanAxis, this.curPos(1))
                
                fprintf("INTSYS: Please turn the OPO on\n");
                pause;
                resOA = this.oa.runOA();
                
                fprintf("INTSYS: Please turn the OPO off\n");
                pause;
                resAO = this.aoS2D.scan2D();
                
                fprintf("INTSYS: Please turn pulser on\n");
                pause;
                resUS = this.peUS.scan();
%                 resUS = this.peUS.runPEUS;
                
                this.putResultsToArrays(resOA, resAO, resUS);
                
%                 this.integrateCurrentPlane();
                
%                 this.sf.stopScanTime("timeTable3D", 'singlePlane', this.curPos);
            end
%             this.fs.saveResultsToDisk(this.res);
%             this.fs.closeFileSystem();
            res = this.res;
            this.sf.printStr(sprintf("Done Integrated System scan.\n"), true);
        end

        function calcSysGrid(this)
            this.grid.aoScanVec     = this.vars.aoS2D.grid.scanVec;
            this.grid.aoScanLen     = this.vars.aoS2D.grid.scanLen;
            this.grid.aoScanIdxLen  = this.vars.aoS2D.grid.scanIdxLen;
            this.grid.aoScanIdx     = this.vars.aoS2D.grid.scanIdx;
            this.grid.aoScanCntr    = this.vars.aoS2D.grid.scanCntr;
            this.grid.aoScanZero    = this.vars.aoS2D.grid.scanZero;
            
            this.grid.aoDepthVec     = this.vars.aoS2D.grid.depthVec;
            this.grid.aoDepthLen     = this.vars.aoS2D.grid.depthLen;
            this.grid.aoDepthIdxLen  = this.vars.aoS2D.grid.depthIdxLen;
            this.grid.aoDepthIdx     = this.vars.aoS2D.grid.depthIdx;
            this.grid.aoDepthCntr    = this.vars.aoS2D.grid.depthCntr;
            this.grid.aoDepthZero    = this.vars.aoS2D.grid.depthZero;
            
            this.grid.oaScanVec    = this.grid.oaScanStart : this.grid.oaScanStride : this.grid.oaScanEnd;
            this.grid.oaScanLen    = abs(this.grid.oaScanStart - this.grid.oaScanEnd);
            this.grid.oaScanIdxLen = length(this.grid.oaScanVec);
            this.grid.oaScanIdx    = 1:1:this.grid.oaScanIdxLen;
            this.grid.oaScanCntr   = this.grid.oaScanVec - min(this.grid.oaScanVec);
            this.grid.oaScanZero   = abs(this.grid.oaScanVec - this.grid.oaScanVec(1));
            
            switch this.grid.aoScanAxis
                case 'X'
                    switch this.grid.oaScanAxis
                        case 'Y'
                            this.grid.mainPlane = 'XZ';
                        case 'Z'
                            this.grid.mainPlane = 'XY';
                    end
                case 'Y'
                    switch this.grid.oaScanAxis
                        case 'X'
                            this.grid.mainPlane = 'YZ';
                        case 'Z'
                            this.grid.mainPlane = 'XY';
                    end   
                case 'Z'
                    switch this.grid.oaScanAxis
                        case 'X'
                            this.grid.mainPlane = 'YZ';
                        case 'Y'
                            this.grid.mainPlane = 'XZ';
                    end
            end 
        end
        
    end
end

%         function integrateCurrentPlane(this)
%             
%             
%         end
        
%         function turnOpoOn(this)
%             questdlg('Please turn on OPO on.',...
%                      'Operator Action Needed', 'OK');
%         end
%         
%         function turnOpoOff(this)
%             questdlg('Please turn OPO off and then Connect US transducer to amplifier.',...
%                      'Operator Action Needed', 'OK');
%         end
%         
%         function connectUSPulser(this)
%             questdlg('Please connect US pulser to US transducer.',...
%                      'Operator Action Needed', 'OK');
%         end
