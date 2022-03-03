classdef contScan < handle
    %SCAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        stages
        daq
        vars
        res
        plots
        data
    end
    
    methods (Static)
        function uVars = createUserVars()
            uVars.numOfCh      = 1;
            uVars.tVec         = [];
            uVars.laserRepRate = 30; %[Hz]
                        
            uVars.binSize  = 1;
            uVars.binType  = 'mid'; % 'mid', 'post'
            uVars.adaptBin = false; % set optimized binning according to velocity
            
            uVars.scanType      = 'cont'; % 'cont', 'disc', 'pos'
            uVars.scanDirection = 'bi';   %'bi', 'uni'
            uVars.stabDist      = 0.1; %[mm]
            uVars.infSinglePos  = false;
            uVars.limSinglePos  = 1;
            uVars.allAtOnce     = false;
            
            uVars.axScan   = 'Y';  %'X', 'Y', 'Z'
            uVars.axDisc1  = 'Z';  %'X', 'Y', 'Z'
            uVars.axDisc2  = 'X';  %'X', 'Y', 'Z'
            uVars.spanType = 'center'; %'center', 'limits'
            
            uVars.axScanSpan  = 1;%[mm]
            uVars.axDisc1Span = 0;%[mm]
            uVars.axDisc2Span = 0;%[mm]
            
            uVars.axScanStride  = 0.1; %[mm]
            uVars.axDisc1Stride = 0.1;%[mm]
            uVars.axDisc2Stride = 0.1;%[mm]
            
            uVars.axScanRef  = [];%[mm]
            uVars.axDisc1Ref = [];%[mm]
            uVars.axDisc2Ref = [];%[mm]
            
            %External Functions
            uVars.doPP        = true;
            uVars.filter      = [];
            uVars.postProc    = [];
            uVars.initPlots   = [];
            uVars.plotResults = [];
            uVars.userAux     = []; %external user data to support user PostProcessing
            
            uVars.plotMode  = 'mid'; % 'max', 'mid', 'user'
            uVars.chToPlot  = 1;
            uVars.posToPlot = 0;
            
            uVars.saveData = false;
            uVars.dir = '';
            uVars.filename = '';
            
            uVars.tg.chatID  = [];
            uVars.tg.text    = false;
            uVars.tg.figure  = false;
            uVars.tg.rep     = 1;
        end
    end
    methods
        function this = contScan(daq, stages)
            this.stages = stages;
            this.daq = daq;
            this.vars.trig.pin = 1;
            this.vars.trig.mode = 'singlePosition';
            this.vars.trig.params = [];
            this.vars.trig.pol = 1;
        end
        
        function curVars = getVars(this)
            curVars = this.vars;
        end
        
        function status = setVars(this, uVars)
            status = true;
            
            this.vars.uVars         = uVars;
            this.vars.numOfCh       = uVars.numOfCh;
            this.vars.tVec          = uVars.tVec;
            this.vars.delay         = this.vars.tVec(1);
            this.vars.fs            = 1/abs(this.vars.tVec(2)-this.vars.tVec(1));
            this.vars.numOfSamples  = length(this.vars.tVec);
            this.vars.laserRepRate  = uVars.laserRepRate;
            
            this.vars.binSize       = uVars.binSize;
            this.vars.binType       = uVars.binType;
            this.vars.adaptBin      = uVars.adaptBin; 
            
            this.vars.scanType      = uVars.scanType; %'cont', 'disc', 'pos'
            this.vars.scanDirection = uVars.scanDirection; % 'uni' / 'bi'
            this.vars.stabDist      = uVars.stabDist; % a distance for stages velocity to stabilize at the beggining of each line.
            this.vars.infSinglePos  = uVars.infSinglePos;
            this.vars.limSinglePos  = uVars.limSinglePos;
            this.vars.allAtOnce     = uVars.allAtOnce; % in case of single point should measure all repetitions at once or measure each trigger and display
            this.vars.stabTime      = uVars.stabTime;
            
            this.vars.axScan   = uVars.axScan;   %'X', 'Y', 'Z'
            this.vars.axDisc1  = uVars.axDisc1;  %'X', 'Y', 'Z'
            this.vars.axDisc2  = uVars.axDisc2;  %'X', 'Y', 'Z'
            this.vars.spanType = uVars.spanType; %'center', 'limits'
            
            this.vars.axScanSpan    = uVars.axScanSpan;
            this.vars.axDisc1Span   = uVars.axDisc1Span;
            this.vars.axDisc2Span   = uVars.axDisc2Span;
            
            this.vars.axScanStride  = uVars.axScanStride;
            this.vars.axDisc1Stride = uVars.axDisc1Stride;
            this.vars.axDisc2Stride = uVars.axDisc2Stride;
            
            this.vars.axScanRef  = uVars.axScanRef;
            this.vars.axDisc1Ref = uVars.axDisc1Ref;
            this.vars.axDisc2Ref = uVars.axDisc2Ref;
            
            switch this.vars.scanType
                case 'cont'
                    this.vars.axScanRes = this.vars.axScanStride/this.vars.binSize;
                    this.vars.scanVel   = this.vars.axScanRes * this.vars.laserRepRate;
                    this.vars.maxVel    = this.stages.getMaxVelocityAx(this.vars.axScan);
                    
                    if this.vars.scanVel > this.vars.maxVel
                        if ~this.vars.adaptBin
                            error("CS: Error: Requested scan velocity is too high.");
                        else
                            this.vars.scanVel = this.vars.maxVel;
                            this.vars.axScanRes = this.vars.scanVel/this.vars.laserRepRate;
                            this.vars.binSize   = ceil(this.vars.axScanStride/this.vars.axScanRes);
                            % we re calculate these parameters since floor might
                            % change the max calculated bin size
                            this.vars.axScanRes = this.vars.axScanStride/this.vars.binSize;
                            this.vars.scanVel   = this.vars.axScanRes * this.vars.laserRepRate;
                            fprintf("CS: Scan velocity too high.\n Bin size was increased to %d, new velocity: %.2f [mm/s]",...
                                this.vars.binSize, this.vars.scanVel);
                        end
                    end

                    switch this.vars.binType
                        case 'mid'
                            preSpan  = this.vars.axScanStride /2;
                            postSpan = this.vars.axScanStride /2; 

                        case 'post'
                            preSpan  = 0;
                            postSpan = this.vars.axScanStride;
                    end
    
                     this.vars.stabTime = max(this.vars.stabDist/this.vars.scanVel-1e-1,0);
                case 'disc'
                    this.vars.axScanRes = this.vars.axScanStride;
                    preSpan  = 0;
                    postSpan = 0;
                case 'pos'
                    this.vars.spanType = 'point';
                    this.vars.scanDirection = 'point';
                    preSpan  = 0;
                    postSpan = 0; 
            end
            
            switch this.vars.spanType
                case 'center'
%                     this.vars.axScanCenter  = uVars.axScanCenter;
%                     this.vars.axDisc1Center = uVars.axDisc1Center;
%                     this.vars.axDisc2Center = uVars.axDisc2Center;

                    this.vars.scanVec = ...
                        [-(this.vars.axScanSpan + preSpan) : this.vars.axScanRes     : (this.vars.axScanSpan + postSpan - this.vars.axScanRes)]     + this.vars.axScanRef;
                    this.vars.disc1Vec = ...
                        [-this.vars.axDisc1Span: this.vars.axDisc1Stride : (this.vars.axDisc1Span - this.vars.axDisc1Stride)] + this.vars.axDisc1Ref;
                    this.vars.disc2Vec = ...
                        [-this.vars.axDisc2Span: this.vars.axDisc2Stride : (this.vars.axDisc2Span - this.vars.axDisc2Stride)] + this.vars.axDisc2Ref;
                    
                    this.vars.disc1Vec = ...
                        [((-this.vars.axDisc1Span: this.vars.axDisc1Stride : (this.vars.axDisc1Span - this.vars.axDisc1Stride))-3.5),...
                        ((-this.vars.axDisc1Span: this.vars.axDisc1Stride : (this.vars.axDisc1Span - this.vars.axDisc1Stride))),...
                        ((-this.vars.axDisc1Span: this.vars.axDisc1Stride : (this.vars.axDisc1Span - this.vars.axDisc1Stride))+3.5)]...
                         + this.vars.axDisc1Ref;
                        
                case 'limits'
                    this.vars.scanVec = ...
                        [-preSpan : sign(this.vars.axScanSpan)*this.vars.axScanRes      : this.vars.axScanSpan + postSpan]  + this.vars.axScanRef;%
                    this.vars.disc1Vec = ...
                        [0 : sign(this.vars.axDisc1Span)*this.vars.axDisc1Stride : this.vars.axDisc1Span] + this.vars.axDisc1Ref;
                    this.vars.disc2Vec = ...
                        [0 : sign(this.vars.axDisc2Span)*this.vars.axDisc2Stride : this.vars.axDisc2Span] + this.vars.axDisc2Ref;
                
                case 'point'
                    if this.vars.infSinglePos
                        this.vars.allAtOnce = false;
                    end
                    if this.vars.allAtOnce
                        this.vars.scanVec  = ones(1, this.vars.binSize)*this.vars.axScanRef;
                        this.vars.disc1Vec = this.vars.axDisc1Ref;
                    else
                        this.vars.scanVec  = ones(1, 1)*this.vars.axScanRef;
                        this.vars.disc1Vec = ones(1, this.vars.binSize)*this.vars.axDisc1Ref;
                    end
                    this.vars.disc2Vec = this.vars.axDisc2Ref;
            end
            
                       
            if this.vars.axScanSpan == 0 && ~strcmp(this.vars.scanType, 'pos')
                this.vars.scanVec = this.vars.axScanRef;
            end
            if this.vars.axDisc1Span == 0
                this.vars.disc1Vec = this.vars.axDisc1Ref;
            end
            if this.vars.axDisc2Span == 0
                this.vars.disc2Vec = this.vars.axDisc2Ref;
            end
            this.vars.scanSize  = [length(this.vars.scanVec),...
                                        length(this.vars.disc1Vec),...
                                        length(this.vars.disc2Vec)];
            
            if strcmp(this.vars.scanType, 'pos')
                this.vars.daqAcqNum = this.vars.binSize;
                this.vars.scanVecBin = this.vars.scanVec;
                
                if this.vars.allAtOnce
                    this.vars.scanSizeBin =  [1,...
                                              1,...
                                              length(this.vars.disc2Vec)];
                else              
                    this.vars.scanSizeBin =  [1,...
                                              this.vars.binSize,...
                                              length(this.vars.disc2Vec)];
                end
                                      
            elseif strcmp(this.vars.scanType, 'disc') 
                this.vars.daqAcqNum = this.vars.binSize;
                this.vars.scanVecBin = this.vars.scanVec;
                
                this.vars.scanSizeBin =  [length(this.vars.scanVecBin),...
                                          length(this.vars.disc1Vec),...
                                          length(this.vars.disc2Vec)];
            else
                this.vars.daqAcqNum = this.vars.scanSize(1);
                this.vars.scanVecBin = this.vars.scanVec(1:this.vars.binSize:end)- preSpan;
                
                this.vars.scanSizeBin =  [length(this.vars.scanVecBin),...
                                          length(this.vars.disc1Vec),...
                                          length(this.vars.disc2Vec)];           
            end
        
            switch this.vars.scanDirection
                case 'uni'
                    this.vars.scanStartVec(1:this.vars.scanSize(2)) = this.vars.scanVec(1)   - sign(this.vars.axScanSpan)*this.vars.stabDist;
                    this.vars.scanEndVec(1:this.vars.scanSize(2)) = this.vars.scanVec(end) + sign(this.vars.axScanSpan)*this.vars.stabDist;
                    this.vars.trigVec(1:this.vars.scanSize(2)) = this.vars.scanVec(1);
                    this.vars.polarVec(1:this.vars.scanSize(2)) = (sign(this.vars.axScanSpan)+1)/2;
                case 'bi'
%                     this.vars.scanStartVec = zeros(1,this.vars.scanSize(2));
                    this.vars.scanStartVec(1:2:this.vars.scanSize(2)) = this.vars.scanVec(1)   - sign(this.vars.axScanSpan)*this.vars.stabDist;
                    this.vars.scanStartVec(2:2:this.vars.scanSize(2)) = this.vars.scanVec(end) + sign(this.vars.axScanSpan)*this.vars.stabDist;

%                     this.vars.scanEndVec = zeros(1,this.vars.scanSize(2));
                    this.vars.scanEndVec(1:2:this.vars.scanSize(2)) = this.vars.scanVec(end) + sign(this.vars.axScanSpan)*this.vars.stabDist;
                    this.vars.scanEndVec(2:2:this.vars.scanSize(2)) = this.vars.scanVec(1)   - sign(this.vars.axScanSpan)*this.vars.stabDist;

%                     this.vars.trigVec = zeros(1,this.vars.scanSize(2));
                    this.vars.trigVec(1:2:this.vars.scanSize(2)) = this.vars.scanVec(1);
                    this.vars.trigVec(2:2:this.vars.scanSize(2)) = this.vars.scanVec(end);

                    %This code is to set the polarity of the trrigger signal 
                    %in the bi-directional scan. Since the scan might be in 
                    %either positive or negative direction, while trigger is 
                    %active only if position is larger or equal to the defined 
                    %trigger position, the trigger signal polarity should be inverted.
                    
%                     this.vars.polarVec = zeros(1,this.vars.scanSize(2));
                    this.vars.polarVec(1:2:this.vars.scanSize(2)) = (sign(this.vars.axScanSpan)+1)/2;
                    this.vars.polarVec(2:2:this.vars.scanSize(2)) = abs((sign(this.vars.axScanSpan)-1)/2);
                case 'point'
                    
            end

            %User Functions
            this.vars.doPP            = uVars.doPP;
            this.vars.filter          = uVars.filter;
            this.vars.userPP          = uVars.postProc;
            this.vars.initUserPlots   = uVars.initPlots;
            this.vars.userPlotResults = uVars.plotResults;
            this.vars.userAux         = uVars.userAux;
            
            % Plots parameters
            this.vars.plotMode  = uVars.plotMode;
            this.vars.chToPlot    = uVars.chToPlot;
            this.vars.posToPlot  = uVars.posToPlot;
            [~, this.vars.posIdxToPlot] = min(abs(this.vars.scanVec - this.vars.posToPlot));
            
            % Save Data
            this.vars.saveData = uVars.saveData;
            this.vars.dir = uVars.dir;
            this.vars.userFilename = uVars.filename;
            if this.vars.saveData && ~exist(this.vars.dir, 'dir')
                fprintf("CS: Error: requested directory does not exist saving to current dir.\n");
                this.vars.dir = '.';
            end
                        
            % Telegram Messages
            this.vars.tg.chatID  = uVars.tg.chatID;
            this.vars.tg.text    = uVars.tg.text;
            this.vars.tg.figure  = uVars.tg.figure;
            this.vars.tg.rep     = uVars.tg.rep;
            
            fprintf("CS: Done setting Vars.\n");                        
        end
        
        function resCs = scan(this)
            switch this.vars.scanType
                case 'cont'
                    resCs = this.scanCont();
                case 'disc'
                    resCs = this.scanDisc();
                case 'pos'
                    resCs = this.scanPos();
            end
        end
        
        function resCs = scanDisc(this)
            this.initPlots();
            this.data= [];
            this.data.resMat = zeros([this.vars.scanSizeBin, this.vars.numOfCh, this.vars.numOfSamples]); 
%             vel = this.stages.getMaxVelocityAx(this.vars.axScan);
%             this.stages.setVelocity(this.vars.axScan, vel);
            this.stages.moveAbsAx(this.vars.axScan, this.vars.scanStartVec(1));
%             this.stages.turnOnIOAll(this.vars.trig.pin);
            this.vars.averageTime = 0;
            for i=1:this.vars.scanSize(3)
                this.stages.moveAbsAx(this.vars.axDisc2, this.vars.disc2Vec(i))
                this.vars.idxs(2) = i;
                for j = 1:this.vars.scanSize(2)
                    this.stages.moveAbsAx(this.vars.axDisc1, this.vars.disc1Vec(j))
                    this.vars.idxs(1) = j;
                    T = tic;
                    for k = 1:this.vars.scanSize(1)                   
                       
                        if strcmp(this.vars.scanDirection, 'bi') && ~mod(j,2)
                            curPos = this.vars.scanVec(end-(k-1));
                        else
                            curPos = this.vars.scanVec(k);
                        end
                        this.stages.moveAbsAx(this.vars.axScan, curPos)
                        pause(this.vars.stabTime)
                        % sigMat should return as [time x Ch x trigger(pos)]
                        sigMat = this.daq.acquire();
                        if ~mod(j,2) && strcmp(this.vars.scanDirection, 'bi')
                            sigMat = flip(sigMat,3);
                        end

                        sigMatBin = mean(reshape(sigMat,...
                                         this.vars.numOfSamples,...
                                         this.vars.numOfCh,...
                                         this.vars.binSize), 3);             

                        this.data.curData = permute(sigMatBin, [3,4,5,2,1]);
    %                     resMat[scan disc1 disc2 ch t]
                        this.data.resMat(k, j, i, :, :) = this.data.curData;
                        this.postProc();
                        this.plotResults();
                    end
                    this.data.curData = this.data.resMat(:, j, i, :, :);
                    this.postProc();
                    this.plotResults();
                    T = toc(T);
                    this.printProgress(T);
                end
            end
%             this.stages.setVelocity(this.vars.axScan, vel);
            resCs = this.data.resMat;
            if this.vars.saveData
                this.saveData();
                if this.vars.tg.text
                    tgprintf(this.vars.tg.chatID, "Data Saved")
                end 
            end
            if this.vars.tg.text
                tgprintf(this.vars.tg.chatID, "Scan Ended")
            end
        end
        
        function resCs = scanPos(this)
            this.initPlots();
            this.data= [];
             
%             vel = this.stages.getMaxVelocityAx(this.vars.axScan);
%             this.stages.setVelocity(this.vars.axScan, vel);
            % Move stages to requested Position
            this.stages.moveAbsAx(this.vars.axScan, this.vars.scanVec(1))
            this.stages.moveAbsAx(this.vars.axDisc1, this.vars.disc1Vec(1))
            this.stages.moveAbsAx(this.vars.axDisc2, this.vars.disc2Vec(1))
%             this.stages.turnOnIOAll(this.vars.trig.pin);
            this.vars.idxs(2) = 1;
            this.vars.averageTime = 0;
            if this.vars.infSinglePos
                this.vars.idxs(1) = 0;
                while true
                    T = tic;
                    % sigMat should return as [time x Ch x trigger(pos)]
                    % for single pos trigeer = binning
                    sigMat = this.daq.acquire(); 
                    sigMatBin = mean(sigMat, 3);             
                                 
                    this.data.curData = permute(sigMatBin, [3,4,5,2,1]);
                    this.postProc();
                    this.plotResults();
                    T = toc(T);
                    this.printProgress(T); 
                end
            else
                this.data.resMat = zeros([this.vars.scanSizeBin, this.vars.numOfCh, this.vars.numOfSamples]);
                for j = 1:this.vars.scanSize(2)
                    this.vars.idxs(1) = j;
                    T = tic;
                    sigMat = this.daq.acquire(); 
                    sigMatBin = mean(reshape(sigMat,...
                                     this.vars.numOfSamples,...
                                     this.vars.numOfCh,...
                                     this.vars.binSize,...
                                     this.vars.scanSizeBin(1)), 3);  
                    this.data.curData = permute(sigMatBin, [4,3,5,2,1]);
%                     resMat[scan disc1 disc2 ch t]
                    this.data.resMat(:, j, 1, :, :) = this.data.curData;
                    this.postProc();
                    this.plotResults();
                    T = toc(T);
                    this.printProgress(T); 
                end
                resCs = this.data.resMat;
                if this.vars.saveData
                    this.saveData();
                    if this.vars.tg.text
                        tgprintf(this.vars.tg.chatID, "Data Saved")
                    end 
                end
                if this.vars.tg.text
                    tgprintf(this.vars.tg.chatID, "Scan Ended")
                end
            end
        end

        function resCs = scanCont(this)
            this.initPlots();
            this.data= [];
            this.data.resMat = zeros([this.vars.scanSizeBin, this.vars.numOfCh, this.vars.numOfSamples]); 
            vel = this.stages.getMaxVelocityAx(this.vars.axScan);
            this.stages.setVelocity(this.vars.axScan, vel);
            this.stages.moveAbsAx(this.vars.axScan, this.vars.scanStartVec(1));
            this.stages.setVelocity(this.vars.axScan, this.vars.scanVel);
            this.stages.turnOnIOAll(this.vars.trig.pin);
            this.vars.averageTime = 0;
            for i=1:this.vars.scanSize(3)
                this.stages.moveAbsAx(this.vars.axDisc2, this.vars.disc2Vec(i))
                this.vars.idxs(2) = i;
                for j = 1:this.vars.scanSize(2)
                    this.vars.idxs(1) = j;
                    T = tic;
                    this.stages.moveAbsAx(this.vars.axDisc1, this.vars.disc1Vec(j))
                    if strcmp(this.vars.scanDirection, 'uni')
                        this.stages.moveAbsAx(this.vars.axScan, this.vars.scanStartVec(j));
                    end
                    trigParams.dest    = this.vars.scanEndVec(j);
                    trigParams.trigPos = this.vars.trigVec(j);
                    this.stages.setTrigger(this.vars.axScan, this.vars.trig.pin,...
                                           this.vars.trig.mode, trigParams,...
                                           ~this.vars.polarVec(j));
%                     this.stages.moveNBAbsAx(this.vars.axScan, this.vars.scanEndVec(j));
                    this.stages.startRoutine(this.vars.axScan);
                    % sigMat should return as [time x Ch x trigger(pos)]
                    pause(this.vars.stabTime);
                    sigMat = this.daq.acquire();
                    if ~mod(j,2) && strcmp(this.vars.scanDirection, 'bi')
                        sigMat = flip(sigMat,3);
                    end

                    sigMatBin = mean(reshape(sigMat,...
                                     this.vars.numOfSamples,...
                                     this.vars.numOfCh,...
                                     this.vars.binSize,...
                                     this.vars.scanSizeBin(1)), 3);             
                                 
                    this.data.curData = permute(sigMatBin, [4,3,5,2,1]);
%                     resMat[scan disc1 disc2 ch t]
                    this.data.resMat(:, j, i, :, :) = this.data.curData; 
                    this.postProc();
                    this.plotResults();
                    this.stages.blockWhileMoving(this.vars.axScan);
                    T = toc(T);
                    this.printProgress(T); 
                end
                this.stages.blockWhileMoving(this.vars.axDisc1);
                this.stages.blockWhileMoving(this.vars.axDisc2);
            end
            this.stages.setVelocity(this.vars.axScan, vel);
            resCs = this.data.resMat;
            if this.vars.saveData
                this.saveData();
                if this.vars.tg.text
                    tgprintf(this.vars.tg.chatID, "Data Saved")
                end 
                fprintf("CS: Data Saved.\n");
            end
            if this.vars.tg.text
                tgprintf(this.vars.tg.chatID, "Scan Ended")
            end
            fprintf("CS: Scan Ended.\n");
        end
        
        function postProc(this)
            i = this.vars.idxs(2);
            j = this.vars.idxs(1);
            
            if j==1
                this.data.p2p = zeros([this.vars.scanSizeBin, this.vars.numOfCh]);
            end
            
            if this.vars.doPP
                %AC Coupling
                this.data.ppData = this.data.curData - mean(this.data.curData,5);

                %Filtering
                if isa(this.vars.filter, 'function_handle')
                    this.data.ppData = this.vars.filter(this.data.ppData);
                end
            else
                this.data.ppData = this.data.curData;
            end
            
            %peak2peak
            this.data.p2p(:,j,i,:) = peak2peak(this.data.ppData,5);
            
            %User PostProcessing
            if isa (this.vars.userPP, 'function_handle')
               this.data.uPP = this.vars.userPP(this.vars.userAux, this.data, this.vars); 
            end
        end
        
        function initPlots(this)
           this.plots.hFig = figure();
            
           this.plots.ax1 = subplot(2,2,1);
           this.plots.ax2 = subplot(2,2,2);
           this.plots.ax3 = subplot(2,1,2);
           
           % Plot Channel vs. Time sinogram for Chosen Position
           if this.vars.numOfCh >1
               this.plots.hP1 = imagesc(this.plots.ax1,...
                                        'XData', this.vars.tVec*1e6,...
                                        'YData', 1:this.vars.numOfCh,...
                                        'CData', zeros(this.vars.numOfSamples, this.vars.numOfCh));
               ylabel(this.plots.ax1, 'Channel[#]')
               this.plots.hTitP1 = title(this.plots.ax1, "Channel Vs. Time at Position");
               colorbar(this.plots.ax1);
           else
               this.plots.hP1 = plot(this.plots.ax1, this.vars.tVec*1e6, zeros(1,this.vars.numOfSamples));
               ylabel(this.plots.ax1, 'Amp [mV]')
               this.plots.hTitP1 = title(this.plots.ax1, "Channel Vs. Time at Position");
           end
           axis(this.plots.ax1, 'tight')
           xlabel(this.plots.ax1, 't[\mus]')
           
           % Plot Position vs. Time sinogram for Chosen Channel
           this.plots.hP2 = imagesc(this.plots.ax2,...
                                    'XData', this.vars.tVec*1e6,...
                                    'YData', this.vars.scanVecBin,...
                                    'CData', zeros(this.vars.scanSizeBin(1),this.vars.numOfSamples));
           this.plots.hTitP2 = title(this.plots.ax2, "Position Vs. Time on Channel");
           axis(this.plots.ax2, 'tight')
           xlabel(this.plots.ax2, 't[\mus]')
           ylabel(this.plots.ax2, 'Scan Pos [mm]')
           colorbar(this.plots.ax2);
           
           % Plot peak2peak Image
           if this.vars.scanSize(2) > 1
               this.plots.hP3 = imagesc(this.plots.ax3,...
                                        'XData', this.vars.scanVecBin,...
                                        'YData', this.vars.disc1Vec,...
                                        'CData', zeros(this.vars.scanSizeBin(1), this.vars.scanSizeBin(2)));
               ylabel(this.plots.ax3, 'Discrete Axis Pos [mm]')
               axis(this.plots.ax3, 'equal');
               colorbar(this.plots.ax3);
           else
               this.plots.hP3 = plot(this.vars.scanVecBin, zeros(1, this.vars.scanSizeBin(1)));
               ylabel(this.plots.ax3, 'P2P Amp [mV]')
           end
           this.plots.hTitP3 = title(this.plots.ax3, "Peak2Peak Map");
           axis(this.plots.ax3, 'tight')
           xlabel(this.plots.ax3, 'Cont Axis Pos [mm]');
           
           if isa (this.vars.initUserPlots, 'function_handle')
              this.plots.userPlots = this.vars.initUserPlots(this.vars);
           end
        end
        
        function plotResults(this)
            %data - [scan, 1, 1, ch, time]
            i = this.vars.idxs(2);
            j = this.vars.idxs(1);
            switch this.vars.plotMode
                case 'max'
                    [~, idxCh]  = max(max(this.data.p2p, 1), 2);
                    [~, idxPos] = max(max(this.data.p2p, 2), 1);
                case 'mid'
                    idxCh  = ceil(this.vars.numOfCh/2);
                    idxPos = ceil(this.vars.scanSizeBin(1)/2);
                case 'user'
                    idxCh  = this.vars.chToPlot;
                    idxPos = this.vars.posIdxToPlot;
            end
            
            if strcmp(this.vars.scanType, 'disc') || strcmp(this.vars.scanType, 'pos')
                idxPos = 1;
            end
            
            % Plot Channel vs. Time sinogram for Chosen Position
            if this.vars.numOfCh >1
               set(this.plots.hP1, 'CData', squeeze(this.data.ppData(idxPos,:,:,:,:))*1e3)
            else
               set(this.plots.hP1, 'YData', squeeze(this.data.ppData(idxPos,:,:,:,:))) 
            end
            set(this.plots.hTitP1, 'String', sprintf("Channel Vs. Time on Position: %d", idxPos));

            % Plot Position vs. Time sinogram for Chosen Channel
            if strcmp(this.vars.scanType, 'cont')
                pData = this.data.ppData(:,:,:,idxCh,:)*1e3;
            else
                pData = this.data.resMat(:,j,i,idxCh,:);
            end
            
            set(this.plots.hP2, 'CData', squeeze(pData));
            set(this.plots.hTitP2, 'String', sprintf("Position Vs. Time on Channel: %d", idxCh));
            
            % Plot peak2peak Image for selected ch
            if this.vars.scanSize(2) >1
                set(this.plots.hP3, 'CData', squeeze(this.data.p2p(:,:,i,idxCh))')
            else
                set(this.plots.hP3, 'YData', squeeze(this.data.p2p(:,:,i,idxCh)))
            end
            
            if isa (this.vars.userPlotResults, 'function_handle')
              this.vars.userPlotResults(this.plots.userPlots, this.data.uPP, this.vars);
            end
            
            if this.vars.tg.figure && ~mod(this.vars.idxs(1)-1, this.vars.tg.rep)
                tgprint(this.vars.tg.chatID, this.plots.hFig, 'photo');
            end
               
            drawnow();
        end
        
        function printProgress(this, T)
           i = this.vars.idxs(2);
           j = this.vars.idxs(1);
           totalNumOfLines = this.vars.scanSize(2)*this.vars.scanSize(3);
           curNumOfLines = (i-1)*this.vars.scanSize(2)+j;
           percent = curNumOfLines/totalNumOfLines * 100;
           this.vars.averageTime = (this.vars.averageTime *(curNumOfLines-1) +T )/curNumOfLines;
           totalTime = totalNumOfLines * this.vars.averageTime / 60;
           timeToNow = curNumOfLines * this.vars.averageTime / 60;
           str = sprintf("Done line %d/%d (%.2f [%%%%]). Progress: %.2f/%.2f mins.\n",...
                            curNumOfLines, totalNumOfLines, percent, timeToNow, totalTime);
           fprintf(str)
           if this.vars.tg.text && ~mod(this.vars.idxs(1)-1, this.vars.tg.rep)
               tgprintf(this.vars.tg.chatID, str);
           end 
        end
        
        function saveData(this)
            resCs = this.data.resMat;
            csVars = this.getVars();
            timeStamp = strrep(datestr(datetime('now')),':','-');
            this.vars.filename = sprintf("%s/%s-%s.mat", this.vars.dir, timeStamp, this.vars.userFilename);
            save(this.vars.filename, 'resCs', 'csVars', '-v7.3');
        end
    end
end

