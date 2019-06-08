classdef scan2d < scanObj
    %SCAN2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function this = scan2d()
            this.strings.scan = "Done Scan for (X,Y,Q) = (%8.3f, %8.3f, %d)\n";
            this.strings.timeTable = "X%8.3fY%8.3f%d";
        end
        
        function setScanUserVars(this, uVars)
            this.scan.speckleTime = uVars.speckleTime;
            this.scan.useQuant    = uVars.useQuant;
            
            if this.scan.useQuant
                this.scan.numOfQuant   = ceil(uVars.timeFrame / uVars.speckleTime);
                this.scan.timeToSample = ones(1, this.scan.numOfFrames) * this.scan.speckleTime;
            else 
                this.scan.numOfQuant   = ones(1, this.scan.numOfFrames);
                this.scan.timeToSample = this.scan.timeFrames;
            end
        end
        
        
        function initResultsArrays(this)
            this.results.phiCh    = zeros(this.stages.xIdxLen, this.stages.yIdxLen, this.scan.zIdxLen, this.scan.numOfQuant, this.scan.channels);
            this.results.phiQuant = zeros(this.stages.xIdxLen, this.stages.yIdxLen, this.scan.zIdxLen, this.scan.numOfQuant);
            this.results.phi      = zeros(this.stages.xIdxLen, this.stages.yIdxLen, this.scan.zIdxLen);
            this.results.lastFFT  = zeros(this.scan.channels, this.scan.samplesPerPos);
        end
            
        function startScan(this,uVars)
            this.setUserVars(uVars)
            this.initResultsArray();

            this.stages.obj.moveStageAbs([this.stages.vars.xVec(1), this.stages.vars.yVec(1)]);
                                      
            for i=1:this.stages.vars.yIdxLen
                
                for j=1:this.stages.vars.xIdxLen
                    this.curPos = ([this.stages.vars.xVec(j), this.stages.vars.yVec(i), 0]);
                    this.stages.obj.moveStageAbs([this.stages.vars.xVec(j), this.stages.vars.yVec(i)]);
                    
                    this.printStr(sprintf("Scaninng on Position: (%.2f, %.2f)\n", this.stages.vars.xVec(j), this.stages.vars.yVec(i)), true);
                    this.startScanTime('singlePos');
                    
                    for k=1:this.scan.numOfQuant(i)
                        this.curPos(3) = k;
    %                         this.printStr(sprintf("Quant number: %d\n", k), true);
                        this.startScanTime('singleQuant');

                        this.startScanTime('netAcoustoOptics');
                        res = this.acoustoOptics.obj.measureAndAnlayse();
                        this.stopScanTime('netAcoustoOptics');

                        this.startScanTime('copyTime');
                        this.results.phiCh(j,i,:,k,:)  = permute(gather(res.phiCh), [2,3,4,5,1]);
                        this.results.phiQuant(j,i,:,k) = gather(res.phi);
                        this.stopScanTime('copyTime');

                        if (this.scan.useQuant)
                            this.shiftSpeckle();
                        end
                        
                        this.stopScanTime('singleQuant');
                        this.printStr(sprintf(this.strings.scan, this.curPos(1), this.curPos(2), this.curPos(3)), true);
                    end
                    this.stopScanTime('singlePos');
                end
            end                            
        end
        
        function calcStagesVars(this)
            this.stages.vars.startPosX = this.uVars.stages.startPosX;
            this.stages.vars.Xstride   = this.uVars.stages.Xstride;
            this.stages.vars.endPosX   = this.uVars.stages.endPosX;
            
            this.stages.vars.startPosY = this.uVars.stages.startPosY;
            this.stages.vars.Ystride   = this.uVars.stages.Ystride;
            this.stages.vars.endPosY   = this.uVars.stages.endPosY;
            
            this.stages.vars.xVec = this.stages.vars.startPosX:this.stages.vars.Xstride:this.stages.vars.endPosX;
            this.stages.vars.yVec = this.stages.vars.startPosY:this.stages.vars.Ystride:this.stages.vars.endPosY;
            this.stages.vars.xLen = abs(this.stages.vars.startPosX - this.stages.vars.endPosX);
            this.stages.vars.yLen = abs(this.stages.vars.startPosY - this.stages.vars.endPosY);
            this.stages.vars.xIdxLen = length(this.stages.vars.xVec);
            this.stages.vars.yIdxLen = length(this.stages.vars.yVec);
            this.stages.vars.xIdx = 1:1:this.stages.vars.xIdxLen;        
            this.stages.vars.yIdx = 1:1:this.stages.vars.yIdxLen;
        end
        
        function setStagesUserVars(this, uVars)
            this.uVars.stages.startPosX = uVars.startPosX;
            this.uVars.stages.Xstride = uVars.Xstride;
            this.uVars.stages.endPosX = uVars.endPosX;
            
            this.uVars.stages.startPosY = uVars.startPosY;
            this.uVars.stages.Ystride = uVars.Ystride;
            this.uVars.stages.endPosY = uVars.endPosY;
        end
        
    end
end

