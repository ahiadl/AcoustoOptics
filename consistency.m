classdef consistency < scanObj
    %CONSISTENCYOBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        curScan
        
    end
    
    methods
        
        function this = consistency()
            this.strings.scan = "Done Scan for (F,S,Q) = (%8.3f, %d, %d)\n";
            this.strings.timeTable = "F%dS%dQ%d";
            this.curScan = zeros(1,3);
            this.graphics = Graphics();         
            this.initPlotRequests()

        end
        
        
        function plotRequests = initPlotRequests(this)
%             this = plotRequest(type, tit, xlabel, ylabel, xlims, ylims, clims, legend, legVals)
            plotRequests.allFrames     = plotRequest('errBar', "All Frames",                         'Frames[s]', 'Mean \Phi', [], [], [], [], []);
            plotRequests.curFrameSets  = plotRequest('errBar', "Frames = %.2f [s] Sets (Current)",   'Set [#]',   'Mean \Phi', [], [], [], [], []);
            plotRequests.curSetsQuant  = plotRequest('stem',   "Set = %d[#] Quants (Current)",       'Frames[s]', 'Mean \Phi', [], [], [], [], []);
            plotRequests.chosenFrame   = plotRequest('errBar', "Frame: %.2f [s] -Total",             'Frames[s]', 'Mean \Phi', [], [], [], [], []);
            plotRequests.chosenSet     = plotRequest('stem',   "Frame: %.2f [s], set #d - Total",    'Frames[s]', 'Mean \Phi', [], [], [] ,[] ,[]);
            plotRequests.chosenFrameCh = plotRequest('errBar', "Frame: %.2f [s], Channels",          'Frames[s]', 'Mean \Phi', [], [], [], 'Ch = %d', 1:4);
            plotRequests.chosenSetCh   = plotRequest('stem',   "Frame: %.2f [s], ser #d - Channels", 'Frames[s]', 'Mean \Phi', [], [], [], 'Ch = %d', 1:4);
            
%             this.graphics.setPlotRequests(plotRequests);
%             acoustoOptics.createGraphicsRequest();
        end
        
        function setScanUserVars(this, uVars)
            this.scan.timeFrames  = uVars.startTime : uVars.stride : uVars.endTime;
            this.scan.numOfFrames = length(this.scan.timeFrames);
            this.scan.speckleTime = uVars.speckleTime;
            this.scan.useQuant    = uVars.useQuant;
            this.scan.numOfSets   = uVars.numOfSets;
            
            if this.scan.useQuant
                this.scan.numOfQuant   = ceil(this.scan.timeFrames / this.scan.speckleTime);
                this.scan.timeToSample = ones(1, this.scan.numOfFrames) * this.scan.speckleTime;
            else 
                this.scan.numOfQuant   = ones(1, this.scan.numOfFrames);
                this.scan.timeToSample = this.scan.timeFrames;
            end
        end
        
        function setStagesUserVars(this, uVars) 
           this.stages.vars.uVars = uVars;
           this.stages.vars.xPos = uVars.xPos;
           this.stages.vars.yPos = uVars.yPos;
        end
  
        function initResultsArrays(this)
            this.results.phiCh       = zeros(this.scan.zIdxLen, this.scan.numOfFrames, this.scan.numOfSets, this.scan.numOfQuant(end), this.scan.channels);
            this.results.phiQuant    = zeros(this.scan.zIdxLen, this.scan.numOfFrames, this.scan.numOfSets, this.scan.numOfQuant(end));
            this.results.phiSets     = zeros(this.scan.zIdxLen, this.scan.numOfFrames, this.scan.numOfSets);
            this.results.phiFrame    = zeros(this.scan.zIdxLen, this.scan.numOfFrames);
            
            this.results.phiSetsStd  = zeros(this.scan.zIdxLen, this.scan.numOfFrames, this.scan.numOfSets);
            this.results.phiFrameStd = zeros(this.scan.zIdxLen, this.scan.numOfFrames);
        end
        
        function startScan(this, uVars)
            this.setUserVars(uVars); % in the consistency space
            this.resetTimeTable();
            this.initResultsArrays();
            this.stages.obj.moveStageAbs(...
                                         [this.stages.vars.xPos,...
                                          this.stages.vars.yPos]);
            this.curScan = zeros(1,3);
            for i = 1:this.scan.numOfFrames
                this.printStr(sprintf("------------------Starting a New Time Frame----------------------\n"), true);
                
                if ~(this.scan.useQuant)
                    this.acoustoOptics.obj.setSamplingTime(this.scan.timeFrames(i));
                end
                
                this.acoustoOptics.vars.current = this.acoustoOptics.obj.getAlgoVars();
                this.curScan(1)=this.scan.timeFrames(i);
                for j=1:this.scan.numOfSets
                    this.printStr(sprintf("------------------Starting a New Set----------------------\n"), true);
                    this.curScan(2) =  j;
                    
                    for k=1:this.scan.numOfQuant(i)
                        this.curScan(3) = k;
                        this.startScanTime('singleQuant');
                        
                        this.startScanTime('netAcoustoOptics');
                        res = this.acoustoOptics.obj.measureAndAnlayse();
                        this.stopScanTime('netAcoustoOptics');
                        
                        this.startScanTime('copyTime');
                        this.results.phiCh(:,i,j,k,:)  = permute(gather(res.phiCh), [2,3,4,5,1]);
                        this.results.phiQuant(:,i,j,k) = gather(res.phi);
                        this.stopScanTime('copyTime');
                        
                        if (this.scan.useQuant)
                            this.shiftSpeckle();
                        end

                        this.stopScanTime('singleQuant');
                        this.storeAcoustoOpricTimeTable();
                        this.printStr(sprintf(this.strings.scan, this.curScan(1), this.curScan(2), this.curScan(3)), true);
                        
                        this.updateGraphics('Quant');
                    end
                    
                    this.startScanTime('setMean');
                    this.results.phiSets(:,i,j)    = mean(this.results.phiQuant(:,i,j,:), 4);
                    this.results.phiSetsStd(:,i,j) = std(this.results.phiQuant(:,i,j,:), 0, 4);
                    this.stopScanTime('setMean');
                    this.updateGraphics('Set');
                    
                end
                
                this.startScanTime('frameMean');
                this.results.phiFrame(:,i)    = mean(this.results.phiSets(:,i,:), 3);
                this.results.phiFrameStd(:,i) = std(this.results.phiQuant(:,i,:), 0, 3);
                this.stopScanTime('frameMean');
                this.updateGraphics('Frame');
                
                %     updatePlots
                
                if ~(this.scan.useQuant)
                     this.acoustoOptics.obj.resetDigitizer(); %release Buufers
                end
  
            end
        end
        
        function updateGraphics(this, graph)
            switch graph
                case 'Quant'
                    this.graphics.plotRes('curSetsQuant',...
                                          this.curScan(2),...
                                          1:this.scan.numOfQuant(this.curScan(1)),...
                                          this.results.phiQuant(:,this.curScan(1),this.curScan(2),:) );
                case 'Set'
                    this.graphics.plotRes('curFrameSets',...
                                          this.curScan(2),...
                                          1:this.scan.numOfSets,...
                                          this.results.phiSets(:,this.curScan(1),:),...
                                          this.results.phiSetsStd(:,this.curScan(1),:) );
                case 'Frame'
                   this.graphics.plotRes('allFrames',...
                                          this.curScan(2),...
                                          1:this.scan.numOfSets,...
                                          this.results.phiFrame,...
                                          this.results.phiFrameStd ); 
            end
        end
        
        
        
%         function plotResults(this)
%             if this.plotReq.measSamples.getIsValid() && isgraphics(this.plotReq.measSamples.getAx())
%                 [~, ch, ~, ~] = this.plotReq.measSamples.getRequest();
%                     this.displayResults(this.plotReq.measSamples,...
%                                         this.timing.tMeasVec*1e6,...
%                                         this.res.rawData(ch, 1:this.samples.samplesPerMeas)',...
%                                         sprintf('Measured Signal, ch: %d', ch),...
%                                         't[\mus]', 'Voltage');
%             end 
%         end
        
    end
end

