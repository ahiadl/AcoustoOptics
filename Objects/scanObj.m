classdef scanObj < handle
    %CONSISTENCYOBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        acoustoOptics
        stages    
        laserCtrl
        owner
        owned
        
        results
        timeTable
        
        scan
        curPos
        curScan
        fileSystem
        
        strings
        graphics
    end
    
    methods
        
        function this = scanObj(acoustoOpticHandle, stagesHandle, owner)
            printStr(this, "----Creating Scan Object----", true);
            
            printStr(this, "1. Acousto Optics Object", true);
            if ~isempty(acoustoOpticHandle)
                this.acoustoOptics.obj = acoustoOpticHandle;
                this.acoustoOptics.intExt = 'ext';
            else
                this.acoustoOptics.obj = acoustoOptics();
                this.acoustoOptics.intExt = 'int';
            end
            
            printStr(this, "2. Stages Object", true);
            if ~isempty(stagesHandle)
                this.stages.obj = stagesHandle;
                this.stages.intExt = 'ext';
            else
                this.stages.obj = stages('COM3');
                this.stages.intExt = 'int';
            end
            
            this.owned = false;
            if ~isempty(owner)
                this.owner = owner;
                this.owned = true;
            end
        end
        
        function init(this)
            printStr(this,  "----Initiating Submodules----", true);
            printStr(this, "1. Connecting to Acousto Optics Peripherals", true);
            this.acoustoOptics.obj.init();
            printStr(this, "2. Connecting to Arduino uController", true);
%             this.laserCtrl.Device = arduino('COM7', 'Uno');
%             this.laserCtrl.ch = 'D3';
            printStr(this, "3. Connecting to Stages", true);
            this.stages.obj.connect();
        end
        
        function setUserVars (this, uVars)
            this.setFileSystemUserVars(uVars.fileSystem);
            this.setAcoustoOpticsUserVars(uVars.acoustoOptics);
            this.setScanUserVars(uVars.scan);
            this.setStagesUserVars(uVars.stages);
            this.updateGeneralScanVars();
            this.setGReq(uVars.gReq);
%             this.setGraphicsDynamicVars();
        end
        
        function setGReq(this, gReq)
            this.graphics.gReq = gReq;
            this.graphics.obj.setGlobalReq(gReq);
        end
        
        function setAcoustoOpticsUserVars(this, uVars)
            this.acoustoOptics.vars.uVars = uVars;
            this.acoustoOptics.obj.setMeasVars(this.acoustoOptics.vars.uVars);
            this.acoustoOptics.vars.algoVars = this.acoustoOptics.obj.getAlgoVars();
            this.acoustoOptics.obj.configPeripherals();
        end
        
        function setScanUserVars(this, uVars)
            % This function is a dummy function and should be overlapped
            this.uVars.scan = uVars.scan;
        end
        
        function setStagesUserVars(this, uVars) 
            % This function is a dummy function and should be overlapped
            this.stages.vars.uVars = uVars;
        end
        
        function updateGeneralScanVars(this)
            this.scan.channels = this.acoustoOptics.vars.uVars.channels;
            this.scan.zIdx     = this.acoustoOptics.vars.algoVars.len.zIdx; 
            this.scan.zIdxLen  = this.acoustoOptics.vars.algoVars.len.zIdxLen;
        end 
        
        function vars = getVars(this)
%            vars.timeTable = this.timeTable;
           vars.acoustoOptics = this.acoustoOptics.vars;
           vars.acoustoOptics.uVars.gReq = [];
           vars.stages = this.stages.vars;
           vars.scan = this.scan;
        end
        
        function results = readResults(this)
           results = this.results;
        end
        
        function shiftSpeckle(this)
            this.startScanTime('speckleShift');
            
            this.startTime('wait1');
            pause(0.01)
            this.stopTime('wait1');
            
            this.startTime('arduinoOn');
            writePWMVoltage(this.laserCtrl.Device, this.laserCtrl.ch, 1);
            this.stopTime('arduinoOn');
            
            this.startTime('wait2');
            pause(0.01)
            this.stopTime('wait2');
            
            this.startTime('arduinoOff');
            writePWMVoltage(this.laserCtrl.Device, this.laserCtrl.ch, 0);
            this.stopTime('arduinoOff');
            
            this.stopScanTime('speckleShift');
        end

        function setFileSystemUserVars(this, uVars)
            this.fileSystem = uVars;
            saveData = uVars.saveFullData || uVars.saveReducedData || uVars.saveFigs || uVars.saveResults;
            if saveData
                dateStr = strrep(['Results-',datestr(datetime('now'))], ':', '-');
                this.fileSystem.resDir = [this.fileSystem.resDirPath, '\', dateStr,'-',this.fileSystem.scanName];
                this.fileSystem.rawDataDir = [this.fileSystem.resDir, '\rawData'];
                this.fileSystem.figsDir = [this.fileSystem.resDir, '\Figures'];
                mkdir(this.fileSystem.resDir);
                mkdir(this.fileSystem.rawDataDir);
                mkdir(this.fileSystem.figsDir);
            end
            this.turnLogFileOn()
        end
        
        
        function saveReducedVars(this, name)
            curVars = this.getVars();
            curVars.acoustoOptics.algoVars.timing.tSigVec = [];
            curVars.acoustoOptics.algoVars.timing.tMeasVec = [];
            curVars.acoustoOptics.algoVars.timing.tAcqVec = [];
            curVars.acoustoOptics.algoVars.timing.tQuantVec = [];
            curVars.acoustoOptics.algoVars.timing.tPosVec = [];
            save([this.fileSystem.rawDataDir, '\', name(1:end-1)], 'curVars', '-v7.3');
        end
        
        function saveData(this, saveReducedVars)
            saveData = this.fileSystem.saveFullData || this.fileSystem.saveReducedData || this.fileSystem.saveFigs;
            if saveData
                startScanTime(this, 'saveData');
                rawDataName = sprintf("F%.2fS%d-rawData.mat", this.scan.timeFrames(this.curScan(1)), this.curScan(2)); 
                varsName    = sprintf("F%.2fS%d-vars.mat", this.scan.timeFrames(this.curScan(1)), this.curScan(2));
                if this.fileSystem.saveFullData
                    res = this.acoustoOptics.obj.rawData;
                    save([this.fileSystem.rawDataDir, '\' , rawDataName(1:end-1)],  'res', '-v7.3');
                elseif this.fileSystem.saveReducedData
                    for i = 1:length(this.fileSystem.fieldsToSave)
                        if reducedSaveMask(i)
                            reduced.(this.fileSystem.fieldsToSave(i)) = eval(spintf('res.%s', this.fileSystem.fieldsToSave(i)));
                        end
                    end
                    save([this.fileSystem.rawDataDir, '\' , rawDataName],  'reduced', '-v7.3');
                    startScanTime(this, 'saveData');
                    this.saveVars(varsName)
                    stopScanTime(this, 'saveData');
                end
                if saveReducedVars
                    this.saveReducedVars(varsName)
                end
                stopScanTime(this, 'saveData');  
            end
        end
        
        function saveResults(this)
            if this.fileSystem.saveResults
                res = this.results;
                this.saveVars([this.fileSystem.resDir, '\Vars.mat'])
                save([this.fileSystem.resDir, '\Results.mat'], 'res', '-v7.3');
                this.turnLogFileOff();  
            end   
        end
        
        function saveVars(this, name)
            curVars = this.getVars();
            save(name, 'curVars', '-v7.3');
        end
        
        function turnLogFileOn(this)
            diary([this.fileSystem.resDirPath, '\log.txt'])
            diary on
        end
        
        function turnLogFileOff(this)
            diary off
        end
        
        function totStr = printStr(this, str, dispStr)
            totStr = sprintf('%s-%s\n', datestr(datetime('now'), 'HH:MM:SS.FFF'), str);
            if dispStr
                fprintf(totStr)
            end
        end
        
        function resetTimeTable(this)
           this.timeTable = struct(); 
        end
        
        function startTime(this, str)
            meas = sprintf("%s", str);
            this.timeTable.(meas).(str) = tic;
        end
        
        function stopTime(this, str)
            meas = sprintf("%s", str);
            this.timeTable.(meas).(str) = toc(this.timeTable.(meas).(str));
        end

        function  startScanTime(this, str)
            meas = sprintf(this.strings.timeTable, this.curScan);
            meas = strrep(meas, '.' , '');
            this.timeTable.scan.(meas).(str) = tic;
        end
        
        function  stopScanTime(this, str)
            meas = sprintf(this.strings.timeTable, this.curScan);
            meas = strrep(meas, '.' , '');
            this.timeTable.scan.(meas).(str) = toc(this.timeTable.scan.(meas).(str));
        end
        
        function storeAcoustoOpricTimeTable(this)
            meas = sprintf(this.strings.timeTable, this.curScan);
            meas = strrep(meas, '.' , '');
            tmp = this.acoustoOptics.obj.getTimeTable();
            this.timeTable.scan.(meas).algo = tmp.algo;
            this.timeTable.scan.(meas).digitizer = tmp.digitizer;
            this.timeTable.scan.(meas).acq = tmp.acq;
        end
               
        function timeTable = getTimeTable(this)
           timeTable = this.timeTable;
        end
        
        function gH = getGraphicsHandle(this)
           gH = this.graphics.obj; 
        end
        
        function setOwnerObjectGraphicsUpdate(this)
            % sets the owner object to update the graphics - performance
            % considerations
            this.ownerGraphUpdate = true;
        end
        
        function setExtOwnerGraphicsUpdate(this)
            % sets the object to update the graphics (instead of owner)
            this.ownerGraphUpdate = false;
        end
        
        function restartStages(this)
           this.stages.obj.connect(); 
        end

    end
end

