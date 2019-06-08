classdef scanObj < handle
    %CONSISTENCYOBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        acoustoOptics
        stages    
        laserCtrl
        owner

        results
        timeTable
        
        scan
        curPos      
        fileSystem
        
        strings
        plotReq
        graphics
    end
    
    methods
        
        function this = scanObj(owner)
            printStr(this, "----Creating Scan Object----", true);
            printStr(this, "1. Acousto Optics Object", true);
            this.acoustoOptics.obj = acoustoOptics();
            printStr(this, "2. Stages Object", true);
            this.stages.obj = stages('COM3');
            if nargin == 1
                this.owner = owner;
                this.owned = 1;
            end
        end
        
        function init(this)
            printStr(this,  "----Initiating Submodules----", true);
            printStr(this, "1. Connecting to Acousto Optics Peripherals", true);
            this.acoustoOptics.obj.init();
            printStr(this, "2. Connecting to Arduino uController", true);
            this.laserCtrl.Device = arduino('COM7', 'Uno');
            this.laserCtrl.ch = 'D3';
            printStr(this, "3. Connecting to Stages", true);
            this.stages.obj.connect();
        end
        
        function setUserVars (this, uVars)
            this.setFileSystemUserVars(uVars.fileSystem);
            this.setScanUserVars(uVars.scan);
            this.setAcoustoOpticsUserVars(uVars.acoustoOptics);
            this.setStagesUserVars(uVars.stages);
            this.updateGeneralScanVars();
        end
        
        function setAcoustoOpticsUserVars(this, uVars)
            this.acoustoOptics.vars.uVars = uVars;
            this.acoustoOptics.vars.uVars.algo.timeToSample = this.scan.timeToSample(1);
            this.acoustoOptics.obj.updateUserVars(this.acoustoOptics.vars.uVars);
            this.acoustoOptics.vars.algoVars = this.acoustoOptics.obj.getAlgoVars();
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
            this.scan.channels = this.acoustoOptics.vars.uVars.digitizer.channels;
            this.scan.zIdx     = this.acoustoOptics.vars.algoVars.len.zIdx; 
            this.scan.zIdxLen  = this.acoustoOptics.vars.algoVars.len.zIdxLen;
        end 
        
        function vars = getVars(this)
%            vars.timeTable = this.timeTable;
           vars.acoustoOptics = this.acoustoOptics.vars;
           vars.stages = this.stages.vars;
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
            saveData = uVars.saveFullData || uVars.saveReducedData || saveFigs || saveResults;
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
        
        function saveData(this, res)
            saveData = uVars.saveFullData || uVars.saveReducedData || saveFigs;
            if saveData
                startScanTime(this, 'saveData');
                rawDataName = printStr('F%.2f-S=%d-Q%d-rawData.mat', this.curScanPos(1), this.curScanPos(2), this.curScanPos(3));
                varsName = printStr('F%.2f-S=%d-Q%d-vars.mat', this.curScanPos(1), this.curScanPos(2), this.curScanPos(3));
                if this.fileSystem.saveFullData
                    save([this.fileSystem.rawDataDir, '\' , rawDataName],  'res', '-v7.3');
                    this.saveVars(varsName)
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
                stopScanTime(this, 'saveData');  
            end
        end
        
        function saveResults(this)
            if uVars.saveResults
                this.saveVars()
                save([this.fileSystem.resDir, '\Results.mat'], 'this.results', '-v7.3');
                this.turnLogFileOff();  
            end   
        end
        
        function saveVars(this, name)
            curVars = this.getVars();
            save(name, 'curVars', '-v7.3');
        end
        
        function turnLogFileOn(this)
            diary([this.fileSystem.resDir, '\log.txt'])
            diary on
        end
        
        function turnLogFileOff(this)
            diary off
        end
        
        function totStr = printStr(this, str, dispStr)
            totStr = sprintf('%s: %s\n', datestr(datetime('now'), 'HH:MM:SS.FFF'), str);
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
            meas = sprintf(this.strings.timeTable, this.curScan(1), this.curScan(2), this.curScan(3));
            meas = strrep(meas, '.' , '');
            this.timeTable.scan.(meas).(str) = tic;
        end
        
        function  stopScanTime(this, str)
            meas = sprintf(this.strings.timeTable, this.curScan(1), this.curScan(2), this.curScan(3));
            meas = strrep(meas, '.' , '');
            this.timeTable.scan.(meas).(str) = toc(this.timeTable.scan.(meas).(str));
        end
        
        function storeAcoustoOpricTimeTable(this)
            meas = sprintf(this.strings.timeTable, this.curScan(1), this.curScan(2), this.curScan(3));
            meas = strrep(meas, '.' , '');
            tmp = this.acoustoOptics.obj.getTimeTable();
            this.timeTable.scan.(meas).algo = tmp.algo;
            this.timeTable.scan.(meas).digitizer = tmp.digitizer;
            this.timeTable.scan.(meas).acq = tmp.acq;
            this.timeTable
        end
               
        function timeTable = getTimeTable(this)
           timeTable = this.timeTable;
        end

    end
end

