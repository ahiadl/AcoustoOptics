classdef scanObj < handle
    %CONSISTENCYOBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %Objects
        acoustoOptics
        graphics
        fileSystem
        owner
        
        %Data
        results
        timeTable
        
        %Control Vars
        scan
        owned
        curPos
        curScan
        strings
        
    end
    
    methods
        
        function this = scanObj(acoustoOpticHandle, owner)
            printStr(this, "----Creating Scan Object----", true);
            
            printStr(this, "1. Acousto Optics Object", true);
            if ~isempty(acoustoOpticHandle)
                this.acoustoOptics.obj = acoustoOpticHandle;
                this.acoustoOptics.intExt = 'ext';
            else
                this.acoustoOptics.obj = acoustoOptics();
                this.acoustoOptics.intExt = 'int';
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
        end
        
        
        function setFigsVars(this, vars)
            this.graphics.figs = vars;
            this.graphics.obj.setGlobalReq(vars);
        end
        
        function setAcoustoOpticsUserVars(this, uVars)
            this.acoustoOptics.vars.uVars = uVars;
            this.acoustoOptics.obj.setMeasVars(this.acoustoOptics.vars.uVars);
            this.acoustoOptics.vars.algoVars = this.acoustoOptics.obj.getAlgoVars();
            this.acoustoOptics.obj.configPeripherals();
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
           vars.scan   = this.scan;
        end
        
        function results = readResults(this)
           results = this.results;
        end
        
        function shiftSpeckle(this)

        end

        % Filesystem Functions
        function setFileSystemUserVars(this, uVars)
%             this.fileSystem = uVars;
%             saveData = uVars.saveFullData || uVars.saveReducedData || uVars.saveFigs || uVars.saveResults || uVars.savePhiChCmplx;
%             if saveData
%                 dateStr = strrep([datestr(datetime('now'))], ':', '-');
%                 this.fileSystem.resDir = [this.fileSystem.resDirPath, '\', dateStr,'-',this.fileSystem.scanName];
%                 this.fileSystem.rawDataDir = [this.fileSystem.resDir, '\rawData'];
%                 this.fileSystem.figsDir = [this.fileSystem.resDir, '\Figures'];
%                 mkdir(this.fileSystem.resDir);
%                 mkdir(this.fileSystem.rawDataDir);
%                 mkdir(this.fileSystem.figsDir);
%             end
%             this.turnLogFileOn()
        end
        
        function saveReducedVars(this, name)
%             curVars = this.getVars();
%             curVars.acoustoOptics.algoVars.timing.tSigVec = [];
%             curVars.acoustoOptics.algoVars.timing.tMeasVec = [];
%             curVars.acoustoOptics.algoVars.timing.tAcqVec = [];
%             curVars.acoustoOptics.algoVars.timing.tQuantVec = [];
%             curVars.acoustoOptics.algoVars.timing.tPosVec = [];
%             save(sprintf("%s%s%s",this.fileSystem.rawDataDir,'\', name), 'curVars', '-v7.3');
        end
        
        function saveData(this, saveReducedVars)
%             saveData = this.fileSystem.saveFullData || this.fileSystem.saveReducedData || this.fileSystem.saveFigs || this.filesystem.savePhiChCmplx;
%             if saveData
%                 startScanTime(this, 'saveData');
%                 rawDataName = sprintf("rawData-F%.2fS%d.mat", this.scan.timeFrames(this.curScan(1)), this.curScan(2)); 
%                 varsName    = sprintf("F%.2fS%d-vars.mat", this.scan.timeFrames(this.curScan(1)), this.curScan(2));
%                 if this.fileSystem.saveFullData
%                     res = this.acoustoOptics.obj.rawData;
%                     save(sprintf("%s%s%s", this.fileSystem.rawDataDir,'\', rawDataName),  'res', '-v7.3');
%                 elseif this.fileSystem.saveReducedData
%                     for i = 1:length(this.fileSystem.fieldsToSave)
%                         if reducedSaveMask(i)
%                             reduced.(this.fileSystem.fieldsToSave(i)) = eval(spintf('res.%s', this.fileSystem.fieldsToSave(i)));
%                         end
%                     end
%                     save([this.fileSystem.rawDataDir, '\' , rawDataName],  'reduced', '-v7.3');
%                     startScanTime(this, 'saveData');
%                     this.saveVars(varsName)
%                     stopScanTime(this, 'saveData');
%                 elseif this.filesystem.savePhiChCmplx
%                     
%                 end
%                 if saveReducedVars
%                     this.saveReducedVars(varsName)
%                 end
%                 stopScanTime(this, 'saveData');  
%             end
        end
        
        function saveResults(this)
%             if this.fileSystem.saveResults
%                 res = this.results;
%                 this.saveVars([this.fileSystem.resDir, '\Vars.mat'])
%                 save([this.fileSystem.resDir, '\Results.mat'], 'res', '-v7.3');
%                 this.turnLogFileOff();  
%             end   
        end
        
        function saveVars(this, name)
            curVars = this.getVars();
            save(name, 'curVars', '-v7.3');
        end
        
        function turnLogFileOn(this)
%             diary([this.fileSystem.resDir, '\log.txt'])
%             diary on
        end
        
        function turnLogFileOff(this)
%             diary off
        end
        
        function totStr = printStr(this, str, dispStr)
            totStr = sprintf('%s-%s\n', datestr(datetime('now'), 'HH:MM:SS.FFF'), str);
            if dispStr
                fprintf(totStr)
            end
        end
        
        % Time Functions
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

        function startScanTime(this, str)
            meas = sprintf(this.strings.timeTable, this.curScan);
            meas = strrep(meas, '.' , '');
            this.timeTable.scan.(meas).(str) = tic;
        end
        
        function stopScanTime(this, str)
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
        
        % Graphics Functions
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

