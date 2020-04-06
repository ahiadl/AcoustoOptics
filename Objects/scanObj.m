classdef scanObj < handle
    %CONSISTENCYOBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ao
        graphics
        fileSystem
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
        curScan     %[1st, 2nd]     [mm]
        curPosIdx   %[r, xIdx, yIdx][#]
        curScanIdx; %[r, 1st, 2nd]  [#]
        
        owned
    end
    
    methods
           
        function this = scanObj(acoustoOpticHandle, owner)
            printStr(this, "----Creating Scan Object----", true);
            
            printStr(this, "1. Acousto Optics Object", true);
            if ~isempty(acoustoOpticHandle)
                this.ao = acoustoOpticHandle;
%                 this.ao.intExt = 'ext';
            else
                this.ao = acoustoOptics();
%                 this.ao.intExt = 'int';
            end
            
            this.owned = false;
            if ~isempty(owner)
                this.owner = owner;
                this.owned = true;
            end
        end
 
        function vars = getVars(this)
%            vars.timeTable = this.timeTable;
           vars.acoustoOptics = this.acoustoOptics.vars;
           vars.acoustoOptics.uVars.gReq = [];
           vars.stages = this.stages.vars;
           vars.scan   = this.scan;
        end
        
        function results = getResults(this)
           results = this.results;
        end

        % Filesystem Functions

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

    end
end

