classdef statsFunctions < handle
    %STATSFUNCTIONS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        objName
        strings
        timeTable
    end
   
    
    methods
        function this = statsFunctions(objName, strings)
            this.objName = objName;
            if nargin > 1
                this.strings.timeTable = strings.timeTable;
                this.strings.scan      = strings.scan;
            end
        end
        
        function setStringModel(this, name, model)
            this.strings.(name) = model;
        end
        
        function printStrModel(this, strName, vars, dispStr)
            str = this.strings.(strName);
            if nargin > 2
                str = sprintf(str, vars);
            end
            totStr = sprintf('%s %s: %s\n', datestr(datetime('now'), 'HH:MM:SS.FFF'),this.objName, str);
            if nargin > 3
                if dispStr
                    fprintf(totStr)
                end
            else
                fprintf(totStr)
            end 
        end
        
        
        function totStr = printStr(this, str, dispStr)
            totStr = sprintf('%s %s: %s\n', datestr(datetime('now'), 'HH:MM:SS.FFF'),this.objName, str);
            if dispStr
                fprintf(totStr)
            end
        end
        
        function startTime(this, str)
            meas = sprintf("%s", str);
            this.timeTable.(meas).(str) = tic;
        end
        
        function stopTime(this, str)
            meas = sprintf("%s", str);
            this.timeTable.(meas).(str) = toc(this.timeTable.(meas).(str));
        end

        function startScanTime(this, str, vars)
            meas = sprintf(this.strings.timeTable, vars);
            meas = strrep(meas, '.' , '');
            this.timeTable.scan.(meas).(str) = tic;
        end
        
        function stopScanTime(this, str, vars)
            meas = sprintf(this.strings.timeTable, vars);
            meas = strrep(meas, '.' , '');
            this.timeTable.scan.(meas).(str) = toc(this.timeTable.scan.(meas).(str));
        end
        
        function resetTimeTable(this)
            this.timeTable = struct();
        end
        
        function timeTable = getTimeTable(this)
            timeTable = this.timeTable;            
        end
        
    end
end

