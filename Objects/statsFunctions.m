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
            this.strings.timeTable = strings.timeTable;
            this.strings.scan      = strings.scan;
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
        
    end
end

