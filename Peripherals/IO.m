classdef IO < handle
    %IO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vars
        deviceString
        s
        hardwareAvailable
    end
    
    methods (Static)
        function uVars = uVarsCreate()
            uVars.port = [];
            uVars.line = [];
            uVars.mode = []; %outputOnly 
        end
    end
    
    methods
        function this = IO()
           this.vars.numOfLines = 8;
           this.vars.numOfPorts = 3;
           this.vars.allocPorts = zeros(this.vars.numOfPorts, this.vars.numOfLines);
           this.vars.portsMode  = zeros(this.vars.numOfPorts, this.vars.numOfLines);
           this.vars.openPorts  = zeros(this.vars.numOfPorts, this.vars.numOfLines);
        end

        
        function connect(this)
            this.deviceString = daq.getDevices;
            if ~isempty(this.deviceString)
                this.s = daq.createSession('ni');
                this.hardwareAvailable = true;
            else
                this.hardwareAvailable = false;
            end
        end
        
        function allocPorts(this, uVars)
            for i =1:length(uVars.port)
                portAndLine = ['Port', num2str(uVars.port(i)), '/Line', num2str(uVars.line(i))];
                if uVars.mode(i)
                    this.s.addDigitalChannel(this.deviceString.ID, portAndLine, 'InputOnly');
                else
                    this.s.addDigitalChannel(this.deviceString.ID, portAndLine, 'OutputOnly');
                end
                this.vars.allocPorts(uVars.port(i), uVars.line(i)) = 1;
                this.vars.portsMode(uVars.port(i), uVars.line(i))  = uVars.mode(i);
            end
        end
        
        function open(this)
            outputSingleScan(this.s, 1);            
            this.vars.openPorts = this.vars.allocPorts;
        end
        
        function close(this)
            outputSingleScan(this.s,0);
            this.vars.openPorts = zeros(this.vars.numOfPorts, this.vars.numOfLines);
        end
                   
    end
        
end

