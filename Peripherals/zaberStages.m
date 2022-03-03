classdef zaberStages < handle
    %STAGES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vars
        
        zPort
        protocol
        devices
        
        connected
    end
    
    methods
        function this = zaberStages(port)
            this.vars.port = port;
            this.connected = false;
        end
        
        function status = connect(this)
            try
                this.zPort = serial(this.vars.port, 'BaudRate',9600);
                set(this.zPort, ...
                    'BaudRate', 9600, ...
                    'DataBits', 8, ...
                    'FlowControl', 'none', ...
                    'Parity', 'none', ...
                    'StopBits', 1);
                fopen(this.zPort);
                this.protocol = Zaber.BinaryProtocol(this.zPort);
                this.devices = this.protocol.finddevices();
                this.vars.numOfDevices = length(this.devices);
                
                this.connected = true;
            catch
                this.connected = false;
            end
            status = this.connected;
        end
        
        function step = getStepSize(this, id)
            deviceName = this.devices(id).Name;
            switch deviceName
                case ('X-VSR40A Stage + unknown peripheral')
                    step = 0.09525;        % [µm]  Lift step size for X-VSR40A
                case ('X-RSW60C Stage + unknown peripheral')
                    step = 0.0009375;      % [deg] Angular step size for X-RSW60C
                case ('T-LSR150B Slide + unknown peripheral')
                    step = 0.49609375;     % [µm]  Linear step size for T-LSR150B
            end
        end
        
        function steps = mm2steps(this, id , range)
            stepSize = this.getStepSize(id);
            deviceName = this.devices(id).Name;
            switch deviceName
                case ('X-VSR40A Stage + unknown peripheral')
                    steps = range*1e3 / stepSize;
                case ('X-RSW60C Stage + unknown peripheral')
                    steps = range / stepSize;
                case ('T-LSR150B Slide + unknown peripheral')
                    steps = range*1e3 / stepSize;
            end
        end
 
        function range = steps2mm(this, id, steps)
            stepSize = this.getStepSize(id);
            deviceName = this.devices(id).Name;
            switch deviceName
                case ('X-VSR40A Stage + unknown peripheral')     
                    range = stepSize * steps * 1e-3;
                case ('X-RSW60C Stage + unknown peripheral')
                    range = stepSize * steps;
                case ('T-LSR150B Slide + unknown peripheral')
                    range = stepSize * steps * 1e-3;
            end
        end

        function moveRel(this, id, rel)
            if id > this.vars.numOfDevices
                fprintf("ZABER: invalid stage id");
            end
            this.devices(id).Protocol.emptybuffer;
            steps = this.mm2steps(id, rel);
            errCode = this.devices(id).moverelative(steps);
            this.devices(id).waitforidle; % Wait for it to get there
            if (errCode) % Error handling
                fprintf('ZABER: stage id: %d received error: %d in moveRel function!\n', id, errCode);
            end
        end
        
        function moveAbs(this, id, pos)
            if id > this.vars.numOfDevices
                fprintf("ZABER: invalid stage id");
                return;
            end
            this.devices(id).Protocol.emptybuffer;
            steps = this.mm2steps(id, pos);
            errCode = this.devices(id).moveabsolute(steps);
            this.devices(id).waitforidle; % Wait for it to get there
            if (errCode) % Error handling
                fprintf('ZABER: stage id: %d received error: %d in moveAbs function!\n', id, errCode);
            end
        end
        
        function posMM = readPosition(this, id)
            if id > this.vars.numOfDevices
                fprintf("ZABER: invalid stage id");
            end
            this.devices(id).Protocol.emptybuffer;
            posSteps = double(this.devices(id).getposition);
            posMM = this.steps2mm(id, posSteps);
        end
        
        function home(this, id)
            if id > this.vars.numOfDevices
                fprintf("ZABER: invalid stage id");
            end
            this.devices(id).Protocol.emptybuffer
            this.devices(id).home;
            this.devices(id).waitforidle; 
        end

        function homeAll(this)
            for i=1:this.vars.numOfDevices
                this.home(i);
            end
        end
        
        function referenceAll(this)
            this.homeAll();
        end
        
        function range = getRange(this)
            %TODO: implement a query object
            range = zeros(1,this.vars.numOfDevices);
            for i = 1:this.vars.numOfDevices
                deviceName = this.devices(i).Name;
                switch deviceName
                    case ('X-VSR40A Stage + unknown peripheral')     
                        range(i) = 40;
                    case ('X-RSW60C Stage + unknown peripheral')
                        range(i) = 60;
                    case ('T-LSR150B Slide + unknown peripheral')
                        range(i) = 150;
                end
            end
        end
        
        function n = getNumOfStages(this)
           n = this.vars.numOfDevices;
        end
        
        function [ax, id] = getDefaultAssign(this)
            ax = [];
            if this.vars.numOfDevices == 1
                ax = 'Z';
                id = 1;
            elseif this.vars.numOfDevices == 2
                ax = ['Z', 'Y'];
                id = [1, 2];
            elseif this.vars.numOfDevices == 3
                ax = ['Z', 'Y', 'X'];
                id = [1, 2, 3];
            end
        end
        
        function [ax, up, low] = getDefaultLimits(this)
            if this.vars.numOfDevices == 1
                ax  = 'X';
                up  = 150;
                low = 0;
            elseif this.vars.numOfDevices == 2
                ax  = ['X', 'Y'];
                up  = [150, 150];
                low = [0, 0];
            elseif this.vars.numOfDevices == 3
                ax  = ['X', 'Y', 'Z'];
                up  = [150 150 150];
                low = [0 0 0];
            end
        end
        
        function setVelocity(this, id, vel)
            % vel in mm/sec
            if id > this.vars.numOfDevices
                fprintf("ZABER: invalid stage id");
            end
            vel = this.devices(id).Units.velocitytonative(vel/1000);
            this.devices(id).Protocol.emptybuffer
            this.devices(id).moveatvelocity(vel);
%             this.devices(id).waitforidle; 
        end
    end
end

