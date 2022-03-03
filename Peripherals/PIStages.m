classdef PIStages < handle
    %STAGES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vars
        
        zPort
        protocol
        devCntl
        
        controller
        availableAxes % assumes one axes per stage
        connected   
        
        stageType
        serialNumbers
    end
    
    methods
        function this = PIStages()%DONE %CHECKED
            addpath('C:\Users\Public\PI\PI_MATLAB_Driver_GCS2');
            if(~isa(this.controller,'PI_GCS_Controller'))
                this.controller = PI_GCS_Controller();
            end
            this.connected = false;
            this.vars.res = 0.001;
        end
        
        function status = connect(this)%CHECKED %DONE
            try
                this.serialNumbers = {'0165500349' ; ... %id 1
                                 '0165500326' ; ... %id 2
                                 '0020550301'};     %id 3
                
                this.stageType = {'M-403.6PD';...
                             'M-403.6PD';...
                             'M-414.3PD'};
                
                 %TODO: query the followinf parameters from the
                 %controller
                this.vars.homePosType = [75, 75, 150];         
                this.vars.maxVec = [10, 10, 100];
                
                this.vars.validId = [];
                this.vars.maxDevices = length(this.serialNumbers);
                this.controller.EnumerateUSBAsArray()
                
                this.vars.homePos = [];
                this.availableAxes = {};
                for i = 1: this.vars.maxDevices
                    try
                        this.devCntl{i} = this.controller.ConnectUSB ( this.serialNumbers{i} );
                        this.vars.connectedMask(i) =  this.devCntl{i}.IsConnected;
                    catch
                        this.vars.connectedMask(i) = false;
                    end
                    
                    if this.vars.connectedMask(i)
                        this.vars.validId = [this.vars.validId, i];
                        this.vars.homePos = [this.vars.homePos, this.vars.homePosType(i)];
                        fprintf("PI: %d: Connected to %s: %s.\n", i, this.stageType{i}, this.serialNumbers{i});
                    else
                        fprintf("PI: %d: Not Connected to %s: %s.\n", i, this.stageType{i}, this.serialNumbers{i});
                    end
                    
                end 

                this.vars.numOfDevices = sum(this.vars.connectedMask);
                 
                for i =1:this.vars.numOfDevices
                    id = this.vars.validId(i);
%                     this.devCntl{id}.InitializeController();
                    this.availableAxes(id) = this.devCntl{id}.qSAIasArray();
                    this.devCntl{id}.SVO(this.availableAxes{id}, 1);
                    this.devCntl{id}.RON(this.availableAxes{id}, 0)
%                     this.devCntl{id}.CST(this.availableAxes{id}, this.stageType{id});
%                     this.reference(id);
                end
                
                if this.vars.numOfDevices > 0
                    this.connected = true;
                else
                    this.connected = false;
                end
                
            catch
                this.connected = false;
            end
            status = this.connected;
        end

        function disconnect(this)
            if this.connected
                for i = 1: this.vars.maxDevices
                        try
                            this.devCntl{i}.CloseConnection;
                            this.vars.connectedMask(i) =  this.devCntl{i}.IsConnected;
                            fprintf("PI: %d: Disonnected from %s: %s.\n", i, this.stageType{i}, this.serialNumbers{i});
                        catch
                            this.vars.connectedMask(i) = true;
                        end
                end
            end
        end
        
        function delete(this)
            this.disconnect()
            this.controller.Destroy;
        end

        function moveRel(this, id, rel)
            if ~this.checkConnectedAndValid(id)
                return;
            end
            
            pos = this.readPosition(id) + rel;
            this.moveAbs(id, pos)
        end %DONE %CHECKED
        
        function moveAbs(this, id, pos) %DONE %CHECKED
            if ~this.checkConnectedAndValid(id)
                return;
            end
            
            this.devCntl{id}.MOV(this.availableAxes{id}, pos);

            while(abs(this.readPosition(id) - ...
                      this.devCntl{id}.qMOV(this.availableAxes{id})) > ...
                      this.vars.res)
                pause(0.001);
            end
        end
        
        function blockWhileMoving(this, id)
            if ~this.checkConnectedAndValid(id)
                return;
            end
            
            while(abs(this.readPosition(id) - ...
                      this.devCntl{id}.qMOV(this.availableAxes{id})) > ...
                      this.vars.res)
                pause(0.001);
            end
        end
        
        function moveRelNB(this, id, rel)
            if ~this.checkConnectedAndValid(id)
                return;
            end
            
            pos = this.readPosition(id) + rel;
            this.moveAbsNB(id, pos)
        end %DONE %CHECKED
        
        function moveAbsNB(this, id, pos) %DONE %CHECKED
            if ~this.checkConnectedAndValid(id)
                return;
            end
            
            this.devCntl{id}.MOV(this.availableAxes{id}, pos);
        end

        function setVelocity(this, id, vel)
            this.devCntl{id}.VEL(this.availableAxes{id}, vel);
        end
        
        function vel = getVelocity(this, id)
            vel =  this.devCntl{id}.qVEL(this.availableAxes{id});
        end
        
        function status = checkConnectedAndValid(this, id)
            if id > this.vars.maxDevices
                fprintf("PI: invalid stage id\n");
                status = false;
                return;
            end
 
            if ~this.vars.connectedMask(id)
                fprintf("PI: notConnected stage id\n");
                status = false;
                return
            end
            status = true;
        end %DONE %CHECKED
        
        function pos = readPosition(this, id) % DONE
            if ~this.checkConnectedAndValid(id)
                return;
            end
      
            pos = this.devCntl{id}.qPOS(this.availableAxes{id});

        end %CHECKED
        
        function referenceAll(this)
            for i = 1: this.vars.maxDevices
                id = this.vars.validId(i);
                this.reference(id);
            end
            
            for i=1: this.vars.maxDevices
                while(abs(this.readPosition(id) - ...
                          this.devCntl{id}.qMOV(this.availableAxes{id})) > ...
                          this.vars.res)
                    pause(0.001);
                end
            end
        end
        
        function reference(this, id)
            if ~this.checkConnectedAndValid(id)
                return;
            end
            
            this.devCntl{id}.SVO(this.availableAxes{id}, 1);
            this.devCntl{id}.FRF(this.availableAxes{id});
            
            while this.devCntl{id}.qFRF(this.availableAxes{id})
                pause(0.1);
            end
            
        end %DONE %CHECKED
        
        function home(this, id)
            if ~this.checkConnectedAndValid(id)
                return;
            end
            this.moveAbs(id, this.vars.homePos(id))
%             this.reference(id);
        end %DONE %CHECKED

        function homeAll(this)
            for i=1:this.vars.numOfDevices
                id = this.vars.validId(i);
                this.home(id);
            end
        end %DONE %CHECKED
        
        function range = getRange(this)
            %TODO: implement a query object
            range = zeros(1,this.vars.numOfDevices);
            for i = 1:this.vars.numOfDevices
                id = this.vars.validId(i);
                range(i) = this.devCntl{id}.qTMX(this.availableAxes{1});
            end
        end %DONE %CHECKED
        
        function n = getNumOfStages(this)
           n = this.vars.numOfDevices;
        end %DONE %CHECKED
        
        function [ax, id] = getDefaultAssign(this)
            ax = [];
            id = [];
            if this.devCntl{1}.IsConnected
                ax = 'X';
                id = [id, 1];
            end
            
            if this.devCntl{2}.IsConnected
                ax = [ax, 'Y'];
                id = [id, 2];
            end
            
            if this.devCntl{3}.IsConnected
                ax = [ax, 'Z'];
                id = [id, 3];
            end
        end %DONE %CHECKED
        
        function [ax, up, low] = getDefaultLimits(this) %DONE %CHECKED
            ax  = [];
            up  = [];
            low = [];
            
            if this.devCntl{1}.IsConnected
                ax = [ax, 'X'];
                up = [up, 150];
                low = [low, 0];
            end
            
            if this.devCntl{2}.IsConnected
                ax  = [ax, 'Y'];
                up  = [up, 150];
                low = [low, 0];
            end
            
            if this.devCntl{3}.IsConnected
                ax = [ax, 'Z'];
                up = [up, 255]; %True limit is 300 but tank limits to 255
                low = [low, 0];
            end
        end
        
        function ids = getIdAssign(this)
            ids = this.vars.validId;
        end
        
        function turnOnIO(this, id, ch)
            this.setTriggerPolarity(id, 0);
            this.devCntl{id}.DIO(ch, 0)
        end
        
        function turnOffIO(this, id, ch)
            this.devCntl{id}.DIO(ch, 1)
        end
        
        function setTrigger(this, id, ch, mode, params) 
            %ch and mode will be used once PI will fix the single Position
            %trigger mode.
            
            % params should include: fields with dest, trigPos
            switch mode
                case 'singlePosition'
                    this.devCntl{id}.VAR('dest', num2str(params.dest));
                    this.devCntl{id}.VAR('trigPos',num2str(params.trigPos));
                    if params.trigPos > params.dest %moving backwards
                        this.vars.routine{id} = 'sBack';
                    else % moving forward
                       this.vars.routine{id} ='sFor';
                    end
            end
        end
        
        function setTriggerPolarity(this, id, pol)
            this.devCntl{id}.CTO(1,7,pol);
        end
        
        function startRoutine(this, id)
             this.devCntl{id}.MAC_START(this.vars.routine{id})
        end
        
        function maxVel = getMaxVelocity(this, id)
            maxVel = this.vars.maxVec(id);
        end
    end
end

