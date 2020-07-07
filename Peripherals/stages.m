classdef stages < handle
    %STAGES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

        devices
        connected
        hardwareAvailable
        
        vars        
        curPos % by device ID
                
        xStage
        yStage
        zStage
        curPosX
        curPosY
        curPosZ
    end
    
    methods
        function this = stages(port)
            this.vars.port = port;
            this.connected = false;
            this.hardwareAvailable = false;
        end
        
        function connect(this)
            try
                [this.vars.ZaberPort, this.devices, this.vars.stepSize] = Zaber_Start(this.vars.port);
                this.vars.numOfStages = length(this.devices);
                this.curPos = zeros(1, this.vars.numOfStages); %this is becuse Zaber_start takes all stages home.
                this.connected = true;
                this.hardwareAvailable = true;
                for i = 1:this.vars.numOfStages
                    this.vars.fullRange(i) = 150; %TODO: create zaber stages object
                end
            catch
                this.connected = false;
                this.hardwareAvailable = false;
                fprintf("STAGES: Can't connect to stages\n")
            end
        end
        
        function moveStageRel(this, Range)
            Zaber_MoveRel(this.devices(1), Range(1));
            Zaber_MoveRel(this.devices(2), Range(2));
            this.curPos(1) = Zaber_ReadPosition(this.devices(1));
            this.curPos(2) = Zaber_ReadPosition(this.devices(2)); 
        end
        
        function moveStageRelId(this, id, rel)
            Zaber_MoveRel(this.devices(id), rel);
            this.curPos(id) = Zaber_ReadPosition(this.devices(id));
        end
        
        function moveStageAbsId(this, id, pos)
            Zaber_MoveAbs(this.devices(id), pos);
            this.curPos(id) = Zaber_ReadPosition(this.devices(id));
        end
        
        function moveStageAbs(this, pos)
            Zaber_MoveAbs(this.devices(1), pos(1));
            Zaber_MoveAbs(this.devices(2), pos(2));
            this.curPos(1) = Zaber_ReadPosition(this.devices(1));
            this.curPos(2) = Zaber_ReadPosition(this.devices(2)); 
        end
        
        function assignStagesAxes(this,ax, stageId)
            axIn = length(ax);
            ax = unique(ax);
            stageId = unique(stageId);
            
            if axIn > length(ax)
                fprintf("Stages: Can't assign 2 stages to the same Axis.\n")
                return
            end
            if axIn > length(stageId)
                fprintf("Stages: Can't assign 1 stage to 2 different Axes.\n")
                    return
            end
            
            this.vars.validAxes(1:3) = false;
            for i = 1:axIn
                switch ax(i)
                    case 'X'
                        this.xStage = this.devices(stageId(i));
                        this.vars.validAxes(1) = true;
                    case 'Y'
                        this.yStage = this.devices(stageId(i));
                        this.vars.validAxes(2) = true;
                    case 'Z'
                        this.zStage = this.devices(stageId(i));
                        this.vars.validAxes(2) = true;
                end
            end
            
            if ~this.vars.validAxes(1)
                this.xStage = [];
            elseif ~this.vars.validAxes(2)
                this.yStage = [];
            elseif ~this.vars.validAxes(3)
                this.zStage = [];
            end
        end
        
        function moveStagesForIdentification(this, stageID)
            range = 10;
            for i=1:this.vars.numOfStages
                if (this.devices(i).DeviceNo == stageID)
                    if this.curPos(i) > this.vars.fullRange(i)-range
                        Zaber_MoveRel(this.devices(i), -range);
                        Zaber_MoveRel(this.devices(i), range);
                    else
                        Zaber_MoveRel(this.devices(i), range);
                        Zaber_MoveRel(this.devices(i), -range);
                    end
                end
            end
        end
        
        function moveStageAxisAbs(this, ax, pos)
            switch ax
                case 'X'
                    Zaber_MoveAbs(this.xStage, pos);
                    this.curPosX = Zaber_ReadPosition(this.xStage);
                case 'Y'
                    Zaber_MoveAbs(this.yStage, pos);
                    this.curPosY = Zaber_ReadPosition(this.yStage);
                case 'Z'
                    Zaber_MoveAbs(this.zStage, pos);
                    this.curPosZ = Zaber_ReadPosition(this.zStage);
            end
        end
        
        function moveHome(this)
            Zaber_HomeAll(this.devices)
        end
        
        function moveStageAbsNonBlocking(this)
            this.devices(1).Protocol.emptybuffer;
            Zaber_Steps = Zaber_Convert_mm2steps(this.devices(1), Range);
            Zaber_Err = this.devices(1).moveabsolute(Zaber_Steps);
            
            if (Zaber_Err) % Error handling
                fprintf('\nSTAGES: %s received error %d while running the Zaber_MoveAbs function!\n\n',Zaber_Device.Name(1:9), Zaber_Err);
            end
            
            this.devices(2).Protocol.emptybuffer;
            Zaber_Steps = Zaber_Convert_mm2steps(this.devices(2), Range);
            Zaber_Err = this.devices(1).moveabsolute(Zaber_Steps);
            
            if (Zaber_Err) % Error handling
                fprintf('\nSTAGES: %s received error %d while running the Zaber_MoveAbs function!\n\n',Zaber_Device.Name(1:9), Zaber_Err);
            end
            
        end
        
        function moveStageRelNonBlocking(this)
            this.devices(1).Protocol.emptybuffer;
            Zaber_Steps = Zaber_Convert_mm2steps(this.devices(1), Range);
            Zaber_Err = this.devices(1).moverelative(Zaber_Steps);
            
            if (Zaber_Err) % Error handling
                fprintf('\nSTAGES: %s received error %d while running the Zaber_MoveAbs function!\n\n',this.devices(1).Name(1:9), Zaber_Err);
            end
            
            this.devices(2).Protocol.emptybuffer;
            Zaber_Steps = Zaber_Convert_mm2steps(this.devices(2), Range);
            Zaber_Err = this.devices(2).moverelative(Zaber_Steps);
            
            if (Zaber_Err) % Error handling
                fprintf('\nSTAGES: %s received error %d while running the Zaber_MoveAbs function!\n\n',this.devices(2).Name(1:9), Zaber_Err);
            end
            
        end
        
        function pos = getPosition(this)
            this.curPos(1) = Zaber_ReadPosition(this.devices(1));
            this.curPos(2) = Zaber_ReadPosition(this.devices(2));
            pos = this.curPos;
        end
        
        function isMoving(this)
            
        end
        
        function ids = getDevicesIDs(this)   
          ids = zeros(1,this.vars.numOfStages);
          for i = 1:this.vars.numOfStages
              ids(i) = this.devices(i).DeviceNo;
          end
        end
    end
end

