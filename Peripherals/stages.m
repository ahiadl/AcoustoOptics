classdef stages < handle
    %STAGES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        curPos
        
        devices
    
        vars
        
        connected
    end
    
    methods
        function this = stages(port)
            this.vars.port = port;
            this.connected = false;
        end
        
        function connect(this)
%             if this.connected
%                 fprintf("Stages are already connected.\n");
%             else
                [this.vars.ZaberPort, this.devices, this.vars.stepSize] = Zaber_Start(this.vars.port);
                this.connected = true;
%             end
        end
        
        function moveStageRel(this, Range)
            Zaber_MoveRel(this.devices(1), Range(1));
            Zaber_MoveRel(this.devices(2), Range(2));
            this.curPos(1) = Zaber_ReadPosition(this.devices(1));
            this.curPos(2) = Zaber_ReadPosition(this.devices(2)); 
        end
        
        function moveStageAbs(this, pos)
            Zaber_MoveAbs(this.devices(1), pos(1));
            Zaber_MoveAbs(this.devices(2), pos(2));
            this.curPos(1) = Zaber_ReadPosition(this.devices(1));
            this.curPos(2) = Zaber_ReadPosition(this.devices(2)); 
        end
        
        function moveHome(this)
            Zaber_HomeAll(this.devices)
        end
        
        function moveStageAbsNonBlocking(this)
            this.devices(1).Protocol.emptybuffer;
            Zaber_Steps = Zaber_Convert_mm2steps(this.devices(1), Range);
            Zaber_Err = this.devices(1).moveabsolute(Zaber_Steps);
            
            if (Zaber_Err) % Error handling
                fprintf('\n%s received error %d while running the Zaber_MoveAbs function!\n\n',Zaber_Device.Name(1:9), Zaber_Err);
            end
            
            this.devices(2).Protocol.emptybuffer;
            Zaber_Steps = Zaber_Convert_mm2steps(this.devices(2), Range);
            Zaber_Err = this.devices(1).moveabsolute(Zaber_Steps);
            
            if (Zaber_Err) % Error handling
                fprintf('\n%s received error %d while running the Zaber_MoveAbs function!\n\n',Zaber_Device.Name(1:9), Zaber_Err);
            end
            
        end
        
        function moveStageRelNonBlocking(this)
            this.devices(1).Protocol.emptybuffer;
            Zaber_Steps = Zaber_Convert_mm2steps(this.devices(1), Range);
            Zaber_Err = this.devices(1).moverelative(Zaber_Steps);
            
            if (Zaber_Err) % Error handling
                fprintf('\n%s received error %d while running the Zaber_MoveAbs function!\n\n',this.devices(1).Name(1:9), Zaber_Err);
            end
            
            this.devices(2).Protocol.emptybuffer;
            Zaber_Steps = Zaber_Convert_mm2steps(this.devices(2), Range);
            Zaber_Err = this.devices(2).moverelative(Zaber_Steps);
            
            if (Zaber_Err) % Error handling
                fprintf('\n%s received error %d while running the Zaber_MoveAbs function!\n\n',this.devices(2).Name(1:9), Zaber_Err);
            end
            
        end
        
        function pos = getPosition(this)
            this.curPos(1) = Zaber_ReadPosition(this.devices(1));
            this.curPos(2) = Zaber_ReadPosition(this.devices(2));
            pos = this.curPos;
        end
        
        function isMoving(this)
            
        end
    end
end

