function [Stages, status] = StagesInit(Stages)
%STAGESINIT Summary of this function goes here
%   Detailed explanation goes here
status = true;
if ~Stages.dontMoveStages
    [Zaber_Port, Zaber_Devices, Zaber_StepSize] = Zaber_Start(Stages.port);

    %bring to the start of the scan
    if (length(Zaber_Devices) == 2)
        Zaber_MoveAbs(Zaber_Devices(2), Stages.startPosY); 
    end
    Zaber_MoveAbs(Zaber_Devices(1), Stages.startPosX); 
    Stages.Zaber_Devices = Zaber_Devices;
    Stages.Zaber_StepSize = Zaber_StepSize;
    Stages.Zaber_Port = Zaber_Port;
end

Stages.Xpos = Stages.startPosX:Stages.Xstride:Stages.endPosX;
Stages.Ypos = Stages.startPosY:Stages.Ystride:Stages.endPosY;
Stages.Xlen = length(Stages.Xpos);
Stages.Ylen = length(Stages.Ypos);

end

