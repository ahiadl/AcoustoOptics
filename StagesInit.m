function [Stages, status] = StagesInit(Stages)
%STAGESINIT Summary of this function goes here
%   Detailed explanation goes here

[Zaber_Port, Zaber_Devices, Zaber_StepSize] = Zaber_Start(Stages.port);

%bring to the start of the scan

Zaber_MoveAbs(Zaber_Devices(2), Stages.startPosY); 
Zaber_MoveAbs(Zaber_Devices(1), Stages.startPosX); 

Stages.Zaber_Devices = Zaber_Devices;
Stages.Zaber_StepSize = Zaber_StepSize;
Stages.Zaber_Port = Zaber_Port;

end

