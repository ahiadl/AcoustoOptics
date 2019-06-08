function [stages] = stagesData(port, startPosX, startPosY, endPosX, endPosY, Xstride, Ystride)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if (exist('port','var')); stages.port = port;
else; stages.port = 'COM3'; end
if (exist('startPosX','var')); stages.startPosX = startPosX;
else; stages.startPosX = 0; end
if (exist('startPosY','var')); stages.startPosY = startPosY;
else; stages.startPosY = 0; end
if (exist('endPosX','var')); stages.endPosX = endPosX;
else; stages.endPosX = 0; end
if (exist('endPosY','var')); stages.endPosY = endPosY;
else; stages.endPosY = 0; end
if (exist('Xstride','var')); stages.Xstride = Xstride;
else; stages.Xstride =0; end
if (exist('Ystride','var')); stages.Ystride = Ystride;
else; stages.Ystride = 0; end
if (exist('Ystride','var')); stages.Ystride = Ystride;
else; stages.Ystride = 0; end
end

