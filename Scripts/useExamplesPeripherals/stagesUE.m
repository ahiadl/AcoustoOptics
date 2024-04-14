% close all;
% clear all;
% clc;
% 
% instrreset;
%% PI
% stages = stagesNew('PI');
% 
% stages.connect();

stages.moveAbsAx('X', 90);
stages.moveAbsAx('Y', 50);

stages.moveRelAx('Z', -10);

ids = stages.getLogicIDs();

stages.moveAbsId(1, 20);
stages.moveRelId(2, -30);

pos = stages.getPosition();

stages.home('X');
stages.home(2);

stages.setLimits(['X', 'Y'], [40 50], [20 10]);

stages.moveAbsAx('X', 75);
stages.moveAbsId(1, 75);

stages.moveAbsAx('X', 10);
stages.moveAbsId(1, 10);

%% Zaber
close all;
clear all;
clc;

instrreset;

stagesZ2 = stages('Zaber', 'COM3');

stagesZ2.connect();
%%
stagesZ2.moveAbsAx('X', 90);
stagesZ2.moveAbsAx('Y', 50);
stagesZ2.moveAbsAx('Z', 10);

stagesZ.moveRelAx('X', -10);

ids = stagesZ.getLogicIDs();

stagesZ.moveAbsId(1, 20);
stagesZ.moveRelId(2, -30);

pos = stagesZ.getPosition();

stagesZ.home('X');
stagesZ.home(2);

stagesZ.setLimits(['X', 'Y'], [40 50], [20 10]);

stagesZ.moveAbsAx('X', 75);
stagesZ.moveAbsId(1, 75);

stagesZ.moveAbsAx('X', 10);
stagesZ.moveAbsId(1, 10);

% TODO: add different allocation operation