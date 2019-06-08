close all; 
clear all;
clc;

a = arduino('COM7', 'Uno');
for i = 1:100
    writePWMVoltage(a,'D3',1);
    pause(0.01)
    writePWMVoltage(a,'D3',0);
    pause(0.5)
end

fprintf("Done")