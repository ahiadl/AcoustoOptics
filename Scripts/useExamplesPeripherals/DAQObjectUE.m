% DAQObjectUseExample
close all
clear all
clc;

addpath(genpath(pwd))

daq = DAQ();

uVars.inputImpedance = 50; %not supported
uVars.activeChannels  = 1:256; %Mask

% Timing Vars
uVars.delay = 0; %not supported
uVars.numOfSamples   = 2030;    %not supported

% Acquisition Vars
uVars.numOfAvg    = 1;
uVars.numOfFrames = 1;
            
daq.setUserVars(uVars);
daq.configAndConnect();

hf = figure;
axRD = axes();
hRD = imagesc(axRD, 'XData', 1:256, 'YData', 1:2030, 'CData', zeros(2030,256));
hTitRD = title(axRD, 'Init');
axis(axRD, 'tight')

for i=1:100
    sigMat = daq.acquire();
    set(hRD, 'CData', log(sigMat));
end

daq.disconnect();