close all;
clear all;
clc;

digi = Digitizer();
digi.connect();

%% TS
digiVars = Digitizer.uVarsCreate();

fs = 5e6;
timeToSample = 0.002; %[s]
% bytesPerSample = 2;

digiVars.mode      = 'TS'; 
digiVars.fs        = fs;
digiVars.useGPU    = false;
digiVars.channels  = 1; 

digiVars.triggerDelay = 0;
digiVars.extClk = true;
digiVars.extTrig = true;
digiVars.draw    = false;

digiVars.bufferSizeBytes = 4*(2^20);
digiVars.samplesPerMeas  = timeToSample*fs;

digiVars.exportCropped = false;
digiVars.croppedSamples = [];

digi.setVars(digiVars);
digi.configure();

res = digi.acquire();
figure();
plot(res)

%% NPT
digiVars = Digitizer.uVarsCreate();

digiVars.mode      = 'NPT'; 
digiVars.fs        = 20e6;
digiVars.useGPU    = false;
digiVars.channels  = 16; 

digiVars.triggerDelay = 0;
digiVars.extClk  = true;
digiVars.extTrig = true;

digiVars.measLenType = 'time';
digiVars.timeToSample = 0.002;
digiVars.avgNum       = 1;
digiVars.numMeas      = 250;

digiVars.draw = false;

digi.setVars(digiVars);
digi.configure();

for i=1:10
    tic
    res = digi.acquire();
    toc
end

figure();
plot(res)



%% Monitor
digiVars = Digitizer.uVarsMonitorCreate();
 
digiVars.fs        = 20e6;
digiVars.channels  = 1; 

digiVars.triggerDelay = 0;
digiVars.extClk  = true;
digiVars.extTrig = true;

digiVars.timeToSample = 0.002;
digiVars.avgNum       = 1;


digi.monitor(digiVars);