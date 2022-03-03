close all;
clear all;
clc;

digi = Digitizer();
digi.connect();

%% TS
digiVars = Digitizer.uVarsCreate();

fs = 5e6;
timeToSample = 50e-6; %[s]
% bytesPerSample = 2;

digiVars.mode      = 'TS'; 
digiVars.fs        = fs;
digiVars.useGPU    = false;
digiVars.channels  = 1; 

digiVars.triggerDelay = 0;
digiVars.extClk = true;
digiVars.extTrig = true;

digiVars.bufferSizeBytes = 8*(2^20);
digiVars.samplesPerMeas  = timeToSample*fs;

digiVars.exportCropped = false;


digi.setVars(digiVars);
digi.configure();

res = digi.acquire();
figure();
plot(res)

%% NPT
digiVars = Digitizer.uVarsCreate();

digiVars.mode      = 'NPT'; 
digiVars.fs        = 100e6;
digiVars.useGPU    = false;
digiVars.channels  = 1; 

digiVars.triggerDelay = 0;
digiVars.extClk  = false;
digiVars.extTrig = true;

digiVars.timeToSample = 300e-6;
digiVars.avgNum       = 4096;
digiVars.numMeas      = 1;

digiVars.draw = true;

digi.setVars(digiVars);
digi.configure();

res = daq.acquire();
% figure();
% plot(res)


%% Monitor
digiVars = Digitizer.uVarsMonitorCreate();
 
digiVars.fs        = 100e6;
digiVars.channels  = 1; 

digiVars.triggerDelay = 0;
digiVars.extClk  = false;
digiVars.extTrig = true;

digiVars.timeToSample = 300e-6;
digiVars.avgNum       = 2048;


daq.monitor(digiVars);





