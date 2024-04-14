close all;
clear all;
clc;

instrreset;

%%
aofg = AOFGen();
aofg.init();

%%
close all
uVars = AOFGen.createUserVars();

uVars.fgClk   = 100e6;

uVars.fSig    = 1e6;
uVars.fSqnc   = 5e3;
uVars.delay   = 0;
uVars.sigType = 'SinePulse';
uVars.cycPerPulse = 1;
uVars.sigPower = 150; %
uVars.fsClk   =20e6;

aofg.setVars(uVars);
vars = aofg.getVars();
aofg.config();

clkData = aofg.vars.clkData;
sigData = aofg.vars.sigData;

figure();
subplot(1,2,1)
plot(sigData)
subplot(1,2,2)
plot(clkData)

%%
aofg.closeOutputs();

%%
fprintf("Connecting to IO\n");
daqIO = IO();
daqIO.connect();

% Configure the IO
IOVars = daqIO.uVarsCreate;

IOVars.mode = 0; % Output Only
IOVars.port = 1;
IOVars.line = 4;

daqIO.allocPorts(IOVars);
daqIO.open();
daqIO.close();

%%
digi = Digitizer();
digi.connect();

digiVars = Digitizer.uVarsCreate();

timeToSample = 2/uVars.fSqnc; %[s]
digiVars.fs       = 100e6;
digiVars.samplesPerMeas  = timeToSample*digiVars.fs ;

digiVars.mode      = 'TS'; 
digiVars.useGPU    = false;

digiVars.triggerDelay = 0;
digiVars.extClk  = false;
digiVars.extTrig = true;
digiVars.draw    = false;

digiVars.bufferSizeBytes = 64*(2^20);

digiVars.exportCropped = false;

digiVars.channels = 1;

digi.setVars(digiVars);
digi.configure();
curRes = digi.acquire();

figure();
plot(curRes)