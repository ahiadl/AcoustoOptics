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

uVars.fgClk   =100e6;

uVars.fSig    = 1.25e6;
uVars.fSqnc   = 20e3;
uVars.delay   = 0;
uVars.sigType = 'Hadamard';
uVars.cycPerPulse = 1;
uVars.sigPower = 20; %
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