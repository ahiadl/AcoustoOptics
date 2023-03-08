close all;
clear all;
clc;

instrreset;


aofg = AORigolGen();
aofg.init();

clc
uVars = AORigolGen.createUserVars();

uVars.fgClk   = 60e6;

uVars.fSig    = 1.25e6;
uVars.fSqnc   = 10e3;
uVars.delay   = 0;
uVars.sigType = 'Hadamard';
uVars.cycPerPulse = 1;
uVars.amp = 0.4; %1 -> 2.5V
uVars.fsClk   = 20e6;

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