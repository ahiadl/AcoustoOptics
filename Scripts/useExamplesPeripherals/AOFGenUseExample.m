close all;
clear all;
clc;

instrreset;

%%
aofg = AOFGen();
aofg.init();

%%
uVars = AOFGen.createUserVars();

uVars.fgClk   =100e6;

uVars.fSig    = 1.25e6;
uVars.fSqnc   = 20e3;
uVars.delay   = 0;
uVars.sigType = 'SinePulse';
uVars.cycPerPulse = 1;

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