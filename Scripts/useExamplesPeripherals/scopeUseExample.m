% close all
% clear all
% clc;

% instrreset 

Scope_visa_address = 'USB0::0x1AB1::0x04CE::DS1ZA172115085::0::INSTR';
scope = Rigol_DS1074z('ni',Scope_visa_address);
%%
scopeVars.ch.id = 1;
scopeVars.ch.coup   = 'AC';
scopeVars.ch.invert = 'OFF';
scopeVars.ch.offset = 0;
scopeVars.ch.scale  = 0.005; %TODO: Logic to set the scale;
scopeVars.ch.fineRes = false;

scopeVars.trigger.mode      = 'EDGE';
scopeVars.trigger.trigCh    = 2;
scopeVars.trigger.edge      = 'POSitive';
scopeVars.trigger.trigLevel = 0; %TODO: should be configured according to pulser trigger level/

scopeVars.timebase.mode   = 'MAIN';
scopeVars.timebase.scale  = 2e-5;
scopeVars.timebase.offset = 0;

scopeVars.waveform.chToSample    = 1;
scopeVars.waveform.WavPointsMode = 'NORMal';
scopeVars.waveform.Format        = 'BYTE';

scopeVars.acq.mode     = 'NORM';
scopeVars.acq.numOfAvg = 128;
scopeVars.acq.waitAfterClear = 2;

scope.setGlobalConf(scopeVars);
scope.setChannelVars(scopeVars.ch);

scope.ConfChannelFromVars(1);
scope.ConfAcq();
scope.ConfTrigger();
scope.ConfTimebase();
scope.ConfWaveform();
%%
res = scope.acquire();

figure();
plot(res.tVec, res.data);