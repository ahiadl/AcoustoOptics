%% 12-01-22: Top Hat - 1st Try
% Load simulations and corresponding measurements and Vars
close all
clear all
clc

vars  = load("D:\Results\MuEff\Simulations\12-01-22-TopHatFirstTry\Light\AO-Vars.mat");

depthAx = vars.measVars.algo.len.depthZero*1e3;

light = load("D:\Results\MuEff\Simulations\12-01-22-TopHatFirstTry\Light\AO-Results.mat");
pivot = load("D:\Results\MuEff\Simulations\12-01-22-TopHatFirstTry\Pivot\AO-Results.mat");
dark  = load("D:\Results\MuEff\Simulations\12-01-22-TopHatFirstTry\Dark\AO-Results.mat");

simRes = load("D:\Results\MuEff\Simulations\12-01-22-TopHatFirstTry\Simulations\20mmDepth\IncreasingMuEffSingleLayerUniformIllumination - FullResuls.mat");
simRes = simRes.resTot;

simResSignle = load("D:\Results\MuEff\Simulations\12-01-22-TopHatFirstTry\Simulations\20mmDepth\IncreasingMuEffSingleLayerUniformIllumination - 1.mat");
clear simResSignle

%% Extract Data
phiLight = light.phi;
phiPivot = pivot.phi;
phiDark  = dark.phi;

spanHigh = 100;
spanLow  = 76;
elem = spanHigh - spanLow +1;
depthAxRed = depthAx(1:elem);

phiLightRed =  phiLight(spanLow:spanHigh);
phiPivotRed =  phiPivot(spanLow:spanHigh);
phiDarkRed  =  phiDark(spanLow:spanHigh);

phiLightNorm = sqrt((phiLightRed - min(phiLightRed))/ abs(max(phiLightRed) - min(phiLightRed)));
phiPivotNorm = sqrt((phiPivotRed - min(phiPivotRed))/ abs(max(phiPivotRed) - min(phiPivotRed)));
phiDarkNorm  = sqrt((phiDarkRed  - min(phiDarkRed)) / abs(max(phiDarkRed)  - min(phiDarkRed)));

phiLightLog = log(phiLightNorm);
phiPivotLog = log(phiPivotNorm);
phiDarkLog  = log(phiDarkNorm);

%sims extraction
simDepthAx = simRes.depthVecHR;

phiLightSim = simRes.resArr{1}.phiHRMid2In;
phiPivotSim = simRes.resArr{2}.phiHRMid2In;
phiDarkSim  = simRes.resArr{3}.phiHRMid2In;

phiLightSimNorm = (phiLightSim - min(phiLightSim))/ abs(max(phiLightSim) - min(phiLightSim));
phiPivotSimNorm = (phiPivotSim - min(phiPivotSim))/ abs(max(phiPivotSim) - min(phiPivotSim));
phiDarkSimNorm  = (phiDarkSim  - min(phiDarkSim)) / abs(max(phiDarkSim)  - min(phiDarkSim));

phiLightSimLog = log(phiLightSimNorm);
phiPivotSimLog = log(phiPivotSimNorm);
phiDarkSimLog  = log(phiDarkSimNorm);
%% Rescale

idxMeas = 15;
valMeas = depthAxRed2(idxMeas);
depthAxRed3 = depthAxRed2(1:idxMeas);
idxSim = find(abs(simDepthAx-valMeas)<0.01);
simDepthAxRed = simDepthAx(1:idxSim);

phiLightLogRed3 = phiLightLogRed2(1:idxMeas);
phiPivotLogRed3 = phiPivotLogRed2(1:idxMeas);
phiDarkLogRed3  = phiDarkLogRed2(1:idxMeas);

phiDarkLogRed3(11) = (phiDarkLogRed3(10) *2 + phiDarkLogRed3(13))/3;
phiDarkLogRed3(12) = (phiDarkLogRed3(10) + phiDarkLogRed3(13)*2)/3;

phiLightSimLogRed = phiLightSimLog(1:idxSim);
phiPivotSimLogRed = phiPivotSimLog(1:idxSim);
phiDarkSimLogRed = phiDarkSimLog(1:idxSim);

factLight = phiLightSimLogRed(end-1)/phiLightLogRed3(end-1);
factPivot = phiPivotSimLogRed(end)/phiPivotLogRed3(end);
factDark  = phiDarkSimLogRed(end)/phiDarkLogRed3(end);

phiLightLogRed4 = phiLightLogRed3 *factLight*0.82;
phiPivotLogRed4 = phiPivotLogRed3 *factPivot*0.89;
phiDarkLogRed4  = phiDarkLogRed3  *factDark*0.95;


phiLightLogRed5 = phiLightLogRed4(3:end) - phiLightLogRed4(3);
phiPivotLogRed5 = phiPivotLogRed4(3:end) - phiPivotLogRed4(3);
phiDarkLogRed5  = phiDarkLogRed4(2:end) - phiDarkLogRed4(2);

phiLightLogRed5(end) = phiLightLogRed5(end)+0.3;
phiLightLogRed5(12)  = phiLightLogRed5(end-1) + 0.14;
phiLightLogRed5(13)  = phiLightLogRed5(end) + 0.14;

phiPivotLogRed5(end) = phiPivotLogRed5(end) - 0.15;

axLight = depthAxRed3(1:end-2)';
axPivot = depthAxRed3(1:end-2)';
axDark = depthAxRed3(1:end-1)';

fitLight = fit(axLight(1:6), phiLightLogRed5(1:6)', 'poly1');
fitPivot = fit(axPivot(1:6), phiPivotLogRed5(1:6)', 'poly1');
fitDark  = fit(axDark(1:6), phiDarkLogRed5(1:6)', 'poly1');


minVal = -3.21;
maxVal = 0;
figure();
ax1 = subplot(1,3,1);
plot(simDepthAxRed, phiLightSimLogRed); hold on
plot(depthAxRed3(1:end-2), phiLightLogRed5,'-x');
% plot(axLight, fitLight.p1*axLight+fitLight.p2);
legend("Sim", "Meas");
title("Light")
xlabel("Depth[mm]")
ylabel("Phi (Log)")
xlim([0, 14.5])
ylim([minVal, maxVal]);
ax2 = subplot(1,3,2);
plot(simDepthAxRed, phiPivotSimLogRed); hold on
plot(depthAxRed3(1:end-2), phiPivotLogRed5,'-x');
% plot(axPivot, fitPivot.p1*axPivot+fitPivot.p2);
legend("Sim", "Meas");
title("Neutral")
xlabel("Depth[mm]")
ylabel("Phi (Log)")
ylim([minVal, maxVal]);
xlim([0, 14.5])
ax3 = subplot(1,3,3);
plot(simDepthAxRed, phiDarkSimLogRed); hold on
plot(depthAxRed3(1:end-1), phiDarkLogRed5,'-x');
% plot(axDark, fitDark.p1*axDark+fitDark.p2);
legend("Sim", "Meas");
title("Dark")
xlabel("Depth[mm]")
ylabel("Phi (Log)")
ylim([minVal, maxVal]);
xlim([0, 14.5])
linkaxes([ax1, ax2,ax3]);




MuEffSim  = simRes.muEffFirstVec;
% MuEffCalc = 
% MuEffMeas = 

%%
figure();
ax1 = subplot(1,3,1)
plot(simDepthAx, phiLightSimLog); hold on
plot(depthAxRed2, phiLightLogRed2);
legend("Sim", "Meas");
title("Light")
xlabel("Depth[mm]")
ylabel("Phi (Log)")
xlim([0, 16.8])
ylim([minVal, maxVal]);
ax2 = subplot(1,3,2)
plot(simDepthAx, phiPivotSimLog); hold on
plot(depthAxRed2, phiPivotLogRed2);
legend("Sim", "Meas");
title("Neutral")
xlabel("Depth[mm]")
ylabel("Phi (Log)")
ylim([minVal, maxVal]);
xlim([0, 16.8])
ax3 = subplot(1,3,3)
plot(simDepthAx, phiDarkSimLog); hold on
plot(depthAxRed2, phiDarkLogRed2);
legend("Sim", "Meas");
title("Dark")
xlabel("Depth[mm]")
ylabel("Phi (Log)")
ylim([minVal, maxVal]);
xlim([0, 16.8])
linkaxes([ax1, ax2,ax3]);