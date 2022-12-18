darkRes   = load('D:\Results\MuEff\05-Jan-2022 12-24-32-Dark\AO-Results.mat');
darkVars  = load('D:\Results\MuEff\05-Jan-2022 12-24-32-Dark\AO-Vars.mat');
midRes    = load('D:\Results\MuEff\05-Jan-2022 13-19-06-Mid\AO-Results.mat');
midVars   = load('D:\Results\MuEff\05-Jan-2022 13-19-06-Mid\AO-Vars.mat');
lightRes  = load('D:\Results\MuEff\05-Jan-2022 13-55-35-Light\AO-Results.mat');
lightVars = load('D:\Results\MuEff\05-Jan-2022 13-55-35-Light\AO-Vars.mat');



%% Phi vs Z

phiDark  = darkRes.phi;
phiMid   = midRes.phi;
phiLight = lightRes.phi;

phiLightLog = real(log(phiLight));
phiMidLog   = real(log(phiMid));
phiDarkLog  = real(log(phiDark));

depthAx = darkVars.measVars.algo.len.depthVec*1e3;
dDepth = depthAx(2) - depthAx(1);

figure();
subplot(1,2,1)
plot(depthAx, phiLight); hold on;
plot(depthAx, phiMid);
plot(depthAx, phiDark);
xlabel("Recond Depth [mm]")
ylabel("Phi [AU]")
legend("light", "mid", "dark", 'Location', 'northwest')
subplot(1,2,2);
plot(depthAx, phiLightLog); hold on;
plot(depthAx, phiMidLog);
plot(depthAx, phiDarkLog);
legend("light", "mid", "dark", 'Location', 'northwest')
xlabel("Recond Depth [mm]")
ylabel("Phi [log]")

figure();
subplot(1,2,1)
plot(phiLight); hold on;
plot(phiMid);
plot(phiDark);
xlabel("Index")
ylabel("Phi [AU]")
legend("light", "mid", "dark", 'Location', 'northwest')
subplot(1,2,2);
plot(phiLightLog); hold on;
plot(phiMidLog);
plot(phiDarkLog);
xlabel("Index")
ylabel("Phi [log]")
legend("light", "mid", "dark", 'Location', 'northwest')


%% Focus on phantom region
idxHig = 105;
idxLow = 81;

depAxRed = depthAx(1:length(idxLow:idxHig));

figure();
subplot(1,2,1)
plot(depAxRed, flip(phiLight(idxLow:idxHig))); hold on;
plot(depAxRed, flip(phiMid(idxLow:idxHig)));
plot(depAxRed, flip(phiDark(idxLow:idxHig)));
xlabel("Recond Depth [mm]")
ylabel("Phi [AU]")
legend("light", "mid", "dark", 'Location', 'northeast')
subplot(1,2,2);
plot(depAxRed, flip(phiLightLog(idxLow:idxHig))); hold on;
plot(depAxRed, flip(phiMidLog(idxLow:idxHig)));
plot(depAxRed, flip(phiDarkLog(idxLow:idxHig)));
legend("light", "mid", "dark", 'Location', 'northeast')
xlabel("Recond Depth [mm]")
ylabel("Phi [log]")



%% Calculate Gradient
idxHigh = 100;
idxLow  = 89;

% Calculate slope over entire region:
% gradBase   = abs(depthAx(idxHigh) - depthAx(idxLow));
% 
% lightDiff = abs(phiLightLog(idxHigh) - phiLightLog(idxLow));
% midDiff   = abs(phiMidLog(idxHigh)   - phiMidLog(idxLow));
% darkDiff  = abs(phiDarkLog(idxHigh)  - phiDarkLog(idxLow));

% lightGrad = lightDiff / gradBase;
% midGrad   = midDiff   / gradBase;
% darkGrad  = darkDiff  / gradBase;

% Calculate average differential:
dz = depthAx(2) - depthAx(1);

lightGrad = mean(gradient(phiLightLog(idxLow:idxHigh)))./(2*dz)*10;
midGrad   = mean(gradient(phiMidLog(idxLow:idxHigh)))./(2*dz)*10;
darkGrad  = mean(gradient(phiDarkLog(idxLow:idxHigh)))./(2*dz)*10;

factor = [lightGrad/midGrad, darkGrad/midGrad];


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

simRes = load("D:\Results\MuEff\Simulations\FirstMeas-TopHat-40x30-DxW-mm\1LayerSlabMesh-40x30-DxW-mm - FullResults.mat");


% simResSignle = load("D:\Results\MuEff\Simulations\12-01-22-TopHatFirstTry\Simulations\20mmDepth\IncreasingMuEffSingleLayerUniformIllumination - 1.mat");
clear simResSignle

%% Extract Data & Processing
Phantoms = ["Light"; "Pivot"; "Dark"];

phiLight = light.phi;
phiPivot = pivot.phi;
phiDark  = dark.phi;

spanHigh = 92;
spanLow  = 76;
elem = spanHigh - spanLow +1;
depthAxRed = depthAx(1:elem);

% figure()
% plot(phiLight); hold on
% plot(phiPivot);
% plot(phiDark);

%Cut Phantom Forward Region
phiLightRedSq =  phiLight(spanLow:spanHigh);
phiPivotRedSq =  phiPivot(spanLow:spanHigh);
phiDarkRedSq  =  phiDark(spanLow:spanHigh);

% figure()
% plot(phiLightRedSq); hold on
% plot(phiPivotRedSq);
% plot(phiDarkRedSq);

spanLight = max(phiLightRedSq) - min(phiLightRedSq);

% Normalize and Square-root - Normalize Individually:
phiLightNorm = sqrt((phiLightRedSq - min(phiLightRedSq))/ abs(max(phiLightRedSq) - min(phiLightRedSq)));
phiPivotNorm = sqrt((phiPivotRedSq - min(phiPivotRedSq))/ abs(max(phiPivotRedSq) - min(phiPivotRedSq)));
phiDarkNorm  = sqrt((phiDarkRedSq  - min(phiDarkRedSq)) / abs(max(phiDarkRedSq)  - min(phiDarkRedSq)));

% % Normalize and Square-root - Normalize to light:
% phiLightNorm = sqrt((phiLightRedSq - min(phiLightRedSq))/ spanLight);
% phiPivotNorm = sqrt((phiPivotRedSq - min(phiPivotRedSq))/ spanLight);
% phiDarkNorm  = sqrt((phiDarkRedSq  - min(phiDarkRedSq)) / spanLight);

% figure()
% plot(phiLightNorm); hold on
% plot(phiPivotNorm);
% plot(phiDarkNorm);

% Log scale:
phiLightLog = log(phiLightNorm);
phiPivotLog = log(phiPivotNorm);
phiDarkLog  = log(phiDarkNorm);

% figure()
% plot(phiLightLog); hold on
% plot(phiPivotLog);
% plot(phiDarkLog);

% Simulated Profile Extraction:
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

%Extract Fluence Net Profiles: 
idxT = length(phiLightLog);
depthAxRed2     = depthAx(1:idxT);
phiLightLogRed2 = flip(phiLightLog);
phiPivotLogRed2 = flip(phiPivotLog);
phiDarkLogRed2  = flip(phiDarkLog);

% figure()
% plot(phiLightLogRed2); hold on
% plot(phiPivotLogRed2);
% plot(phiDarkLogRed2);

% Calculate Measurements MuEff:
idxHig = 10;
idxLow = 3;
dz = depthAxRed(2) - depthAxRed(1);

depthAxRed(idxLow)
depthAxRed(idxHig)

muEffLight = (abs(mean(gradient(phiLightLogRed2(idxLow:idxHig))))/dz) *10;
muEffPivot = (abs(mean(gradient(phiPivotLogRed2(idxLow:idxHig))))/dz) *10;
muEffDark  = (abs(mean(gradient(phiDarkLogRed2(idxLow:idxHig))))/dz) *10;
muEffMeas  = [muEffLight; muEffPivot; muEffDark];

ratioMeas = muEffMeas ./ muEffMeas(2);

% Simulated result (2in illumination):
idxLow = 10;
idxHig = 50;
dzSim  = simDepthAx(2) - simDepthAx(1);

% figure();
% plot(phiLightSimLog); hold on
% plot(phiPivotSimLog);
% plot(phiDarkSimLog);

muEffLightSim = (abs(mean(gradient(phiLightSimLog(idxLow:idxHig))))/dzSim) *10;
muEffPivotSim = (abs(mean(gradient(phiPivotSimLog(idxLow:idxHig))))/dzSim) *10;
muEffDarkSim  = (abs(mean(gradient(phiDarkSimLog(idxLow:idxHig))))/dzSim) *10;
muEffSim      = [muEffLightSim; muEffPivotSim; muEffDarkSim];

ratioSim  = muEffSim ./ muEffSim(2);

% Ground Truth:
muaGT   = simRes.muaVec'*10;
musPGT  = simRes.musVec'*10;
muEffGT = simRes.muEffVec'*10;
ratioGT = muEffGT./muEffGT(2);

Tprop = table(Phantoms, muaGT, musPGT, muEffGT);

% Calc Errors:
ErrMuEff = abs(muEffMeas-muEffSim)./muEffSim *100;
ErrRatio = abs(ratioMeas-ratioSim)./ratioSim * 100;

% Simulation and Measurement Profiles:
dataVec =[[phiLightSimLog;  phiPivotSimLog;  phiDarkSimLog]', ...
          [phiLightLogRed2, phiPivotLogRed2, phiDarkLogRed2]];


Tres = table(Phantoms, muEffGT, muEffSim, muEffMeas, ErrMuEff, ratioGT, ratioSim, ratioMeas, ErrRatio);

% Relations
relMeasLight = phiLightLogRed2./phiPivotLogRed2;
relMeasDark  = phiPivotLogRed2./phiDarkLogRed2;
relSimLight  = phiLightSimLog./phiPivotSimLog;
relSimDark   = phiPivotSimLog./phiDarkSimLog;

% figure();
% subplot(1,2,1)
% plot(relMeasLight); hold on
% plot(relSimLight);
% subplot(1,2,2)
% plot(relMeasDark); hold on
% plot(relSimDark);

idxHig = 13;
idxLow = 3;
depthAxRed2(idxLow)
depthAxRed2(idxHig)
fitModelMeas1 = fit(depthAxRed2(idxLow:idxHig)', relMeasLight(idxLow:idxHig)', 'poly1');
fitModelMeas2 = fit(depthAxRed2(idxLow:idxHig)', relMeasDark(idxLow:idxHig)', 'poly1');
idxHig = 30;
idxLow = 6;
simDepthAx(idxLow)
simDepthAx(idxHig)
fitModelSim1 = fit(simDepthAx(idxLow:idxHig)',relSimLight(idxLow:idxHig), 'poly1');
fitModelSim2 = fit(simDepthAx(idxLow:idxHig)',relSimDark(idxLow:idxHig), 'poly1');

vec = simDepthAx(idxLow):0.05:simDepthAx(idxHig);
fitMeas1 = fitModelMeas1.p1*vec+fitModelMeas1.p2;
fitMeas2 = fitModelMeas2.p1*vec+fitModelMeas2.p2;
fitSim1 = fitModelSim1.p1*vec+fitModelSim1.p2;
fitSim2 = fitModelSim2.p1*vec+fitModelSim2.p2;

% figure();
% subplot(1,2,1);
% plot(vec, fitMeas1); hold on
% plot(vec, fitSim1);
% subplot(1,2,2)
% plot(vec, fitMeas2); hold on
% plot(vec, fitSim2);

%% Plotting

%% Measurements:
% Phi ^2 vs Z
figure();
subplot(1,2,1)
plot(depthAxRed,phiLightNorm); hold on
plot(depthAxRed,phiPivotNorm);
plot(depthAxRed,phiDarkNorm);
legend("Light", "Pivot", "Dark", 'Location', 'northwest')
ylabel("Phi [AU]")
xlabel("Penetration Depth [mm]")
% xlim([90, 115])
subplot(1,2,2)
plot(depthAxRed,phiLightLog); hold on
plot(depthAxRed,phiPivotLog);
plot(depthAxRed,phiDarkLog);
legend("Light", "Pivot", "Dark", 'Location', 'northwest')
ylabel("Phi (log)")
xlabel("Penetration Depth [mm]")

% Phi^2 vs Idx
figure();
subplot(1,2,1)
plot(phiLightNorm); hold on
plot(phiPivotNorm);
plot(phiDarkNorm);
legend("Light", "Pivot", "Dark", 'Location', 'northwest')
ylabel("Phi^2 (log)")
xlabel("Index")
% xlim([78, 98])
subplot(1,2,2)
plot(phiLightLog); hold on
plot(phiPivotLog);
plot(phiDarkLog);
legend("Light", "Pivot", "Dark", 'Location', 'northwest')
ylabel("Phi^2 (log)")
xlabel("Index")
% xlim([78, 98])

%% Simulated Phi
figure();
subplot(1,2,1)
plot(simDepthAx, phiLightSimNorm); hold on
plot(simDepthAx, phiPivotSimNorm); hold on
plot(simDepthAx, phiDarkSimNorm);  hold on
subplot(1,2,2)
plot(simDepthAx, phiLightSimLog); hold on
plot(simDepthAx, phiPivotSimLog); hold on
plot(simDepthAx, phiDarkSimLog);  hold on

figure();
subplot(1,2,1)
plot(phiLightSimNorm); hold on
plot(phiPivotSimNorm); hold on
plot(phiDarkSimNorm);  hold on
subplot(1,2,2)
plot(phiLightSimLog); hold on
plot(phiPivotSimLog); hold on
plot(phiDarkSimLog);  hold on


%% Simulated Phi and Measured Phi
minVal = mink(dataVec, 7);
minVal = minVal(7);
minVal = -4;
maxVal = max(dataVec);
xlimVal = 14.4;

figure();
subplot(1,3,1)
plot(simDepthAx, phiLightSimLog); hold on
plot(depthAxRed2, phiLightLogRed2);
legend("Sim", "Meas");
title("Light")
xlabel("Depth[mm]")
ylabel("Phi (Log)")
xlim([0, xlimVal])
ylim([minVal, maxVal]);
subplot(1,3,2)
plot(simDepthAx, phiPivotSimLog); hold on
plot(depthAxRed2, phiPivotLogRed2);
legend("Sim", "Meas");
title("Neutral")
xlabel("Depth[mm]")
ylabel("Phi (Log)")
ylim([minVal, maxVal]);
xlim([0, xlimVal])
subplot(1,3,3)
plot(simDepthAx, phiDarkSimLog); hold on
plot(depthAxRed2, phiDarkLogRed2);
legend("Sim", "Meas");
title("Dark")
xlabel("Depth[mm]")
ylabel("Phi (Log)")
ylim([minVal, maxVal]);
xlim([0, xlimVal])

%% Relations and Substraction:
% Simulated Phi relations vs Measured phi relations
% Phi relations = phiLight/PhiPivot

xlimVal = 18;

figure();
ax1 = subplot(1,2,1);
plot(simDepthAx, phiLightSimLog./phiPivotSimLog); hold on
plot(depthAxRed2, phiLightLogRed2./phiPivotLogRed2);
plot(vec, fitMeas1);
legend("Sim", "Meas", "Meas-fit",'Location', 'southwest');
xlabel("Depth[mm]")
ylabel("Phi Light(Log) / Phi Pivot(Log)")
xlim([0, xlimVal])
ax2 = subplot(1,2,2);
plot(simDepthAx, phiPivotSimLog./phiDarkSimLog); hold on
plot(depthAxRed2, phiPivotLogRed2./phiDarkLogRed2);
plot(vec, fitMeas2);
legend("Sim", "Meas", "Meas-fit",'Location', 'southwest');
xlabel("Depth[mm]")
ylabel("Phi Pivot(Log) / Phi Dark(Log)")
xlim([0, xlimVal])
linkaxes([ax1, ax2])

% Simulated Phi substraction vs Measured phi subtraction
% Phi substraction = phiLight - PhiPivot
figure();
subplot(1,2,1)
plot(simDepthAx, phiLightSimLog-phiPivotSimLog); hold on
plot(depthAxRed2, 2*(phiLightLogRed2-phiPivotLogRed2));
legend("Sim", "Meas",'Location', 'northwest');
xlabel("Depth[mm]")
ylabel("2*(Phi Light(Log)- Phi Pivot(Log))")
xlim([0, xlimVal])
subplot(1,2,2)
plot(simDepthAx, phiPivotSimLog - phiDarkSimLog); hold on
plot(depthAxRed2, 2*(phiPivotLogRed2 - phiDarkLogRed2));
legend("Sim", "Meas",'Location', 'northwest');
xlabel("Depth[mm]")
ylabel("2*(Phi Light(Log)- Phi Pivot(Log))")
xlim([0, xlimVal])

%% 04-02-2022: 3D Scan 5cm depth: Simulations & Measurements
close all
clear all
clc

dirName = 'D:\Results\MuEff\5Phantoms-1Layer-3DScan';
for i = 1:5
    str = sprintf("%s/Measurements/Phantom-%d/3DScan-Results.mat", dirName, i);
    res(i) = load(str);
end

str = sprintf("%s/Measurements/Phantom-%d/3DScan-Vars.mat", dirName, i);
vars = load(str);


str = sprintf("%s/Simulations/1LayerSlabMesh-50x30-DxW-mm - FullResults.mat", dirName);
resSim = load(str);


%% Extract Measurements
close all;

numOfPhantoms = length(resSim.muaVec);
pivIdx = floor(numOfPhantoms/2)+1;

for i = 1:numOfPhantoms
    Phantoms(i,1) = sprintf("%d",i);
end

% ----------------------------------------------
% Ground Truth
% ----------------------------------------------

% Extract Ground Truth:
muaGT   = resSim.muaVec'*10;
musPGT  = resSim.musVec'*10;
muEffGT = resSim.muEffVec'*10;
ratioGT = muEffGT./muEffGT(pivIdx);

Tprop = table(Phantoms, muaGT, musPGT, muEffGT)

% ----------------------------------------------
% Measurement
% ----------------------------------------------

% Extract measurements Variables and Results:
depthAx  = vars.grid.depthVec * 1e3;       % in mm
dx       = (depthAx(2) - depthAx(1)) / 10; % in cm
numOfAvg = size(res(1).phi,1) * size(res(1).phi,2);
phiLen   = size(res(1).phi,3);

idxHig = 102;
idxLow = 70;
clear phi phiRed phiNorm avgPhi
for i = 1:numOfPhantoms
    phi(i,:,:)      = reshape(res(i).phi, numOfAvg, phiLen);
    phiRed(i, :, :) = phi(i,:, idxLow:idxHig); 
    phiNorm(i,:,:)  = normMatf(phiRed(i,:,:), 3);
end

% Align All Measurements
[~,I]  = max(phiNorm,[],3);
M = max(I(:));
refIdx = mode(I(:));
minIdx = 3;
maxIdx = M+3;
idxVec = minIdx:1:maxIdx;
noi = length(idxVec);
shiftMat = I - refIdx;

depthAxAlign = depthAx(1:noi);

for i=1:5
    for j=1:9
        tmpIdx       = idxVec + shiftMat(i,j);
        phiAlign(i,j,:) = flip(phiNorm(i,j,tmpIdx));
    end
end

% Average Measurements
avgPhi = squeeze(mean(sqrt(phiAlign),2));

phiNetSq  = sqrt(phiAlign);
phiNetLog = log(phiNetSq);
phiNetLog(phiNetLog == -inf) = NaN;
phiAvgLog = squeeze(mean(phiNetLog,2,'omitnan'));

% Calculate MuEff and Ratio
idxHig = 12;
idxLow = 10;
idxPeak = 6;

dz = depthAx(2) - depthAx(1);
depthAxRed = depthAx(idxPeak:idxHig);

[gx, ~] = gradient(phiAvgLog(:, idxLow:idxHig));
muEffMeas = abs((mean(gx,2) / dz) *10);
ratioMeas = muEffMeas ./ muEffMeas(3);

% ----------------------------------------------
% Simulations
% ----------------------------------------------

% Extract Simulations:
depthVecSim = resSim.depthVec;
dzSim = depthVecSim(2) - depthVecSim(1);

for i=1:numOfPhantoms
    phiSim(i,:) = resSim.resArr{i}.phiHRMid2In;
end

% Normalize & Log:
phiSimNorm = normMatf(phiSim,2);
phiSimLog  = log(phiSimNorm);

% figure();
% plot(phiSimLog')

% Calculate MuEff and Ratio
idxHig = 50;
idxLow = 5;

[gxSim, ~] = gradient(phiSimLog(:, idxLow:idxHig));
muEffSim = abs( (mean(gxSim,2) / dzSim) *10 );
ratioSim = muEffSim / muEffSim(3);

% ----------------------------------------------
% Errors and Analysis
% ----------------------------------------------

% Calc Error
ErrMuEff = abs(muEffMeas-muEffSim)./muEffSim *100;
ErrRatio = abs(ratioMeas-ratioSim)./ratioSim * 100;

Tres = table(Phantoms, muEffGT, muEffSim, muEffMeas, ErrMuEff, ratioGT, ratioSim, ratioMeas, ErrRatio)

% Compare:
idxHig = 12;
idxLow = 6;
depthVecMeasRed = depthAx(idxLow:idxHig) - depthAx(idxLow); 
phiLogMeasRed = phiAvgLog(:, idxLow:idxHig);

idxHig = 16;
idxLow = 1;
depthVecSimRed = depthVecSim(idxLow:idxHig);
phiLogSimRed   =  phiSimLog(:, idxLow:idxHig);

% Phi Relations
idxT = 3;
for i=1:5
    phiRelMeas(i,:) =  phiLogMeasRed(i,idxT:end) ./ phiLogMeasRed(3,idxT:end);
    phiRelSim(i,:)  =  phiLogSimRed(i,idxT:end) ./ phiLogSimRed(3,idxT:end);
    
    fitModelMeas = fit(depthVecMeasRed(idxT:end)', phiRelMeas(i,:)',   'poly1');
    fitModelSim  = fit(depthVecSimRed(idxT:end)',  phiRelSim(i,:)', 'poly1');
    
    fitMeas(i,:) = fitModelMeas.p1*depthVecMeasRed+fitModelMeas.p2;
end

%%
idxT = 5;

h1 = figure();
subplot(2,3,1)
plot(depthAx, squeeze(phi(idxT,:,:))');  
xlabel("Depth [mm]");
ylabel("Phi^2[AU]");
title("Raw")

subplot(2,3,2)
for j=1:9        
    plot(squeeze(phiNorm(i,j,:))); hold on
end
title("Normalization");
xlabel("Depth [mm]");
ylabel("Phi^2[AU]");

subplot(2,3,3)
plot(depthAxAlign, squeeze(phiAlign(idxT,:,:))')
title("Alignment")
xlabel("Depth [mm]");
ylabel("Phi^2[AU]");

subplot(2,3,4)
plot(depthAxAlign, squeeze(avgPhi(idxT,:)))
title("Averaging")
xlabel("Depth [mm]");
ylabel("Phi^2[AU]");

subplot(2,3,5)
plot(depthAxAlign, squeeze(phiAvgLog(idxT,:)))
title("Log")
xlabel("Depth [mm]");
ylabel("Phi[Log]");

for i=1:numOfPhantoms
    legStr{i} = sprintf("Phantom - %d", i);
end

figure();
plot(depthAxAlign, phiAvgLog')
legend(legStr);
xlabel("Depth [mm]");
ylabel("Phi[Log]");

figure()
for i=1:numOfPhantoms
    ax(i) = subplot(2,3,i);
    plot(depthVecSimRed, phiLogSimRed(i,:)); hold on
    plot(depthVecMeasRed, phiLogMeasRed(i,:), '-+');
    legend("Sim", "Meas")
    xlabel("Depth [mm]");
    ylabel("Phi[Log]");
    title(sprintf("Phantom-%d", i));
end
linkaxes(ax);
%%
figure()
for i=1:numOfPhantoms
    ax(i) = subplot(2,3,i);
    plot(depthVecSimRed(idxT:end), phiRelSim(i,:)); hold on
    plot(depthVecMeasRed(idxT:end), phiRelMeas(i,:), '-+');
    plot(depthVecMeasRed, fitMeas(i,:));
    legend("Sim", "Meas", "Fit", 'Location', 'southeast')
    xlabel("Depth [mm]");
    ylabel("Relations to Pivot");
    title(sprintf("Phantom-%d", i));
end

%% Load Transducer Pulse
resPulseObjX = load("../Measurements/Transducer Pressure Field/07-Apr-2022 10-59-12-FocalAxis-FullScan-1.25MHz.mat");
resPulseObjY = load("../Measurements/Transducer Pressure Field/07-Apr-2022 11-35-35-YAxis-FullScan-1.25MHz.mat");
resPulseObjZ = load("../Measurements/Transducer Pressure Field/07-Apr-2022 12-26-26-ZAxis-FullScan-1.25MHz.mat");

scanXVars = resPulseObjX.csVars;
scanXMat  = squeeze(resPulseObjX.resCs);

scanYVars = resPulseObjY.csVars;
scanYMat  = squeeze(resPulseObjY.resCs);

scanZVars = resPulseObjZ.csVars;
scanZMat  = squeeze(resPulseObjZ.resCs);

%%
tVec     = scanXVars.tVec;

%X:
scanVecX1 = scanXVars.scanVec;
p2pXScan = peak2peak(scanXMat,2);
p2pXNorm = normMatf(p2pXScan,1);
p2pXdB = db(p2pXNorm);
[maxValX, idxX] = max(p2pXScan);
sigXF = scanXMat(idxX,:);
sigXAC = sigXF - mean(sigXF);
sigXFNorm = sigXAC/max(abs(sigXAC));

[envUp, envDo] = envelope(sigXFNorm, 10, 'peaks');

% Extract Speed of Sound & Align X-Vec:
dx = scanVecX1(2) - scanVecX1(1);
sigXF2 = scanXMat(idxX+1,:);
for i=1:10
    [~, idxT(i)] = max(abs(scanXMat(idxX +(i-1),:)));
end
dTpeak =  mean(tVec(idxT(2:end)) - tVec(idxT(1:end-1)));
c = dx/dTpeak;
c = 1500e3;
shiftFactor = round(abs(tVec(idxT(1))*c - scanVecX1(idxX)),5);
scanVecX = scanVecX1 + shiftFactor;
xAxisMeas = tVec*c;

%Calculate rayleigh length
pos6dB      = scanVecX1(p2pXdB >= -6);
focalLen    = abs(pos6dB(1) -  pos6dB(end));
focalPoint  = scanVecX(idxX);

%Pulse Parameters:
envNorm    = envUp./max(envUp);
envdB      = db(envNorm);
pos6dB     = xAxisMeas(envdB >= -6);
pulseDepth = abs(pos6dB(1) -  pos6dB(end));

% Y:
scanVecY1 = scanYVars.scanVec;
p2pYScan = peak2peak(scanYMat,2);
p2pYNorm = normMatf(p2pYScan,1);
p2pYdB = db(p2pYNorm);
[maxValY, idxY] = max(p2pYScan);
sigYF = scanYMat(idxY,:);
scanVecY = scanVecY1 - scanVecY1(idxY);

% Z:
scanVecZ1 = scanZVars.scanVec;
p2pZScan = peak2peak(scanZMat,2);
p2pZNorm = normMatf(p2pZScan,1);
p2pZdB = db(p2pZNorm);
[maxValZ, idxZ] = max(p2pZScan);
sigZF = scanZMat(idxZ,:);
scanVecZ = scanVecZ1 - scanVecZ1(idxZ);

% Waist Size:
pos6dBY = scanVecY(p2pYdB >= -6);
pos6dBZ = scanVecZ(p2pZdB >= -6);
waistSize = mean ([abs(pos6dBY(1) -  pos6dBY(end)), abs(pos6dBZ(1) -  pos6dBZ(end))]);

% Max Pressure:
maxPressure = mean([maxValX, maxValY, maxValZ])/837*1e3; %[KPa]

Parameters = ["Focal Point"; "Focal Len"; "Waist Size"; "Pulse Depth"; "Max Pressure"];
Units = ["mm"; "mm"; "mm"; "mm"; "KPa"];
Values = [focalPoint; focalLen; waistSize; pulseDepth; maxPressure];
TransducerParams = table(Parameters, Values, Units)


% Temporal Response along focal axis
x = 60:10:140;
for i=1:length(x)
    [~, idxs(i)] = min(abs(scanVecX - x(i)));
end
sampleSig = scanXMat(idxs,:);

x2 = 54:5:94;
for i=1:length(x)
    [~, idxs2(i)] = min(abs(scanVecX - x2(i)));
end
sampleSig2 = scanXMat(idxs2,:);

%%

% Pressure fiel profiles (Linear)
figure();
subplot(3,1,1)
plot(scanVecX, p2pXNorm);
xlabel("X[mm]")
ylabel("Normalized Pressure [AU]");
subplot(3,1,2)
plot(scanVecY, p2pYNorm);
xlabel("Y[mm]")
ylabel("Normalized Pressure [AU]");
subplot(3,1,3)
plot(scanVecZ, p2pZNorm);
xlabel("Z[mm]")
ylabel("Normalized Pressure [AU]");

% Pressure fiel profiles (dB)
figure();
subplot(3,1,1)
plot(scanVecX, p2pXdB); hold on
plot(scanVecX, -6*ones(1,length(scanVecX)));
xlabel("X[mm]")
ylabel("Normalized Pressure [dB]");
subplot(3,1,2)
plot(scanVecY, p2pYdB); hold on
plot(scanVecY, -6*ones(1,length(scanVecY)));
xlabel("Y[mm]")
ylabel("Normalized Pressure [dB]");
subplot(3,1,3)
plot(scanVecZ, p2pZdB); hold on
plot(scanVecZ, -6*ones(1,length(scanVecZ)));
xlabel("Z[mm]")
ylabel("Normalized Pressure [dB]");

% Temporal shape of pulse at focus:
figure();
subplot(1,2,1)
plot(tVec*1e6, sigXFNorm);
xlim([41, 75])
ylim([-1.1, 1.1])
xlabel('T [\mus]');
ylabel('Pressure [MPa]');
subplot(1,2,2)
plot(xAxisMeas, sigXFNorm); hold on,
plot(xAxisMeas, envUp);
xlim([65.2, 115.3])
ylim([-1.1, 1.1])
xlabel('X [mm]');
ylabel('Pressure [MPa]');

% Pulse shape divergence for different delays:
figure();
for i=1:9
    plot(xAxisMeas, sampleSig(i,:)); hold on
end
xlabel('X [mm]')
ylabel('Pressure [MPa]');
xlim([55, 150])

figure();
for i=1:9
    plot(xAxisMeas, sampleSig2(i,:)); hold on
end
xlabel('X [mm]')
ylabel('Pressure [MPa]');
xlim([55, 100])

% Sinogram on X:
figure();
imagesc(xAxisMeas, scanVecX, scanXMat);
axis tight equal

%% Convolve Simulations:
dxt        = 0.005;
alignPoint = focalPoint-10;

xtVec   = xAxisMeas(1):dxt:xAxisMeas(end);
sigFInt = interp1(xAxisMeas, sigXFNorm, xtVec);

[envupInt, ~] = envelope(sigFInt, 10, 'peaks');

[~, pulsePivotIdx]    = max(sigFInt);
pulsePivotPos         = xtVec(pulsePivotIdx);
pulseNegDist          = 2;
pulsePosDist          = 15;
pulseNegPos           = pulsePivotPos - pulseNegDist;
pulsePosPos           = pulsePivotPos + pulsePosDist;
[~, pulseNegWidthIdx] = min(abs(xtVec - pulseNegPos));
[~, pulsePosWidthIdx] = min(abs(xtVec - pulsePosPos));
pulseFullWidthIdx     = pulsePosWidthIdx - pulseNegWidthIdx + 1;
pulseFullLen          = pulseFullWidthIdx*dxt;

figure();
subplot(1,2,1)
plot(xtVec, sigFInt);
xlim([pulseNegPos, pulsePosPos])
ylim([-1.1, 1.1])
xlabel('X [mm]');
ylabel('Pressure [MPa]');
subplot(1,2,2)
plot(sigFInt); hold on,
plot(envupInt);
xlim([pulseNegWidthIdx, pulsePosWidthIdx])
ylim([-1.1, 1.1])
xlabel('X [mm]');
ylabel('Pressure [MPa]');

%-----------------------------------------------------------
% Interpulate Simulations to match temporal Resolution:
%------------------------------------------------------------
% Interpolate simulation for 0.01 mm res
simDepth =  abs(depthVecSim(end)- depthVecSim(1));
depthVecSimInt = depthVecSim(1) : dxt: depthVecSim(end);
phIdxVec = 1:1:numOfPhantoms;

[X, Y]    = meshgrid(depthVecSim, phIdxVec);
[Xq, Yq]  = meshgrid(depthVecSimInt, phIdxVec);
phiSimInt = interp2(X,  Y, phiSimNorm, Xq, Yq);

% Align to transducer direction
depthVecSimAlign = flip(abs(depthVecSimInt - alignPoint));
phiSimAlign      = flip(phiSimInt,2);

%------------------------------------------------------------
% Pad with zeros for convolution
%------------------------------------------------------------
%simulation
phiSimIntPad   = [ zeros(numOfPhantoms,pulseFullWidthIdx), phiSimAlign,  zeros(numOfPhantoms,pulseFullWidthIdx)]';
depthVecSimPad = [ dxt*(-pulseFullWidthIdx:1:-1)+depthVecSimAlign(1), depthVecSimAlign,  dxt*(1:pulseFullWidthIdx)+depthVecSimAlign(end)];

% if 

%-------------------------------------------------------------
% Interpolate field to match spatial resolution:
%-------------------------------------------------------------
% Interpolate1: interpolate to a dxt res on time axis:
[X,Y]       = meshgrid(xAxisMeas, scanVecXPad);
[Xq, Yq]    = meshgrid(xtVec,scanVecXPad);
scanXMatInt = interp2(X,Y, scanXMatPad, Xq, Yq);

% Cut Relevant Indices On Time Axis:
[xtVec(1), xtVec(end)]
[depthVecSimPad(1), depthVecSimPad(end)]
mask      = (xtVec >= depthVecSimPad(1)) & (xtVec <= (depthVecSimPad(end)+dxt));
xtRed     = xtVec(mask);
scanXMatRed1 = scanXMatInt(:, mask);

% Interpolate2: interpolate to a dxt res on scan axis:
[X,Y]        = meshgrid(xtRed, scanVecXPad);
scanVecInt   = scanVecXPad(1):dxt:scanVecXPad(end);
[Xq, Yq]     = meshgrid(xtRed, scanVecInt);
scanXMatInt2 = interp2(X,Y, scanXMatRed1, Xq, Yq, 'cubic');

% Cut Relevant Indices On Scan Axis:
[scanVecInt(1), scanVecInt(end)]
[depthVecSimPad(1), depthVecSimPad(end)]

mask        = (scanVecInt >= depthVecSimPad (1)) & (scanVecInt <= (depthVecSimPad (end)+dxt));
scanVecRed  = scanVecInt(mask);
scanMat     = scanXMatInt2(mask,:)';

% figure(); imagesc(scanVecRed, xtRed, scanMat); axis tight equal

scanMatNorm = scanMat/max(scanMat(:));
scanMatNorm = scanMatNorm - mean(scanMatNorm(:));
%% convolve: 
[envMat, ~]  = envelope(scanMatNorm', 100, 'peaks');
envMat = envMat';
envMatNorm = envMat/max(envMat(3000*size(envMat,1):((size(envMat,1)^2-size(envMat,1)*1000))));

% figure(); 
% plot(envMatNorm(1:1000:end,:)'); hold on
% plot(phiSimIntPad(:,1))


phiSimConv = envMatNorm*phiSimIntPad;
phiSimConvNorm = normMatf(phiSimConv,1);
phiSimConvLog = log(phiSimConvNorm);

% figure();
% subplot(1,2,1)
% plot(phiSimConv); hold on
% subplot(1,2,2)
% plot(phiSimConvLog)


figure()
for i=1:numOfPhantoms
    subplot(3,2,i)
    plot(depthVecSimPad, log(phiSimIntPad(:,i))); hold on
    plot(scanVecRed, phiSimConvLog(:,i))
    xlim([46,86])
    legend("Clean Sim", "Conv Sim", 'Location', 'northwest')
end



    


