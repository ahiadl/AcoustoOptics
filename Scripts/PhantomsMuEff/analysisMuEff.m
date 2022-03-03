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
legend("light", "mid", "dark")
subplot(1,2,2);
plot(depthAx, phiLightLog); hold on;
plot(depthAx, phiMidLog);
plot(depthAx, phiDarkLog);
legend("light", "mid", "dark")

figure();
subplot(1,2,1)
plot(phiLight); hold on;
plot(phiMid);
plot(phiDark);
legend("light", "mid", "dark")
subplot(1,2,2);
plot(phiLightLog); hold on;
plot(phiMidLog);
plot(phiDarkLog);
legend("light", "mid", "dark")

%% Calculate Gradient
idxHigh = 131;
idxLow  = 139;
gradBase   = abs(depthAx(idxHigh) - depthAx(idxLow));

lightDiff = abs(phiLightLog(idxHigh) - phiLightLog(idxLow));
midDiff   = abs(phiMidLog(idxHigh)   - phiMidLog(idxLow));
darkDiff  = abs(phiDarkLog(idxHigh)  - phiDarkLog(idxLow));

lightGrad = lightDiff / gradBase;
midGrad   = midDiff   / gradBase;
darkGrad  = darkDiff  / gradBase;

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

%% Plotting

% Phi ^2 vs Z
figure();
subplot(1,2,1)
plot(depthAxRed,phiLightNorm); hold on
plot(depthAxRed,phiPivotNorm);
plot(depthAxRed,phiDarkNorm);
legend("Light", "Pivot", "Dark", 'Location', 'northwest')
title("Phi")
% xlim([90, 115])
subplot(1,2,2)
plot(depthAxRed,phiLightLog); hold on
plot(depthAxRed,phiPivotLog);
plot(depthAxRed,phiDarkLog);
legend("Light", "Pivot", "Dark", 'Location', 'northwest')
title("Phi (log)")

% Phi^2 vs Idx
figure();
subplot(1,2,1)
plot(phiLightNorm); hold on
plot(phiPivotNorm);
plot(phiDarkNorm);
legend("Light", "Pivot", "Dark", 'Location', 'northwest')
title("Phi^2")
% xlim([78, 98])
subplot(1,2,2)
plot(phiLightLog); hold on
plot(phiPivotLog);
plot(phiDarkLog);
legend("Light", "Pivot", "Dark", 'Location', 'northwest')
title("Phi^2 (log)")
% xlim([78, 98])

% Simulated Phi
figure();
subplot(1,2,1)
plot(simDepthAx, phiLightSim); hold on
plot(simDepthAx, phiPivotSim); hold on
plot(simDepthAx, phiDarkSim);  hold on
subplot(1,2,2)
plot(simDepthAx, phiLightSimLog); hold on
plot(simDepthAx, phiPivotSimLog); hold on
plot(simDepthAx, phiDarkSimLog);  hold on

% Simulated Phi and Measured Phi
idx = 17;
depthAxRed2 = depthAx(1:idx);
phiLightLogRed2 = flip(phiLightLog(1:idx));
phiPivotLogRed2 = flip(phiPivotLog(1:idx));
phiDarkLogRed2  = flip(phiDarkLog(1:idx));

figure();
subplot(1,3,1)
plot(simDepthAx, phiLightSimLog); hold on
plot(depthAxRed2, phiLightLogRed2);
legend("Sim", "Meas");
title("Light")
xlabel("Depth[mm]")
ylabel("Phi (Log)")
subplot(1,3,2)
plot(simDepthAx, phiPivotSimLog); hold on
plot(depthAxRed2, phiPivotLogRed2);
legend("Sim", "Meas");
title("Pivot")
xlabel("Depth[mm]")
ylabel("Phi (Log)")
subplot(1,3,3)
plot(simDepthAx, phiDarkSimLog); hold on
plot(depthAxRed2, phiDarkLogRed2);
legend("Sim", "Meas");
title("Dark")
xlabel("Depth[mm]")
ylabel("Phi (Log)")


% Simulated Phi relations vs Measured phi relations
% Phi relations = phiLight/PhiPivot
figure();
subplot(1,2,1)
plot(simDepthAx, phiLightSimLog./phiPivotSimLog); hold on
plot(depthAxRed2, phiLightLogRed2./phiPivotLogRed2);
legend("Sim", "Meas");
title("Light")
xlabel("Depth[mm]")
ylabel("Phi Light(Log) / Phi Pivot(Log)")
subplot(1,2,2)
plot(simDepthAx, phiPivotSimLog./phiDarkSimLog); hold on
plot(depthAxRed2, phiPivotLogRed2./phiDarkLogRed2);
legend("Sim", "Meas");
title("Dark")
xlabel("Depth[mm]")
ylabel("Phi Pivot(Log) / Phi Dark(Log)")

% Simulated Phi substraction vs Measured phi subtraction
% Phi substraction = phiLight - PhiPivot
figure();
subplot(1,2,1)
plot(simDepthAx, phiLightSimLog-phiPivotSimLog); hold on
plot(depthAxRed2, phiLightLogRed2-phiPivotLogRed2);
legend("Sim", "Meas");
title("Light")
xlabel("Depth[mm]")
ylabel("Phi Light(Log)- Phi Pivot(Log)")
subplot(1,2,2)
plot(simDepthAx, phiPivotSimLog - phiDarkSimLog); hold on
plot(depthAxRed2, phiPivotLogRed2 - phiDarkLogRed2);
legend("Sim", "Meas");
title("Dark")
xlabel("Depth[mm]")
ylabel("Phi Light(Log)- Phi Pivot(Log)")


highIdx = 15;
lowIdx  = 8;

dDepth = depthAx(highIdx) - depthAx(lowIdx);

gLight = (phiLightLog(highIdx) - phiLightLog(lowIdx))/ 2 / dDepth * 10; % /2 = phi^2 => phi ; * 10 = (1/mm => 1/cm) 
gPivot = (phiPivotLog(highIdx) - phiPivotLog(lowIdx))/ 2 / dDepth * 10; % /2 = phi^2 => phi ; * 10 = (1/mm => 1/cm) 
gDark  = (phiDarkLog(highIdx) - phiDarkLog(lowIdx))  / 2 / dDepth * 10; % /2 = phi^2 => phi ; * 10 = (1/mm => 1/cm) 

muEffVec   = [gLight, gPivot, gDark];
muEffRatio = [muEffVec./muEffVec(2)];

simulatedRatio = [0.78072, 1, 1.361];

ratioErr = abs((simulatedRatio - muEffRatio)./simulatedRatio *100);

figure();
plot(simulatedRatio, muEffRatio, '-+');


%% 04-02-2022: 3D Scan 5cm depth: Simulations & Measurements
% close all
% clear all
% clc

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

depthAx = vars.grid.depthVec * 1e3;       % in mm
dx      = (depthAx(2) - depthAx(1)) / 10; % in cm

figure();
for i = 1:5
%     figure();
    phi(i,:,:) = reshape(res(i).phi, 9, []);
    subplot(2,3,i)
    imagesc(squeeze(phi(i,:,:)));
%     subplot(1,2,1)
    for j=1:9        
%         plot(phi(j,:)); hold on
    end
%     title(num2str(i));
    avgPhi(i, :) = mean(phi(i,2:end, :),2);
%     subplot(1,2,2)
%     plot(avgPhi(i,:))
end

%% Normalization
redPhi = phi(:,:,60:103);
figure();
for i=1:5
    subplot(2,3,i)
    imagesc(squeeze(redPhi(i,:,:)));
end

minMat = min(redPhi,[],3);
maxMat = max(redPhi,[],3);
spanMat = repmat(maxMat - minMat, 1,1,size(redPhi,3));

phiNorm = (redPhi - minMat) ./ spanMat;

figure();
for i=1:5
    subplot(2,3,i)
    imagesc(squeeze(phiNorm(i,:,:)));
end

figure();
for i=1:5
    subplot(2,3,i)
    imagesc(log(squeeze(phiNorm(i,:,:))));
end

%% Calibration
[~,I]  = max(phiNorm,[],3);
M = max(I(:));
refIdx = mode(I(:));
minIdx = 3;
maxIdx = M+3;
idxVec = minIdx:1:maxIdx;
noi = length(idxVec);
shiftMat = I - refIdx;

for i=1:5
    for j=1:9
        tmpIdx      = idxVec + shiftMat(i,j);
        phiNet(i,j,:) = flip(phiNorm(i,j,tmpIdx));
    end
end

phiAvg  = squeeze(mean(phiNet,2));
figure();
plot(phiAvg')

% This code is to align after averaging
% [~,I]  = max(phiAvg,[],2);
% M = min(I(:));
% refIdx = mode(I(:));
% minIdx = M-2;
% maxIdx = noi-2;
% idxVec = minIdx:1:maxIdx;
% noi = length(idxVec);
% shiftMat = I - refIdx;
% 
% for i=1:5
%     tmpIdx      = idxVec + shiftMat(i);
%     phiAlign(i,:) = phiAvg(i,tmpIdx);
% end

phiLog  = log(sqrt(phiAvg));

depthAxRed = depthAx(1:noi);

figure();
for i =1:5
    plot(depthAxRed, phiLog(i,:)); hold on
end
legend

figure();
for i =1:5
    plot(phiLog(i,:)); hold on
end
legend

%%
% figure();
% for i=1:5
%     subplot(2,3,i)
%     imagesc(squeeze(phiNet(i,:,:)));
% end
% 
% figure();
% for i=1:5
%     subplot(2,3,i)
%     imagesc(log(squeeze(phiNet(i,:,:))));
% end

% peak = 95;
% minIdx = 69;
% maxIdx = 98;
% idxsVec = minIdx:1:maxIdx;
% noi = length(idxsVec);
% peakVec = [95, 94, 95, 96, 93 ];
% shiftVec = peakVec - peak;
% 
% for i =1:5
%    tmpIdx = idxsVec + shiftVec(i);
%    phiNet(i,:) = flip(avgPhi(i,tmpIdx));
% end
% 
% depthAx = depthAx(1:noi);
% 
% minPhi =  min(phiNet, [], 2);
% maxPhi =  max(phiNet, [], 2);
% spanPhi = repmat(abs(maxPhi - minPhi), [1, noi]) ;
% phiNorm = (phiNet - minPhi) ./ spanPhi;
% phiLog  = log(sqrt(phiNorm));
% 
% figure();
% for i =1:5
%     plot(depthAx, phiLog(i,:)); hold on
% end
% legend


%% Calculate gradient
maxIdx =  12;
minIdx =  8;
slope = phiLog(:, minIdx:maxIdx);
muEffVec = abs(mean(diff(slope, 1, 2) / dx, 2));

muEffRatioVec = muEffVec/muEffVec(3);
gtMuEff = [0.7472    1.0627    1.5197    2.1960    3.2337];
gtMuEffRatio = gtMuEff / gtMuEff(3);

figure();
subplot(1,2,1)
plot (1:5, gtMuEffRatio); hold on
plot(1:5, muEffRatioVec);
legend("GT", "Meas")
title("Ratio")
subplot(1,2,2)
plot (1:5, gtMuEff); hold on
plot(1:5, muEffVec);
legend("GT", "Meas")
title("MuEff")

%% Extract Simulation

simDepthAx = resSim.depthVec;
dxSim = (simDepthAx(2) - simDepthAx(1)) / 10;

for i = 1:5
    phiSim(i,:) = resSim.resArr{i}.phiHRMid2In;
end

noi = length(simDepthAx);

minPhiSim =  min(phiSim, [], 2);
maxPhiSim =  max(phiSim, [], 2);
spanPhiSim = repmat(abs(maxPhiSim - minPhiSim), [1, noi]) ;
phiNormSim = (phiSim - minPhiSim) ./ spanPhiSim;
phiLogSim  = log(phiNormSim);

% figure();
% plot(phiLogSim');

maxIdx =  50;
minIdx =  1;
slopeSim = phiLogSim(:, minIdx:maxIdx);
muEffVecSim = abs(mean(diff(slopeSim, 1, 2) / dxSim, 2));

muEffRatioVecSim = muEffVecSim/muEffVecSim(3);

gtMuEff = [0.7472    1.0627    1.5197    2.1960    3.2337];
gtMuEffRatio = gtMuEff / gtMuEff(3);

figure()
for i = 1:5
    subplot(2,3,i)
    plot(simDepthAx, phiLogSim(i,:)); hold on
    plot(depthAxRed(1:end-5), phiLog(i,6:end))
    legend("Sim", "Meas", 'Location', 'southwest')
end

figure();
subplot(1,2,1)
plot (1:5, gtMuEffRatio); hold on
plot(1:5, muEffRatioVecSim);
plot(1:5, muEffRatioVec);
legend("GT", "Sim", "Meas")
title("Ratio")
subplot(1,2,2)
plot (1:5, gtMuEff); hold on
plot(1:5, muEffVecSim);
plot(1:5, muEffVec);
legend("GT", "Sim", "Meas")
title("MuEff")


%%
depthAx = s3D.grid.depthVec;
res = s3D.res;
phi3D = res.phi;

figure();
for i = 1:1
    phi(i,:,:) = reshape(res(i).phi, 9, []);
    subplot(2,3,i)
    imagesc(squeeze(phi(i,:,:)));
end

phi2 = squeeze(phi);

figure();
for i=1:9
    ax(i) = subplot(3,3,i);
    plot(phi2(i,:))
end
linkaxes(ax)


% Normalization
redPhi = phi(:,:,60:103);
figure();
for i=1:1
    subplot(2,3,i)
    imagesc(squeeze(redPhi(i,:,:)));
end

minMat = min(redPhi,[],3);
maxMat = max(redPhi,[],3);
spanMat = repmat(maxMat - minMat, 1,1,size(redPhi,3));

phiNorm = (redPhi - minMat) ./ spanMat;

figure();
for i=1:1
    subplot(2,3,i)
    imagesc(squeeze(phiNorm(i,:,:)));
end

figure();
for i=1:9
    ax(i) = subplot(3,3,i);
    plot(log(squeeze(phiNorm(:,i,:))));
end
linkaxes(ax)

figure();
for i=1:9
    ax(i) = subplot(3,3,i);
    plot(phi2(i,:))
end




%% Calibration
[~,I]  = max(phiNorm,[],3);
M = max(I(:));
refIdx = mode(I(:));
minIdx = 3;
maxIdx = M+3;
idxVec = minIdx:1:maxIdx;
noi = length(idxVec);
shiftMat = I - refIdx;

for i=1:1
    for j=1:9
        tmpIdx      = idxVec + shiftMat(i,j);
        phiNet(i,j,:) = flip(phiNorm(i,j,tmpIdx));
    end
end

depthAxRed = depthAx(1:noi)*1e3;

figure();
plot(depthAxRed, log(squeeze(phiNet))')


phiAvg  = squeeze(mean(phiNet,2));
figure();
plot(phiAvg)

% This code is to align after averaging
% [~,I]  = max(phiAvg,[],2);
% M = min(I(:));
% refIdx = mode(I(:));
% minIdx = M-2;
% maxIdx = noi-2;
% idxVec = minIdx:1:maxIdx;
% noi = length(idxVec);
% shiftMat = I - refIdx;
% 
% for i=1:5
%     tmpIdx      = idxVec + shiftMat(i);
%     phiAlign(i,:) = phiAvg(i,tmpIdx);
% end

phiLog  = log(sqrt(phiAvg'));

depthAxRed = depthAx(1:noi);

figure();
for i =1:1
    plot(depthAxRed, phiLog(i,:)); hold on
end
legend

figure();
for i =1:1
    plot(phiLog(i,:)); hold on
end
legend
