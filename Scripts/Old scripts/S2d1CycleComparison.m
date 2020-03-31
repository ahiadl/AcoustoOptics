close all
clear all
clc;
%%
% singlePulseRes = load('D:\Results\17-Dec-2019 13-51-41-SinglePulse-1CycRes-full2dScan\Results.mat');
% multiRes = load('D:\Results\18-Dec-2019 09-04-45-Multiplexed-1CycRes-fulls2dScan\Results.mat');

% res1.phi = singlePulseRes.res.phi;
% res1.phiStd = singlePulseRes.res.phiStd; 
% save('D:\Results\17-Dec-2019 13-51-41-SinglePulse-1CycRes-full2dScan\phi.mat', 'res1', '-v7.3')

% res2.phi = multiRes.res.phi;
% res2.phiStd = multiRes.res.phiStd; 
% save('D:\Results\18-Dec-2019 09-04-45-Multiplexed-1CycRes-fulls2dScan\phi.mat', 'res2', '-v7.3')
%%
% load Results
singlePulseRes = load('D:\Results\17-Dec-2019 13-51-41-SinglePulse-1CycRes-full2dScan\phi.mat');
multiRes = load('D:\Results\18-Dec-2019 09-04-45-Multiplexed-1CycRes-fulls2dScan\phi.mat');
singlePulseVars = load('D:\Results\17-Dec-2019 13-51-41-SinglePulse-1CycRes-full2dScan\Vars.mat');
multiVars = load('D:\Results\18-Dec-2019 09-04-45-Multiplexed-1CycRes-fulls2dScan\Vars.mat');
%% Rearrrange for display
yAxis = multiVars.curVars.stages.firstVec;
xAxis = multiVars.curVars.stages.secondVec;

yAxisNorm = yAxis - yAxis(ceil(length(yAxis)/2));
xAxisNorm = -(xAxis - xAxis(1));

zIdx = 29;
% 
% xS1 = multiVars.curVars.stages.firstIdx(yAxisNorm == 28.5);
% yS1 = multiVars.curVars.stages.secondIdx(xAxisNorm == 12.5);
% xS2 = multiVars.curVars.stages.firstIdx(yAxisNorm == 26.5);
% yS2 = multiVars.curVars.stages.secondIdx(xAxisNorm == -9);
% 
% xM1 = multiVars.curVars.stages.firstIdx(yAxisNorm == 6); 
% yM1 = multiVars.curVars.stages.secondIdx(xAxisNorm == 26);
% xM2 = multiVars.curVars.stages.firstIdx(yAxisNorm == -14);  
% yM2 = multiVars.curVars.stages.secondIdx(xAxisNorm == 30);

imgSPOrig = singlePulseRes.res1.phi(:,:,zIdx);
imgMulOrig = multiRes.res2.phi(:,:,zIdx);

imgSP = imgSPOrig;
% imgSP(xS1, yS1) = imgSP(xS1-1, yS1+1);
% imgSP(xS2, yS2) = imgSP(xS2+1, yS2+1);
% 
imgMul = imgMulOrig;
% imgMul(xM1, yM1) = imgMul(xM1+1, yM1+1);
% imgMul(xM2, yM2) = imgMul(xM2+1, yM2+1);


imgSPNorm = (imgSP - min(min(imgSP))) /  abs(max(max(imgSP)) - min(min(imgSP)) );
imgSPLog = log10(imgSPNorm); 

imgMulNorm = (imgMul - min(min(imgMul))) /  abs(max(max(imgMul)) - min(min(imgMul)) );
imgMulLog = log10(imgMulNorm);

%% Display
figure();
subplot(1,2,1)
imagesc('XData', yAxisNorm, ...
        'YData', xAxisNorm, ...
        'CData', imgSPLog);
axis equal
axis tight
title('Single Pulse')
xlabel('Y[mm]')
ylabel('X[mm]')
colorbar
subplot(1,2,2)
imagesc('XData', yAxisNorm, ...
        'YData', xAxisNorm, ...
        'CData', imgMulLog);
axis equal
axis tight
title('Multiplexed')
xlabel('Y[mm]')
ylabel('X[mm]')
colorbar

figure();
subplot(1,2,1)
imagesc('XData', yAxisNorm, ...
        'YData', xAxisNorm, ...
        'CData', imgSPNorm);
axis equal
axis tight
title('Single Pulse')
xlabel('Y[mm]')
ylabel('X[mm]')
colorbar
subplot(1,2,2)
imagesc('XData', yAxisNorm, ...
        'YData', xAxisNorm, ...
        'CData', imgMulNorm);
axis equal
axis tight
title('Multiplexed')
xlabel('Y[mm]')
ylabel('X[mm]')
colorbar

%% Display For Paper:

% hFig = figure();
% set(hFig, 'Position', [680,42,300,850]


