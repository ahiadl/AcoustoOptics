clear all;
close all;
clc

%% Preliminary Results - 2D scan 4Cycles in pulse
hadScan = load('D:\ResultsToKeep\Scan-100mW-2msQ-2s-wHad\Results.mat');
hadVars = load('D:\ResultsToKeep\Scan-100mW-2msQ-2s-wHad\Vars.mat');

naiveScan = load('D:\ResultsToKeep\Scan-100mW-2msQ-2s-woHad\Results.mat');
naiveVars = load('D:\ResultsToKeep\Scan-100mW-2msQ-2s-woHad\Vars.mat');

zIdx = 8;
hadNorm   =  hadScan.res.phiRep(:,:,zIdx) - min(min(hadScan.res.phiRep(:,:,zIdx)));
hadNorm   = hadNorm/max(max(hadNorm));
hadNormLog = log10(hadNorm);
climsHad = [];
naiveNorm   =  naiveScan.res.phiRep(:,:,zIdx) - min(min(naiveScan.res.phiRep(:,:,zIdx)));
naiveNorm = naiveNorm/max(max(naiveNorm));
naiveNormLog = log10(naiveNorm);
climsNaive = [];

%% Figures for publication
% axes Pos = [left, bottom, width,height] = [0.13,0.58,0.67,0.34]
left = 0.13;
bottom = 0.58;
width = 0.67;
height = 0.34;
vGap = 0.07;

labelsFontSize = 24;
axesFontSize = 22;

hFig = figure();
set(hFig, 'Position', [119,42,936,954]);

% titAx = subplot(1,1,1);
% hCBTit = title("Normalized Intensity");
% a.Title.Rotation = 90;

ax = subplot(2,2,2);
set(ax, 'Position' , [left, bottom, width,height])
imagesc('XData', naiveVars.curVars.stages.yVec - mean(naiveVars.curVars.stages.yVec),...
        'YData', -1*(naiveVars.curVars.stages.xVec - max(naiveVars.curVars.stages.xVec)),...
        'CData', naiveNorm);
hCBar = colorbar;
hCBarTit = ylabel(hCBar, 'Normalized Intensity');
hCBarTit.FontSize = labelsFontSize;
set(ax, 'FontSize', axesFontSize)
hLabel = ylabel( 'Y[mm]' );
set(hLabel, 'FontSize', labelsFontSize);
axis equal
axis tight


bottom = bottom - height - vGap;
ax = subplot(2,2,4);
set(ax, 'Position' , [left, bottom, width,height])
imagesc('XData', hadVars.curVars.stages.yVec - mean(hadVars.curVars.stages.yVec),...
        'YData', -1*(hadVars.curVars.stages.xVec - max(hadVars.curVars.stages.xVec)),...
        'CData', hadNorm);
hCBar = colorbar;
hCBarTit = ylabel(hCBar, 'Normalized Intensity');
hCBarTit.FontSize = labelsFontSize;
set(ax, 'FontSize', axesFontSize)
hLabel = xlabel( 'X[mm]' );
set(hLabel, 'FontSize', labelsFontSize);
hLabel = ylabel( 'Y[mm]' );
set(hLabel, 'FontSize', labelsFontSize);
axis equal
axis tight

%% Figures for publication
% axes Pos = [left, bottom, width,height] = [0.13,0.58,0.67,0.34]
left = 0.13;
bottom = 0.58;
width = 0.67;
height = 0.34;
vGap = 0.09;

labelsFontSize = 24;
axesFontSize = 22;

hFig = figure();
set(hFig, 'Position', [119,42,936,954]);

% titAx = subplot(1,1,1);
% hCBTit = title("Normalized Intensity");
% a.Title.Rotation = 90;

ax = subplot(2,1,2);
set(ax, 'Position' , [left, bottom, width,height])
imagesc('XData', naiveVars.curVars.stages.yVec - mean(naiveVars.curVars.stages.yVec),...
        'YData', -1*(naiveVars.curVars.stages.xVec - max(naiveVars.curVars.stages.xVec)),...
        'CData', naiveNormLog);
hCBar = colorbar;
hCBarTit = ylabel(hCBar, {'Log Normalized' ; 'Intensity'});
hCBarTit.FontSize = labelsFontSize;
set(ax, 'FontSize', axesFontSize)
hLabel = xlabel( 'X[mm]' );
set(hLabel, 'FontSize', labelsFontSize);
hLabel = ylabel( 'Y[mm]' );
set(hLabel, 'FontSize', labelsFontSize);
axis equal
axis tight


bottom = bottom - height - vGap;
ax = subplot(2,1,2);
set(ax, 'Position' , [left, bottom, width,height])
imagesc('XData', hadVars.curVars.stages.yVec - mean(hadVars.curVars.stages.yVec),...
        'YData', -1*(hadVars.curVars.stages.xVec - max(hadVars.curVars.stages.xVec)),...
        'CData', hadNormLog);
hCBar = colorbar;
hCBarTit = ylabel(hCBar, {'Log Normalized' ; 'Intensity'});
hCBarTit.FontSize = labelsFontSize;
set(ax, 'FontSize', axesFontSize)
hLabel = xlabel( 'X[mm]' );
set(hLabel, 'FontSize', labelsFontSize);
hLabel = ylabel( 'Y[mm]' );
set(hLabel, 'FontSize', labelsFontSize);
axis equal
axis tight

% setupIm = imread('..\images\AOSetupIllustration.jpg');
% subplot(1,2,1)
% imshow(setupIm)

%% Main Plain

%Hadamard
% Original Scale
figure();
subplot(1,2,1)
imagesc('XData', hadVars.curVars.stages.yVec,...
        'YData', hadVars.curVars.stages.xVec,...
        'CData', hadScan.res.phiRep(:,:,zIdx));
colorbar;
title('Hadamard');
axis equal
axis tight

subplot(1,2,2)
imagesc('XData', naiveVars.curVars.stages.yVec,...
        'YData', naiveVars.curVars.stages.xVec,...
        'CData', naiveScan.res.phiRep(:,:,zIdx));
colorbar;
title('Naive');
axis equal
axis tight

% subplot(1,2,2)
% imagesc('XData', hadVars.curVars.stages.yVec,...
%         'YData', hadVars.curVars.stages.xVec,...
%         'CData', hadNormLog);
% colorbar;
% title('Hadamard (log)');
% axis equal
% axis tight

figure()
%Naive
% Original Scale
subplot(1,2,1)



% Logarithmic Scale (Normalized)
subplot(1,2,2)
imagesc('XData', naiveVars.curVars.stages.yVec,...
        'YData', naiveVars.curVars.stages.xVec,...
        'CData', naiveNormLog);
colorbar;
title('Naive (log)');
axis equal
axis tight



%% Secondary Plane
zIdx = 7;

hadNormSec   =  hadScan.res.phiRep(:,:,zIdx) - min(min(hadScan.res.phiRep(:,:,zIdx)));
hadNormSec   = hadNormSec/max(max(hadNormSec));
hadNormLogSec = log10(hadNormSec);
climsHad = [];
naiveNormSec   =  naiveScan.res.phiRep(:,:,zIdx) - min(min(naiveScan.res.phiRep(:,:,zIdx)));
naiveNormSec = naiveNormSec/max(max(naiveNormSec));
naiveNormLogSec = log10(naiveNormSec);
climsNaive = [];

figure();
% Original Scale
subplot(1,2,1)
imagesc('XData', hadVars.curVars.stages.yVec,...
        'YData', hadVars.curVars.stages.xVec,...
        'CData', hadScan.res.phiRep(:,:,zIdx));
colorbar;
title('Hadamard');
axis equal
axis tight

% Logarithmic Scale (Normalized)
subplot(1,2,2)
imagesc('XData', hadVars.curVars.stages.yVec,...
        'YData', hadVars.curVars.stages.xVec,...
        'CData', hadNormLogSec);
colorbar;
title('Hadamard (log)');
axis equal
axis tight


figure();
% Original Scale
subplot(1,2,1)
imagesc('XData', naiveVars.curVars.stages.yVec,...
        'YData', naiveVars.curVars.stages.xVec,...
        'CData', naiveScan.res.phiRep(:,:,zIdx));
colorbar;
title('Naive');
axis equal
axis tight

% Logarithmic Scale (Normalized)
subplot(1,2,2)
imagesc('XData', naiveVars.curVars.stages.yVec,...
        'YData', naiveVars.curVars.stages.xVec,...
        'CData', naiveNormLogSec);
colorbar;
title('Naive (log)');
axis equal
axis tight



%% Fellgett's Advantage - Dynamic Range
cyc = [4];
orders = {[7, 11, 15, 19]};

for c = 1:length(cyc)
    for i = 1:length(orders{c})
        for k=1:2
            if k == 1        
                dirName = sprintf("D:/Results/AcoustoOptics-%dcyc-woHad-%dn-20mW/", cyc(c), orders{c}(i));
                cyc4VarNaive(c,i) = load(sprintf("%sVars.mat", dirName));
                fprintf("Loading: Cycles: %d, Order: %d, Hadamard: %d\n", cyc(c), orders{c}(i), cyc4VarNaive(c,i).vars.algo.uVars.useHadamard);
                for j =1:30
                    filename = sprintf("%s%d-Results",dirName, i);
                    cyc4Naive(c,i,j) = load(filename);
                    cyc4Naive(c,i,j)= load(filename);
                    cyc4PhiNaive{c,i}(j,:) = gather(cyc4Naive(c,i,j).res.phi);
                end
            else
                dirName = sprintf("D:/Results/AcoustoOptics-%dcyc-wHad-%dn-20mW/", cyc(c),  orders{c}(i));
                cyc4VarHad(c,i) = load(sprintf("%sVars.mat", dirName));
                fprintf("Loading: Cycles: %d, Order: %d, Hadamard: %d\n", cyc(c), orders{c}(i), cyc4VarHad(c,i).vars.algo.uVars.useHadamard);
                for j =1:30
                    filename = sprintf("%s%d-Results",dirName, i);
                     cyc4Had(c,i,j)= load(filename);
                     cyc4PhiHad{c,i}(j,:) = gather(cyc4Had(c,i,j).res.phi);
                end
            end
        end
    end
end

hFig = figure();
ax = axes();
for c = 1:length(cyc)
    for i=1:length(orders{c})
        meanSigHad = mean(cyc4PhiHad{c,i}, 1);
        meanSigNaive = mean(cyc4PhiNaive{c,i}, 1);
        figure()
        plot(meanSigHad); hold on;
        plot(meanSigNaive);
        
        hadAVGMax(c,i) = mean(max(cyc4PhiHad{c,i}, [], 2));
        hadAVGMin(c,i) = mean(mean(mink(cyc4PhiHad{c,i}, 3, 2),2));
        hadAVGDynamicRange =  abs(hadAVGMax(c,i) - hadAVGMin(c,i))/hadAVGMin(c,i);
        
        naiveAVGMax(c,i) = mean(max(cyc4PhiNaive{c,i}, [], 2));
        naiveAVGMin(c,i) = mean(mean(mink(cyc4PhiNaive{c,i}, 3, 2), 2));
        naiveAVGDynamicRange =  abs(naiveAVGMax(c,i) - naiveAVGMin(c,i))/naiveAVGMin(c,i);
        
        dynamicRangeFactor{c}(i) = hadAVGDynamicRange/naiveAVGDynamicRange;
    
    end
    set(0, 'currentfigure',  hFig);
    subplot(1,length(cyc), c);
    plot(orders{c}, sqrt(orders{c}), '-+'); hold on;
    plot(orders{c}, dynamicRangeFactor{c}, '-v');hold off
    legend('Theory', 'Measurement');
    title(sprintf("%d Cycles In Pulse", cyc(c)));
    xlabel('Hadamard Order')
    if(c==1)
        ylabel('Dynamic Range Factor')
    end
end