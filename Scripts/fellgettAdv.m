clear all;
close all;

naiveDir = 'D:\ResultsToKeep\AcoustoOptics-woHad-20mW\';
hadDir = 'D:\ResultsToKeep\AcoustoOptics-wHad-20mW\';

for i = 1:30
    filename = sprintf("%s%d-Results.mat", hadDir, i);
    load(filename);
    hadRes(i) = res;
    
    hadPhiCh(i,:,:,:) = res.phiCh;
    hadPhi(i,:) = res.phi;
end

for i = 1:30
    filename = sprintf("%s%d-Results.mat", naiveDir, i);
    load(filename);
    naiveRes(i) = res;
    
    naivePhiCh(i,:,:,:) = gather(res.phiCh);
    naivePhi(i,:) = gather(res.phi);
    
end

%%
hadPhiMean = mean(hadPhi, 1);
hadPhiStd  = std(hadPhi, 0, 1);
naivePhiMean = mean(naivePhi, 1);
naivePhiStd  = std(naivePhi, 0, 1);
 
SNRHad = (hadPhiMean(8)-  mean(hadPhiMean(11:end)))/ mean(hadPhiMean(11:end));
SNRNaive = (naivePhiMean(8) - mean(naivePhiMean(11:end)))  / mean(naivePhiMean(11:end));

factor = SNRHad / SNRNaive;

SNRHad2   = hadPhiMean(7).^2/hadPhiStd(7).^2;
SNRnaive2 = naivePhiMean(7).^2/naivePhiStd(7).^2;

hadPhiNormOffset   = mean(hadPhiMean(12:end));
naivePhiNormOffset = mean(naivePhiMean(12:end));

hadPhiNorm   = hadPhiMean   - hadPhiNormOffset;
naivePhiNorm = naivePhiMean - naivePhiNormOffset;

figure()
subplot(1,2,1)
errorbar(1:19, hadPhiMean, hadPhiStd)
title('hadamard')
subplot(1,2,2)
errorbar(1:19, naivePhiMean, naivePhiStd)
title('naive')

%%
hadScan = load('D:\ResultsToKeep\Scan-100mW-2msQ-2s-wHad\Results.mat');
hadVars = load('D:\ResultsToKeep\Scan-100mW-2msQ-2s-wHad\Vars.mat');

naiveScan = load('D:\ResultsToKeep\Scan-100mW-2msQ-2s-woHad\Results.mat');
naiveVars = load('D:\ResultsToKeep\Scan-100mW-2msQ-2s-woHad\Vars.mat');

%% Main Plain

zIdx = 8;
hadNorm   =  hadScan.res.phiRep(:,:,zIdx) - min(min(hadScan.res.phiRep(:,:,zIdx)));
hadNorm   = hadNorm/max(max(hadNorm));
hadNormLog = log10(hadNorm);
climsHad = [];
naiveNorm   =  naiveScan.res.phiRep(:,:,zIdx) - min(min(naiveScan.res.phiRep(:,:,zIdx)));
naiveNorm = naiveNorm/max(max(naiveNorm));
naiveNormLog = log10(naiveNorm);
climsNaive = [];


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
imagesc('XData', hadVars.curVars.stages.yVec,...
        'YData', hadVars.curVars.stages.xVec,...
        'CData', hadNormLog);
colorbar;
title('Hadamard (log)');
axis equal
axis tight

figure()
%Naive
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