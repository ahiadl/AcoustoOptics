%% Number of Channels

dirName = "D:\ResultsToKeep\SNR\Channels";
measName = "Live-%sCH-2s";

ch = [1, 2, 4, 8, 16];


for i = 1:length(ch)
    measDir = sprintf("%s/Live-%dCH-2s/LiveAOResults", dirName, ch(i));
    for j = 1:30
        fileName = sprintf("%s/AO-Results-%d.mat", measDir, j);
        curRes = load(fileName);
        
        SNR(i,j)    = curRes.SNR.val;
        SNRRaw(i,j) = curRes.SNR.valRaw;
    end
end

avgSNR = mean(SNR,2);
stdSNR = std(SNR,0,2);

avgSNRRaw = mean(SNRRaw,2);
stdSNRRaw = std(SNRRaw,0,2);

% figure()
% errorbar(ch, avgSNR./avgSNR(1), stdSNR); hold on
% plot(ch, sqrt(ch), '-+')

hFig1 = figure();
ax = subplot(1,2,1)
errorbar(ch, avgSNRRaw./avgSNRRaw(1), stdSNRRaw/avgSNRRaw(1)); hold on
plot(ch, sqrt(ch), '-+')
xlim([0, 17])
xlabel("Number Of Channels")
ylabel("SNR Fold")
hleg = legend("Measurement", "theory");
set(hleg, 'Location', 'northwest')
set(ax, 'FontSize', 14)

%% Time Of Sampling
dirName = "D:\ResultsToKeep\SNR\Time";
measName = "Live-16CH-%ds";

t = [1, 2, 4, 8, 16];

for i = 1:length(t)
    measDir = sprintf("%s/Live-16CH-%ds/LiveAOResults", dirName, t(i));
    for j = 1:30
        fileName = sprintf("%s/AO-Results-%d.mat", measDir, j);
        curRes = load(fileName);
        
        SNR(i,j)    = curRes.SNR.val;
        SNRRaw(i,j) = curRes.SNR.valRaw;
        peakVal(i,j) = curRes.SNR.peak;
        peakValRaw(i,j) = curRes.SNR.peakRaw;
    end
end

avgSNR = mean(SNR,2);
stdSNR = std(SNR,0,2);

avgSNRRaw = mean(SNRRaw, 2);
stdSNRRaw = std(SNRRaw, 0, 2);

avgPeak = mean(peakVal, 2);
stdPeak = std(peakVal,0,2);

% figure()
% errorbar(t, avgSNR./avgSNR(1), stdSNR/stdSNR(1)); hold on
% plot(t, sqrt(t), '-+')

ax = subplot(1,2,2)
errorbar(t, avgSNRRaw./avgSNRRaw(1), stdSNRRaw/avgSNRRaw(1)); hold on
plot(t, sqrt(t), '-+')
xlim([0, 17])
xlabel("Time of Sampling [s]")
hleg = legend("Measurement", "theory");
set(hleg, 'Location', 'northwest')
set(ax, 'FontSize', 14)
% figure()
% errorbar(t, avgPeak, stdPeak); hold on

%% No Quants No Wl
dirName = "D:\ResultsToKeep\SNR\noQuants-NoWL";
measName = "Live-noWL-noQuants-%ds";

t = [1, 2, 4, 8, 16];

for i = 1:length(t)
    measDir = sprintf("%s/Live-noWL-noQuants-%ds/LiveAOResults", dirName, t(i));
    for j = 1:30
        fileName = sprintf("%s/AO-Results-%d.mat", measDir, j);
        curRes = load(fileName);
        
        SNR(i,j)     = curRes.SNR.val;
        SNRRaw(i,j)  = curRes.SNR.valRaw;
        peakVal(i,j) = curRes.SNR.peak;
        peakValRaw(i,j) = curRes.SNR.peakRaw;
    end
end

avgSNR = mean(SNR,2);
stdSNR = std(SNR,0,2);

avgSNRRaw = mean(SNRRaw,2);
stdSNRRaw = std(SNRRaw,0,2);

avgPeak = mean(peakVal,2);
stdPeak = std(peakVal,0,2);

figure()
errorbar(t, avgSNR./avgSNR(1), stdSNR/avgSNR(1)); hold on
plot(t, sqrt(t), '-+')
xlim([0,17])
ylim([0, 4.5])
title("SNR Fold vs theory")

figure()
errorbar(t, avgSNRRaw./avgSNRRaw(1), stdSNRRaw/avgSNRRaw(1)); hold on
plot(t, sqrt(t), '-+')
xlim([0,17])
ylim([0, 4.5])
title("SNR Fold (Raw) vs. theory")

figure()
errorbar(t, avgPeak, stdPeak); hold on
title("Peak Stability")
xlim([0,17])

%% No WL Scan
dirName = "D:\ResultsToKeep\SNR\NoWL";
measName = "Live-noWL-noQuants-%ds";

t = [1, 2, 4, 8, 16];

for i = 1:length(t)
    measDir = sprintf("%s/Live-noWL-%ds/LiveAOResults", dirName, t(i));
    for j = 1:30
        fileName = sprintf("%s/AO-Results-%d.mat", measDir, j);
        curRes = load(fileName);
        
        SNR(i,j)     = curRes.SNR.val;
        SNRRaw(i,j)  = curRes.SNR.valRaw;
        peakVal(i,j) = curRes.SNR.peak;
        peakValRaw(i,j) = curRes.SNR.peakRaw;
    end
end

avgSNR = mean(SNR,2);
stdSNR = std(SNR,0,2);

avgSNRRaw = mean(SNRRaw,2);
stdSNRRaw = std(SNRRaw,0,2);

avgPeak = mean(peakVal,2);
stdPeak = std(peakVal,0,2);

figure()
errorbar(t, avgSNR./avgSNR(1), stdSNR/stdSNR(1)); hold on
plot(t, sqrt(t), '-+')

figure()
errorbar(t, avgSNRRaw./avgSNRRaw(1), stdSNRRaw/stdSNRRaw(1)); hold on
plot(t, sqrt(t), '-+')

figure()
errorbar(t, avgPeak, stdPeak); hold on
