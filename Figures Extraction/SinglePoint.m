close all;
% clear all;
clc;
  
resDir   = 'C:\Users\sahiadl\OneDrive - Technion\Masters\AcoustoOpticSystem\Code\Measurements\'; 
simDir1  = 'Results-04-Apr-2019 11-41-31-50mWMotorized-singlePos-10rep\';
simDir2  = 'Results-04-Apr-2019 16-10-47-74mWNoMotorized-singlePos-10rep\';
fileName = 'paramsAndTotalResults.mat';

% load ('C:\Users\sahiadl.EED\OneDrive - Technion\Masters\AcoustoOpticSystem\Code\Measurements\Results-24-Feb-2019 11-11-37-Single Position - with Speckle Shifting\paramsAndTotalResults.mat');
% load ('D:\Results\Results-04-Apr-2019 11-41-31-50mWMotorized-singlePos-10rep\paramsAndTotalResults.mat');
load([resDir, simDir1,fileName]);
titleFontSize = 18;
labelFontSize = 18;
legendFontSize = 18;
axesFontSize = 18;

timeFrameVec = 2:2:40;
timeFrameStep = 2;

xlims = [timeFrameVec(1)-timeFrameStep, timeFrameVec(end)+timeFrameStep];

zCoor = 10;
% figure()
% stem(timeFrameVec, I(10, :, 1))

figure()
subplot
stdIVec=[];
for i =1:length(timeFrameVec)
    stdIVec(i) = std(I(zCoor,i,:));
end
avgIVec = mean(I(zCoor,:,:), 3);

ylimI = [0.5*min(avgIVec - stdIVec), 1.1*max(avgIVec + stdIVec)];
errorbar(timeFrameVec, avgIVec, stdIVec);
set(gca, 'FontSize', axesFontSize);
title({'Single Point Power Spectrum vs.'; 'Duration Time Of Measurement'}, 'FontSize', titleFontSize );
xlabel('Time Fram Duration [s]','FontSize', labelFontSize );
ylabel('Average Power Spectrum', 'FontSize', labelFontSize );
xlim(xlims)
ylim(ylimI)

for j = 1:4
    for i =1:length(timeFrameVec)
        stdVMat(j,i) = std(V(zCoor,j,i,:));
    end
    avgVMat(j,:) = squeeze(mean(V(zCoor,j,:,:), 4));
end

ylimVCH = [0.5*min(min(avgVMat-stdVMat)), 1.1*max(max(avgVMat+stdVMat))];


figure()
for j = 1:4
    subplot(2,2,j)
%     for i =1:length(timeFrameVec)
%         stdVec(i) = std(V(zCoor,j,i,:));
%     end
%     avgVVec = squeeze(mean(V(zCoor,j,:,:), 4))
    errorbar(timeFrameVec, avgVMat(j,:), stdVMat(j,:));
    hold on
    set(gca, 'FontSize', axesFontSize);
    title({'Single Point Power Spectrum vs. Duration Time Of Measurement'; ['Channel ' num2str(j)]}, 'FontSize', titleFontSize );
    xlabel('Time Fram Duration [s]','FontSize', labelFontSize );
    ylabel('Average Power Spectrum', 'FontSize', labelFontSize );
    xlim(xlims)
    ylim(ylimVCH)
end



% ---------------------------------------------------------------------------------%


% load ('C:\Users\sahiadl.EED\OneDrive - Technion\Masters\AcoustoOpticSystem\Code\Measurements\Results-24-Feb-2019 12-37-48-Single Position - no Speckle Shifting\paramsAndTotalResults.mat');
% load ('D:\Results\Results-04-Apr-2019 16-10-47-74mWNoMotorized-singlePos-10rep\paramsAndTotalResults.mat');
load([resDir, simDir2, fileName]);

% I - (z, T, r)
% V - (z, ch, T, r)

titleFontSize = 18;
labelFontSize = 18;
legendFontSize = 18;
axesFontSize = 18;

timeFrameVec = 2:2:40;
timeFrameStep = 2;
xlims = [timeFrameVec(1)-timeFrameStep, timeFrameVec(end)+timeFrameStep];

zCoor = 10;
% figure()
% stem(timeFrameVec, I(10, :, 1))

figure()
subplot
stdIVec=[];
for i =1:length(timeFrameVec)
    stdIVec(i) = std(I(zCoor,i,:));
end
avgIVec = mean(I(zCoor,:,:), 3);

ylimI = [0.5*min(avgIVec - stdIVec), 1.1*max(avgIVec + stdIVec)];
errorbar(timeFrameVec, avgIVec, stdIVec);
set(gca, 'FontSize', axesFontSize);
title({'Single Point Power Spectrum vs.'; 'Duration Time Of Measurement'}, 'FontSize', titleFontSize );
xlabel('Time Fram Duration [s]','FontSize', labelFontSize );
ylabel('Average Power Spectrum', 'FontSize', labelFontSize );
xlim(xlims)
ylim(ylimI)

for j = 1:4
    for i =1:length(timeFrameVec)
        stdVMat(j,i) = std(V(zCoor,j,i,:));
    end
    avgVMat(j,:) = squeeze(mean(V(zCoor,j,:,:), 4));
end

ylimVCH = [0.5*min(min(avgVMat-stdVMat)), 1.1*max(max(avgVMat+stdVMat))];


figure()
for j = 1:4
    subplot(2,2,j)
%     for i =1:length(timeFrameVec)
%         stdVec(i) = std(V(zCoor,j,i,:));
%     end
%     avgVVec = squeeze(mean(V(zCoor,j,:,:), 4))
    errorbar(timeFrameVec, avgVMat(j,:), stdVMat(j,:));
    hold on
    set(gca, 'FontSize', axesFontSize);
    title({'Single Point Power Spectrum vs. Duration Time Of Measurement'; ['Channel ' num2str(j)]}, 'FontSize', titleFontSize );
    xlabel('Time Fram Duration [s]','FontSize', labelFontSize );
    ylabel('Average Power Spectrum', 'FontSize', labelFontSize );
    xlim(xlims)
    ylim(ylimVCH)
end

%-------Comparison--------%
% noShift = load ('D:\Results\Results-04-Apr-2019 16-10-47-74mWNoMotorized-singlePos-10rep\paramsAndTotalResults.mat');
% Shift = load ('D:\Results\Results-04-Apr-2019 11-41-31-50mWMotorized-singlePos-10rep\paramsAndTotalResults.mat');

noShift = load([resDir, simDir2,fileName]);
Shift = load([resDir, simDir1,fileName]);

figure()
subplot
stdIVec=[];
for i =1:length(timeFrameVec)
    stdIVecNoShift(i) = std(noShift.I(zCoor,i,:));
    stdIVecShift(i) = std(Shift.I(zCoor,i,:));
end
avgIVecNoShift = mean(noShift.I(zCoor,:,:), 3);
avgIVecShift = mean(Shift.I(zCoor,:,:), 3);

ylimINoShift = [0.5*min(avgIVecNoShift - stdIVecNoShift), 1.1*max(avgIVecNoShift + stdIVecNoShift)];
ylimIShift = [0.5*min(avgIVecShift - stdIVecShift), 1.1*max(avgIVecShift + stdIVecShift)];

ylimI = [min([ylimINoShift(1) ylimIShift(1)]), max([ylimINoShift(2) ylimIShift(2)])];
errorbar(timeFrameVec, avgIVecNoShift, stdIVecNoShift);hold on;
errorbar(timeFrameVec, avgIVecShift, stdIVecShift);
set(gca,'yscale','log')
set(gca, 'FontSize', axesFontSize);
title({'Single Point Power Spectrum vs.'; 'Duration Time Of Measurement'}, 'FontSize', titleFontSize );
xlabel('Time Fram Duration [s]','FontSize', labelFontSize );
ylabel('Average Power Spectrum', 'FontSize', labelFontSize );
xlim(xlims)
ylim(ylimI)
legend('No Shift', 'Shift')

% %-------- FFT ------%
% 
% fBarNoSpeckle = load('C:\Users\sahiadl.EED\OneDrive - Technion\Masters\AcoustoOpticSystem\Code\Measurements\fBarNoShift.mat');
% fBarNoSpeckle = fBarNoSpeckle.A;
% fBarWSpeckle = load('C:\Users\sahiadl.EED\OneDrive - Technion\Masters\AcoustoOpticSystem\Code\Measurements\fBar.mat');
% fBarWSpeckle = fBarWSpeckle.A;
% fftResNoSpeckle = load('C:\Users\sahiadl.EED\OneDrive - Technion\Masters\AcoustoOpticSystem\Code\Measurements\fftResNoShift.mat');
% fftResNoSpeckle = fftResNoSpeckle.A;
% fftResWSpeckle = load('C:\Users\sahiadl.EED\OneDrive - Technion\Masters\AcoustoOpticSystem\Code\Measurements\fftRes.mat');
% fftResWSpeckle = fftResWSpeckle.A;
% 
% figure()
% subplot(1,2,1)
% plot(fBarNoSpeckle, abs(fftResNoSpeckle(11,:,1)))
% set(gca, 'FontSize', 18)
% title({"FFT Around 1.25[MHz]"; "no Speckle Shift"})
% 
% subplot(1,2,2)
% plot(fBarWSpeckle, abs(fftResWSpeckle(11,:,1)))
% set(gca, 'FontSize', 18)
% title({"FFT Around 1.25[MHz]";"With Speckle Shift"})
% 
% figure()
% subplot(1,2,1)
% plot(fBarNoSpeckle, squeeze(abs(fftResNoSpeckle(11,:,:))))
% set(gca, 'FontSize', 18)
% title({"FFT Around 1.25[MHz]";"no Speckle Shift"})
% legend('Ch 1', 'Ch 2', 'Ch 3', 'Ch 4')
% 
% subplot(1,2,2)
% plot(fBarWSpeckle, squeeze(abs(fftResWSpeckle(11,:,:))))
% set(gca, 'FontSize', 18)
% title({"FFT Around 1.25[MHz]";"With Speckle Shift"})
% legend('Ch 1', 'Ch 2', 'Ch 3', 'Ch 4')
