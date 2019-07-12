close all;
clear all;
clc;

twoSec = load ('D:\Results\Results-11-Jun-2019 13-04-11-Scan-100mW-500Hz-20VPV-0.5mm\paramsAndTotalResults.mat');

% User parameters:

fourSec = load('D:\Results\Results-11-Jun-2019 15-38-22-Scan-100mW-500Hz-20VpV-0.5mm-4s\paramsAndTotalResults.mat');


upLim =  max(  [max(max(max(fourSec.I))), max(max(max(twoSec.I)))]);
dLim =  min(  [min(min(min(fourSec.I))), min(min(min(twoSec.I)))]);
limits = [dLim, upLim];

zCoor = 10;

axesFontSize = 18;
titleFontSize = 18;
labelFontSize = 18;

depthX = abs(fourSec.params.Stages.Xpos(end)-fourSec.params.Stages.Xpos(1));
lX =length(fourSec.params.Stages.Xpos);
[M, idx] = min(abs(twoSec.params.Stages.Xpos-(twoSec.params.Stages.Xpos(1)-depthX)));

figure();
imagesc('XData', fourSec.params.Stages.Ypos, 'YData', fourSec.params.Stages.Xpos, ...
        'CData', permute(fourSec.I(:,:,zCoor), [2,1,3]));
set(gca, 'FontSize', axesFontSize);
title("4 Seconds of Sampling", 'FontSize', titleFontSize)
caxis(limits)
axis equal
axis tight
colorbar;

zCoor = 10;

figure()
imagesc('XData', twoSec.params.Stages.Ypos, 'YData', twoSec.params.Stages.Xpos(1:idx), ...
        'CData', permute(twoSec.I(:,1:idx,zCoor), [2,1,3]));
set(gca, 'FontSize', axesFontSize);
title("2 Seconds Of Sampling")
caxis(limits)
colorbar;
axis equal
axis tight
% 
% figure()
% imagesc('XData', twoSec.params.Stages.Ypos, 'YData', twoSec.params.Stages.Xpos(1:idx), ...
%         'CData', permute(twoSec.I(:,1:idx,zCoor), [2,1,3]));
% set(gca, 'FontSize', axesFontSize);
% title("With Speckle Shift - Scaled Colors")
% % caxis(limits)
% colorbar;
% axis equal
% axis tight

% figure()
% imagesc('XData', wSpeckleShiftNarrow.params.Stages.Ypos, 'YData', wSpeckleShiftNarrow.params.Stages.Xpos, ...
%         'CData', permute(wSpeckleShiftNarrow.I(:,:,zCoor), [2,1,3]));
% set(gca, 'FontSize', axesFontSize);
% title("With Speckle Shift - Fine Scan")
% colorbar;
% axis equal;
% axis tight;
