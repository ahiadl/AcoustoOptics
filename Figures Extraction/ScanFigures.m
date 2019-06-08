close all;
clear all;
clc;

wSpeckleShift = load ('C:\Users\sahiadl.EED\OneDrive - Technion\Masters\AcoustoOpticSystem\Code\Measurements\Results-20-Feb-2019 16-33-08-Scan-SpeckleShifting-WideArea\paramsAndTotalResults.mat');

% User parameters:

noSpeckleShift = load('C:\Users\sahiadl.EED\OneDrive - Technion\Masters\AcoustoOpticSystem\Code\Results\Results-17-Feb-2019 16-05-20.mat');

wSpeckleShiftNarrow = load('C:\Users\sahiadl.EED\OneDrive - Technion\Masters\AcoustoOpticSystem\Code\Measurements\Results-21-Feb-2019 13-51-42-ShiftingSpeckle-NarrowArea\paramsAndTotalResults.mat');

upLim =  max(  [max(max(max(noSpeckleShift.I))), max(max(max(wSpeckleShift.I)))]);
dLim =  min(  [min(min(min(noSpeckleShift.I))), min(min(min(wSpeckleShift.I)))]);
limits = [dLim, upLim];

zCoor = 10;

axesFontSize = 18;
titleFontSize = 18;
labelFontSize = 18;

depthX = abs(noSpeckleShift.params.Stages.Xpos(end)-noSpeckleShift.params.Stages.Xpos(1));
lX =length(noSpeckleShift.params.Stages.Xpos);
[M, idx] = min(abs(wSpeckleShift.params.Stages.Xpos-(wSpeckleShift.params.Stages.Xpos(1)-depthX)));

figure();
imagesc('XData', noSpeckleShift.params.Stages.Ypos, 'YData', noSpeckleShift.params.Stages.Xpos, ...
        'CData', permute(noSpeckleShift.I(:,:,zCoor), [2,1,3]));
set(gca, 'FontSize', axesFontSize);
title("No Speckle Shift", 'FontSize', titleFontSize)
caxis(limits)
axis equal
axis tight
colorbar;

zCoor = 11;

figure()
imagesc('XData', wSpeckleShift.params.Stages.Ypos, 'YData', wSpeckleShift.params.Stages.Xpos(1:idx), ...
        'CData', permute(wSpeckleShift.I(:,1:idx,zCoor), [2,1,3]));
set(gca, 'FontSize', axesFontSize);
title("With Speckle Shift")
caxis(limits)
colorbar;
axis equal
axis tight

figure()
imagesc('XData', wSpeckleShift.params.Stages.Ypos, 'YData', wSpeckleShift.params.Stages.Xpos(1:idx), ...
        'CData', permute(wSpeckleShift.I(:,1:idx,zCoor), [2,1,3]));
set(gca, 'FontSize', axesFontSize);
title("With Speckle Shift - Scaled Colors")
% caxis(limits)
colorbar;
axis equal
axis tight

figure()
imagesc('XData', wSpeckleShiftNarrow.params.Stages.Ypos, 'YData', wSpeckleShiftNarrow.params.Stages.Xpos, ...
        'CData', permute(wSpeckleShiftNarrow.I(:,:,zCoor), [2,1,3]));
set(gca, 'FontSize', axesFontSize);
title("With Speckle Shift - Fine Scan")
colorbar;
axis equal;
axis tight;
