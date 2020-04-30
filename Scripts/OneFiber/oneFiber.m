close all
clear all
clc

res  = load('D:\Results\16-Mar-2020 16-10-25-2DScan-1Fiber\Results.mat');
vars = load('D:\Results\16-Mar-2020 16-10-25-2DScan-1Fiber\Vars.mat');
res2 = load('D:\Results\18-Dec-2019 09-04-45-Multiplexed-1CycRes-fulls2dScan\phi.mat');
vars2 = load('D:\Results\18-Dec-2019 09-04-45-Multiplexed-1CycRes-fulls2dScan\Vars.mat');
%%
yAxis = vars.curVars.stages.firstVec;
xAxis = vars.curVars.stages.secondVec;

yAxisNorm = yAxis - yAxis(ceil(length(yAxis)/2));
xAxisNorm = -(xAxis - xAxis(1));

yAxis2 = vars2.curVars.stages.firstVec;
xAxis2 = vars2.curVars.stages.secondVec;

yAxisNorm2  = yAxis2 - yAxis2(ceil(length(yAxis2)/2));
xAxisNorm2 = -(xAxis2 - xAxis2(1));

zIdx1 = 29;
zIdx2 = 29;

% Extract Original Images
imgOrig1 = res.res.phi(:,:,zIdx1);
imgOrig2 = res2.res2.phi(:,:,zIdx2);

% Change Scale: Normalize and db
imgNorm1 = (imgOrig1 - min(min(imgOrig1))) /  abs(max(max(imgOrig1)) - min(min(imgOrig1)) );
imgLog1  = db(imgNorm1); 
imgNorm2 = (imgOrig2 - min(min(imgOrig2))) /  abs(max(max(imgOrig2)) - min(min(imgOrig2)) );
imgLog2 = db(imgNorm2);


minLog1 = mink(imgLog1(:),2);
maxLog1 = maxk(imgLog1(:),2);
minLog2 = mink(imgLog2(:),2);
maxLog2 = maxk(imgLog2(:),2);

climsLog = [min([minLog2(2),minLog1(2)]), max([maxLog1(1),maxLog2(1)])];


textFontSize  = 30;
labelFontSize = 28;
axesFontSize  = 28;
codeFontSize  = 28;
%% Linear
hFig = figure();
ax = subplot(1,2,1);
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
imagesc(ax, 'XData', yAxisNorm, ...
            'YData', xAxisNorm,...
            'CData', imgNorm1);
axis equal
axis tight
colorbar;

ax2 = subplot(1, 2, 2);
imagesc(ax2, 'XData', yAxisNorm2, ...
             'YData', xAxisNorm2, ...
             'CData', imgNorm2); hold on
axis equal
axis tight

hCB1 = colorbar();
hCB1L = get(hCB1, 'label');
set(hCB1, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
set(ax2, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')

%% Logarithmic
hFig = figure();
ax = subplot(1,2,1);
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
imagesc(ax, 'XData', yAxisNorm, ...
            'YData', xAxisNorm,...
            'CData', imgLog1);
axis equal
axis tight
caxis(ax, climsLog)

ax2 = subplot(1, 2, 2);
imagesc(ax2, 'XData', yAxisNorm2, ...
             'YData', xAxisNorm2, ...
             'CData', imgLog2); hold on
axis equal
axis tight

hCB1 = colorbar();
hCB1L = get(hCB1, 'label');
set(hCB1, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
set(ax2, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
caxis(ax2, climsLog)

%% Mu effective analysis

figure()
ax = subplot(1,2,1);
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
imagesc(ax, 'XData', yAxisNorm, ...
            'YData', xAxisNorm,...
            'CData', imgNorm1);
axis equal
axis tight
xlabel("Y[mm]")
ylabel("X[mm]")
colorbar;

ax = subplot(1,2,2);
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
imagesc(ax, imgNorm1);
axis equal
axis tight
xlabel("Y[mm]")
ylabel("X[mm]")
colorbar;

OFimg = imgNorm1;

phi_sq = OFimg(:, 25);
phi = sqrt(phi_sq);
grad_phi = gradient(phi);
dx = xAxisNorm(2)-xAxisNorm(1);

figure()
ax=axes();
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
plot(ax, xAxisNorm, phi_sq)
% title("Single Mid-Line")
xlabel("X(Depth)(mm)")
ylabel('\phi^2')

figure()
ax=axes();
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
plot(ax, xAxisNorm, phi)
% title("Single Mid-Line")
xlabel("X(Depth)(mm)")
ylabel('\phi')

figure()
ax=axes();
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
plot(ax, xAxisNorm, -grad_phi/(dx))
% title("Single Mid-Line")
xlabel("X(Depth)(mm)")
ylabel("\nabla\phi")

figure()
ax=axes();
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')
plot(ax, xAxisNorm, grad_phi./(dx*phi))
% title("Single Mid-Line")
xlabel("X(Depth)(mm)")
ylabel('\d(phi)/(dx*phi)')

avg_muEff = mean(grad_phi./(dx*phi));
std_mueff = std(grad_phi./(dx*phi));

phi_sq_avg = OFimg(:, 20:30);
phi_avg = sqrt(phi_sq_avg);
[~, grad_phi_avg] = gradient(phi_avg);

mu_eff_avg = -grad_phi_avg./(dx*phi_avg);
mu_eff_tot = mean(mu_eff_avg,2);
mu_eff_tot_std = std(mu_eff_avg,0,2);

mu_eff_avg_val = mean(mu_eff_tot);
mu_eff_avg_std_val = std(mu_eff_tot);

figure()
ax=axes();

errorbar(ax, xAxisNorm, mu_eff_tot, mu_eff_tot_std)
title("Average \mu_{eff}")
xlabel("X(Depth)(mm)")
ylabel('\mu_{eff}')
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')

phi_sq_avg = mean(OFimg(:, 20:30),2);
phi_avg = sqrt(phi_sq_avg);
[grad_phi_avg] = gradient(phi_avg);

mu_eff_avg = -grad_phi_avg./(dx*phi_avg);
mu_eff_tot = mean(mu_eff_avg);
mu_eff_tot_std = std(mu_eff_avg);

figure()
ax=axes();

plot(ax, xAxisNorm, mu_eff_tot)
title("Average Mid-Lines")
xlabel("X(Depth)(mm)")
ylabel('\mu_{eff}')
set(ax, 'FontSize', axesFontSize, 'TickLabelInterpreter', 'latex')



