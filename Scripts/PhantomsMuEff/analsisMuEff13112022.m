close all;
clear all;
clc;

%% Load Simulations:
resSim    = load("13-Nov-2022 14-22-49-FocusedTrans.mat");
phiSim    = resSim.simData.phiAO;
phiSimMid = resSim.simData.midConv;
varsSim   = resSim.simVars;
display(varsSim.table.T)

muEffSim  = varsSim.pp.simMuEff;
muEffCalc = resSim.simData.calcMuEff;

for i = 1:4
    ratSim(i) = muEffSim(i+1) / muEffSim(i);
end

% figure();
% plot(phiSimMid');

[~, idxSim] = max(phiSimMid,[],2);
offset = 113 - idxSim;
axVecSimAligned = varsSim.cnv.axVecConv - varsSim.cnv.axVecConv(113);

for i = 1:5
    alignSim(i,:) = circshift(phiSimMid(i,:), offset(i));
end

figure();
plot(axVecSimAligned, alignSim');

varsNames =  ["Sample No. [#]"; "Sim. MuEff [1/mm]"; "Calc. MuEff[1/mm]"; "Meas. MuEff[1/mm]"];
%% Load Measurements:
for i=1:5
    if i==1
        filename = "D:/MuEff/Uniform/Focused/Phantom1-Focused/2DResults/Z-109.00/AOResults/AO-R-1-Y-52.00-Results.mat";
    else
        filename = sprintf("D:/MuEff/Uniform/Focused/Phantom%d-Focused/AO-Results.mat",i);
    end
    resMeas(i) = load(filename);
    phiMeas(i,:) = resMeas(i).phi;
    phiNorm(i,:) = resMeas(i).phiNorm;
    phiMeasLog(i,:) = resMeas(i).phiLog;
    phiRaw(i,:) = resMeas(i).rawPhi;
end

varsMeas = load(sprintf("D:/MuEff/Uniform/Focused/Phantom%d-Focused/AO-Vars.mat",i));

%% Cut Phi:
lowIdx = 1;
highIdx = 103;
tailIdx = 43;

phiLogCut  = flip(phiMeasLog(:,lowIdx:highIdx),2);
phiCut     = flip(phiMeas(:,lowIdx:highIdx),2);
phiCutNorm = flip(phiNorm(:,lowIdx:highIdx),2);
dX = varsMeas.measVars.algo.len.dDepth*1e3;
xVecRaw = varsMeas.measVars.algo.len.depthZero(1:length(phiLogCut))*1e3;
xVecRawCut = xVecRaw(1:size(phiCut,2));

maxVals = max(phiCut, [], 2);
minVals = mean(phiCut(:,tailIdx:end),2);
span = maxVals-minVals;
phiCutNorm = (phiCut-minVals)./span;

figure()
subplot(1,3,1)
plot(xVecRawCut,phiCut')
xlabel("X[mm]")
ylabel("Fluence [AU]")
legend("1", "2", "3", "4", "5")
subplot(1,3,2)
plot(xVecRawCut,phiCutNorm')
xlabel("X[mm]")
ylabel("Fluence [AU]")
subplot(1,3,3)
plot(xVecRawCut,phiLogCut')
xlabel("X[mm]")
ylabel("Fluence [log]")

%% Interpolation & Alignment:
xVecInt   = (xVecRaw(1)):varsSim.us.dAx:(xVecRaw(end));
phiLogInt = normMatf(interp1(xVecRaw, phiLogCut', xVecInt, 'spline')', 2);
phiInt    = interp1(xVecRaw, phiCut', xVecInt, 'spline')';

[~, idx] = max(phiInt, [], 2);
offset = 103-idx;

if exist('phiIntAligned', 'var'); clear phiIntAligned; end
for i =1:5
    phiIntAligned(i,:) = circshift(phiInt(i,:), offset(i));
end

xVecIntAligned = xVecInt - xVecInt(103);

figure();
plot(xVecIntAligned, phiIntAligned)
xlabel("X[mm]")
ylabel("Fluence [AU]")
xlim([-5,40])

%% Type 1 Normalization:
% Normalization:
maxVals = max(phiIntAligned, [], 2);
minVals = min(phiIntAligned, [], 2);
span = maxVals-minVals;
phiIntNorm1Aligned = (phiIntAligned-minVals)./span;

% Log:
phiIntLog1 = log(phiIntNorm1Aligned);

% Calc Slope:
figure();
plot(phiIntLog1')

slopeStart = 120;
slopeEnd   = 159;

vec      = slopeStart:slopeEnd;
xVecFit  = xVecInt(1:length(vec));
gradVals1 = mean(abs(gradient(phiIntLog1(:, slopeStart:slopeEnd))),2) /2 /varsSim.us.dAx;

for i=1:5
    fitModel = fit(xVecFit', phiIntLog1(i,slopeStart:slopeEnd)', 'poly1');
    gradValsFit1(i) = abs(fitModel.p1)/2;
end


T1 = table((1:varsSim.numOfPhantoms)', varsSim.pp.simMuEff', resSim.simData.calcMuEff',  gradVals1, 'VariableNames', varsNames);
display(T1);

figure();
subplot(1,2,1)
plot(xVecIntAligned, phiIntNorm1Aligned)
xlabel("X[mm]")
ylabel("Fluence [AU]")
xlim([-5,40])
subplot(1,2,2)
plot(xVecIntAligned, phiIntLog1)
xlabel("X[mm]")
ylabel("Fluence [log]")
xlim([-5,40]);
ylim([-5,0]);

figure();
for i=1:5
    subplot(2,3,i)
    plot(xVecIntAligned, phiIntLog1(i,:)); hold on
    plot(axVecSimAligned, alignSim(i,:))
    xlabel("X[mm]")
    ylabel("Fluence [AU]")
    legend("Sim", "Meas")
    title(sprintf("Phantom: %d", i))
    xlim([-5,40]);
    ylim([-5,0]);
end

%% Type 2 Normalization:
close all
% Normalization.
tailIdx = 304;
maxVals = max(phiIntAligned, [], 2);
minVals = mean(phiIntAligned(:,tailIdx:end),2);
span = maxVals-minVals;
phiIntNorm2Aligned = (phiIntAligned-minVals)./span;

% Log:
phiIntLog2 = log(abs(phiIntNorm2Aligned));

% Calc Slope:
% figure();
% plot(phiIntLog2')

idxVec = 101:241;
clear fitIdxVec;
clear xVecFit;
clear gradVals2;
clear phiVecFit;

for i=1:5
    [~, idxStart(i)] = min(abs(phiIntLog2(i,idxVec)+1.7)) ;
    idxStart(i) = idxStart(i) + idxVec(1) -1;
    mask = phiIntLog2(i,idxVec)<-2.5;
    cunSumVec = cumsum(mask);
    idxEnd(i) = idxVec(cunSumVec == 1);
    
    fitIdxVec{i}   = idxStart(i):idxEnd(i);
    xVecFit{i}     = xVecIntAligned(fitIdxVec{i});
    phiVecFit{i}   = phiIntLog2(i, fitIdxVec{i});
    gradVals2(i) = mean(abs(gradient(phiVecFit{i})),2) /2 /varsSim.us.dAx;
end


T2 = table((1:varsSim.numOfPhantoms)', varsSim.pp.simMuEff', resSim.simData.calcMuEff',  gradVals2', 'VariableNames', varsNames);
display(T2);

figure();
subplot(1,2,1)
plot(xVecIntAligned, phiIntNorm2Aligned)
xlabel("X[mm]")
ylabel("Fluence [AU]")
xlim([-5,40])
subplot(1,2,2)
plot(xVecIntAligned, phiIntLog2)
xlabel("X[mm]")
ylabel("Fluence [log]")
xlim([-5,40]);
ylim([-6,0]);

figure();
for i=1:5
    subplot(2,3,i)
    plot(xVecIntAligned, phiIntLog2(i,:)); hold on
    plot(axVecSimAligned, alignSim(i,:))
    plot(xVecFit{i}, phiVecFit{i}); hold on
    xlabel("X[mm]")
    ylabel("Fluence [AU]")
    legend("Sim", "Meas")
    title(sprintf("Phantom: %d", i))
    xlim([-5,40]);
    ylim([-7,0]);
end

figure()
plot(1:5, gradVals2, '-+'); hold on
plot(1:5, muEffSim,    '-+');
plot(1:5, muEffCalc,   '-+');
xlim([0, 6])
xlabel("Phantom Idx")
ylabel("MuEff")
legend("Meas", "Sim", "Calc")
%% Repeat Process for rawPhi:
% Cut Phi:
lowIdx = 12;
highIdx = 114;
tailIdx = 43;

phiRawCut     = flip(phiRaw(:,lowIdx:highIdx),2);

figure();
subplot(1,2,1)
plot(phiRaw')
subplot(1,2,2)
plot(phiRawCut')

dX = varsMeas.measVars.algo.len.dDepth*1e3;
xVecRaw = varsMeas.measVars.algo.len.depthZero(1:length(phiRawCut))*1e3;
xVecRawCut = xVecRaw(1:size(phiRawCut,2));

maxVals = max(phiRawCut, [], 2);
minVals = mean(phiRawCut(:,tailIdx:end),2);
span = maxVals-minVals;
phiRawCutNorm = (phiRawCut-minVals)./span;

phiRawLogCut = log(abs(phiRawCutNorm));


figure()
subplot(1,3,1)
plot(xVecRawCut,phiRawCut')
xlabel("X[mm]")
ylabel("Fluence [AU]")
legend("1", "2", "3", "4", "5")
subplot(1,3,2)
plot(xVecRawCut,phiRawCutNorm')
xlabel("X[mm]")
ylabel("Fluence [AU]")
subplot(1,3,3)
plot(xVecRawCut,phiRawLogCut')
xlabel("X[mm]")
ylabel("Fluence [log]")

% Interpolation & Alignment:
xVecRawInt   = (xVecRaw(1)):varsSim.us.dAx:(xVecRaw(end));
% phiRawLogInt = normMatf(interp1(xVecRaw, phiRawLogCut', xVecRawInt, 'spline')', 2);
phiRawInt    = interp1(xVecRaw, phiRawCut', xVecRawInt, 'spline')';

[~, idx] = max(phiRawInt, [], 2);
offset = 103-idx;

if exist('phiIntAligned', 'var'); clear phiIntAligned; end
for i =1:5
    phiRawIntAligned(i,:) = circshift(phiRawInt(i,:), offset(i));
end

xVecRawIntAligned = xVecRawInt - xVecRawInt(103);

figure();
plot(xVecRawIntAligned, phiRawIntAligned)
xlabel("X[mm]")
ylabel("Fluence [AU]")
xlim([-5,40])

% Calc MuEff via Raw Exponential Fitting:
phiRawIntAlignedExp = phiRawIntAligned(:,116:end);
xVecRawIntAlignedExp = xVecRawIntAligned(116:end) - xVecRawIntAligned(116);
phiRawIntAlignedExp1 = phiRawIntAlignedExp(1,:);
figure()
plot(xVecRawIntAlignedExp, phiRawIntAlignedExp')

for i =1:5
    
end



% Type 2 Normalization:
close all
% Normalization.
tailIdx = 304;
maxVals = max(phiRawIntAligned, [], 2);
minVals = mean(phiRawIntAligned(:,tailIdx:end),2);
span = maxVals-minVals;
phiRawIntNorm2Aligned = (phiRawIntAligned-minVals)./span;

% Log:
phiRawIntLog2 = log(abs(phiRawIntNorm2Aligned));

% Calc Slope:
figure();
plot(phiRawIntLog2')

%%
idxVec = 103:241;

clear fitRawIdxVec;
clear xVecRawFit;
clear gradValsRaw;
clear phiRawVecFit;

for i=1:5
    [~, idxStart(i)] = min(abs(phiRawIntLog2(i,idxVec)+2)) ;
    idxStart(i) = idxStart(i) + idxVec(1) -1;
    mask = phiRawIntLog2(i,idxVec)<-2.8;
    cumSumVec = cumsum(mask);
    idxEnd(i) = idxVec(cunSumVec == 1);
    
    fitRawIdxVec{i}   = idxStart(i):idxEnd(i);
    xVecRawFit{i}     = xVecRawIntAligned(fitRawIdxVec{i});
    phiRawVecFit{i}   = phiRawIntLog2(i, fitRawIdxVec{i});
    gradValsRaw(i)    = mean(abs(gradient(phiRawVecFit{i})),2) /2 /varsSim.us.dAx;
end
%%

T3 = table((1:varsSim.numOfPhantoms)', varsSim.pp.simMuEff', resSim.simData.calcMuEff',  gradValsRaw', 'VariableNames', varsNames);
display(T3);

figure();
for i=1:5
    subplot(2,3,i)
    plot(xVecRawIntAligned, phiRawIntLog2(i,:)); hold on
    plot(axVecSimAligned, alignSim(i,:))
    plot(xVecRawFit{i}, phiRawVecFit{i}); hold on
    xlabel("X[mm]")
    ylabel("Fluence [AU]")
    legend("Meas", "Sim")
    title(sprintf("Phantom: %d", i))
    xlim([-5,40]);
    ylim([-7,0]);
end

figure()
plot(1:5, gradValsRaw, '-+'); hold on
plot(1:5, muEffSim,    '-+');
plot(1:5, muEffCalc,   '-+');
xlim([0, 6])
xlabel("Phantom Idx")
ylabel("MuEff")
legend("Meas", "Sim", "Calc")


%% Maximal Parallel Region
figure();
plot(phiRawIntLog2')

levelEnd = [3.354, 3.36, 2.92, 3.196, 2.8];
idxVec = 103:263;

clear fitRawIdxVec;
clear xVecRawFit;
clear gradValsRaw;
clear phiRawVecFit;

for i=1:5
    [~, idxStart(i)] = min(abs(phiRawIntLog2(i,idxVec)+2)) ;
    idxStart(i) = idxStart(i) + idxVec(1) -1;
    mask = phiRawIntLog2(i,idxVec)<-levelEnd(i);
    cumSumVec = cumsum(mask);
    idxEnd(i) = idxVec(cumSumVec == 1);
    
    fitRawIdxVec{i}   = idxStart(i):idxEnd(i);
    xVecRawFit{i}     = xVecRawIntAligned(fitRawIdxVec{i});
    phiRawVecFit{i}   = phiRawIntLog2(i, fitRawIdxVec{i});
    gradValsRaw(i)    = mean(abs(gradient(phiRawVecFit{i})),2) /2 /varsSim.us.dAx;
end

T4 = table((1:varsSim.numOfPhantoms)', varsSim.pp.simMuEff', resSim.simData.calcMuEff',  gradValsRaw', 'VariableNames', varsNames);
display(T4);

figure();
for i=1:5
    subplot(2,3,i)
    plot(xVecRawIntAligned, phiRawIntLog2(i,:)); hold on
    plot(axVecSimAligned, alignSim(i,:))
    plot(xVecRawFit{i}, phiRawVecFit{i}); hold on
    xlabel("X[mm]")
    ylabel("Fluence [AU]")
    legend("Meas", "Sim")
    title(sprintf("Phantom: %d", i))
    xlim([-5,40]);
    ylim([-7,0]);
end

figure()
plot(1:5, gradValsRaw, '-+'); hold on
plot(1:5, muEffSim,    '-+');
plot(1:5, muEffCalc,   '-+');
xlim([0, 6])
xlabel("Phantom Idx")
ylabel("MuEff")
legend("Meas", "Sim", "Calc")

%% Fluence Simulation Example:
phiSrc   = normMatf(res.simData.phiSrc{1});
phiDet   = normMatf(res.simData.phiDet{1});
phiLight = normMatf(res.simData.phiLight{1});

xVec = res.simData.srcVars.xVec;
yVec = res.simData.srcVars.yVec;
zVec = res.simData.srcVars.zVec;

figure()
subplot(3,3,1)
imagesc(yVec, xVec, log(squeeze(phiSrc(:,:,121))))
axis tight equal
xlabel("Y[mm]")
ylabel("X[mm]")
subplot(3,3,4)
imagesc(zVec, yVec, log(squeeze(phiSrc(121,:,:))))
axis tight equal
xlabel("Z[mm]")
ylabel("Y[mm]")
subplot(3,3,7)
hs=slice(log(phiSrc), [1,241],[1, 241],[1 241]);
set(hs,'linestyle','none');
axis tight equal
xlabel("Y[mm]")
ylabel("X[mm]")
zlabel("Z[mm]")

subplot(3,3,2)
imagesc(yVec, xVec, log(squeeze(phiDet(:,:,121))))
axis tight equal
xlabel("Y[mm]")
ylabel("X[mm]")
subplot(3,3,5)
imagesc(zVec, yVec, log(squeeze(phiDet(121,:,:))))
axis tight equal
xlabel("Z[mm]")
ylabel("Y[mm]")
subplot(3,3,8)
hs=slice(log(phiDet), [1,241],[1, 241],[1 241]);
set(hs,'linestyle','none');
axis tight equal
xlabel("Y[mm]")
ylabel("X[mm]")
zlabel("Z[mm]")

subplot(3,3,3)
imagesc(yVec, xVec, log(squeeze(phiLight(:,:,121))))
axis tight equal
xlabel("Y[mm]")
ylabel("X[mm]")
subplot(3,3,6)
imagesc(zVec, yVec, log(squeeze(phiLight(121,:,:))))
axis tight equal
xlabel("Z[mm]")
ylabel("Y[mm]")
subplot(3,3,9)
hs=slice(log(phiLight), [1,241],[1, 241],[1 241]);
set(hs,'linestyle','none');
axis tight equal
xlabel("Y[mm]")
ylabel("X[mm]")
zlabel("Z[mm]")
