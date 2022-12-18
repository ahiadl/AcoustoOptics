%% 04-02-2022: 3D Scan 5cm depth: Simulations & Measurements
close all
clear all
clc

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

numOfPhantoms = length(resSim.muaVec);
pivIdx = floor(numOfPhantoms/2)+1;

for i = 1:numOfPhantoms
    Phantoms(i,1) = sprintf("%d",i);
end

% ----------------------------------------------
% Ground Truth
% ----------------------------------------------

% Extract Ground Truth:
muaGT   = resSim.muaVec'*10;
musPGT  = resSim.musVec'*10;
muEffGT = resSim.muEffVec'*10;
ratioGT = muEffGT./muEffGT(pivIdx);

Tprop = table(Phantoms, muaGT, musPGT, muEffGT)

% ----------------------------------------------
% Measurement
% ----------------------------------------------

% Extract measurements Variables and Results:
depthAx  = vars.grid.depthVec * 1e3;       % in mm
dx       = (depthAx(2) - depthAx(1)) / 10; % in cm
numOfAvg = size(res(1).phi,1) * size(res(1).phi,2);
phiLen   = size(res(1).phi,3);

idxHig = 102;
idxLow = 70;
clear phi phiRed phiNorm avgPhi
for i = 1:numOfPhantoms
    phi(i,:,:)      = reshape(res(i).phi, numOfAvg, phiLen);
    phiRed(i, :, :) = phi(i,:, idxLow:idxHig); 
    phiNorm(i,:,:)  = normMatf(phiRed(i,:,:), 3);
end

% Align All Measurements
[~,I]  = max(phiNorm,[],3);
M = max(I(:));
refIdx = mode(I(:));
minIdx = 3;
maxIdx = M+3;
idxVec = minIdx:1:maxIdx;
noi = length(idxVec);
shiftMat = I - refIdx;

depthAxAlign = depthAx(1:noi);

for i=1:5
    for j=1:9
        tmpIdx       = idxVec + shiftMat(i,j);
        phiAlign(i,j,:) = flip(phiNorm(i,j,tmpIdx));
    end
end

% Average Measurements
avgPhi = squeeze(mean(sqrt(phiAlign),2));

phiNetSq  = sqrt(phiAlign);
phiNetLog = log(phiNetSq);
phiNetLog(phiNetLog == -inf) = NaN;
phiAvgLog = squeeze(mean(phiNetLog,2,'omitnan'));

% Calculate MuEff and Ratio
idxHig = 12;
idxLow = 10;
idxPeak = 6;

dz = depthAx(2) - depthAx(1);
depthAxRed = depthAx(idxPeak:idxHig);

[gx, ~] = gradient(phiAvgLog(:, idxLow:idxHig));
muEffMeas = abs((mean(gx,2) / dz) *10);
ratioMeas = muEffMeas ./ muEffMeas(3);

% ----------------------------------------------
% Simulations
% ----------------------------------------------

% Extract Simulations:
depthVecSim = resSim.depthVec;
dzSim = depthVecSim(2) - depthVecSim(1);

for i=1:numOfPhantoms
    phiSim(i,:) = resSim.resArr{i}.phiHRMid2In;
end

% Normalize & Log:
phiSimNorm = normMatf(phiSim,2);
phiSimLog  = log(phiSimNorm);

% figure();
% plot(phiSimLog')

% Calculate MuEff and Ratio
idxHig = 50;
idxLow = 5;

[gxSim, ~] = gradient(phiSimLog(:, idxLow:idxHig));
muEffSim = abs( (mean(gxSim,2) / dzSim) *10 );
ratioSim = muEffSim / muEffSim(3);

% ----------------------------------------------
% Errors and Analysis
% ----------------------------------------------

% Calc Error
ErrMuEff = abs(muEffMeas-muEffSim)./muEffSim *100;
ErrRatio = abs(ratioMeas-ratioSim)./ratioSim * 100;

Tres = table(Phantoms, muEffGT, muEffSim, muEffMeas, ErrMuEff, ratioGT, ratioSim, ratioMeas, ErrRatio)

% Compare:
idxHig = 12;
idxLow = 6;
depthVecMeasRed = depthAx(idxLow:idxHig) - depthAx(idxLow); 
phiLogMeasRed = phiAvgLog(:, idxLow:idxHig);

idxHig = 16;
idxLow = 1;
depthVecSimRed = depthVecSim(idxLow:idxHig);
phiLogSimRed   =  phiSimLog(:, idxLow:idxHig);

% Phi Relations
idxT = 3;
for i=1:5
    phiRelMeas(i,:) =  phiLogMeasRed(i,idxT:end) ./ phiLogMeasRed(3,idxT:end);
    phiRelSim(i,:)  =  phiLogSimRed(i,idxT:end) ./ phiLogSimRed(3,idxT:end);
    
    fitModelMeas = fit(depthVecMeasRed(idxT:end)', phiRelMeas(i,:)',   'poly1');
    fitModelSim  = fit(depthVecSimRed(idxT:end)',  phiRelSim(i,:)', 'poly1');
    
    fitMeas(i,:) = fitModelMeas.p1*depthVecMeasRed+fitModelMeas.p2;
end

%% Load Transducer Pulse
resPulseObjX = load("../Measurements/Transducer Pressure Field/07-Apr-2022 10-59-12-FocalAxis-FullScan-1.25MHz.mat");
resPulseObjY = load("../Measurements/Transducer Pressure Field/07-Apr-2022 11-35-35-YAxis-FullScan-1.25MHz.mat");
resPulseObjZ = load("../Measurements/Transducer Pressure Field/07-Apr-2022 12-26-26-ZAxis-FullScan-1.25MHz.mat");

scanXVars = resPulseObjX.csVars;
scanXMat  = squeeze(resPulseObjX.resCs);

scanYVars = resPulseObjY.csVars;
scanYMat  = squeeze(resPulseObjY.resCs);

scanZVars = resPulseObjZ.csVars;
scanZMat  = squeeze(resPulseObjZ.resCs);

%%
tVec     = scanXVars.tVec;

%X:
scanVecX1 = scanXVars.scanVec;
p2pXScan = peak2peak(scanXMat,2);
p2pXNorm = normMatf(p2pXScan,1);
p2pXdB = db(p2pXNorm);
[maxValX, idxX] = max(p2pXScan);
sigXF = scanXMat(idxX,:);
sigXAC = sigXF - mean(sigXF);
sigXFNorm = sigXAC/max(abs(sigXAC));

[envUp, envDo] = envelope(sigXFNorm, 10, 'peaks');

% Extract Speed of Sound & Align X-Vec:
dx = scanVecX1(2) - scanVecX1(1);
sigXF2 = scanXMat(idxX+1,:);
for i=1:10
    [~, idxT(i)] = max(abs(scanXMat(idxX +(i-1),:)));
end
dTpeak =  mean(tVec(idxT(2:end)) - tVec(idxT(1:end-1)));
c = dx/dTpeak;
c = 1500e3;
shiftFactor = round(abs(tVec(idxT(1))*c - scanVecX1(idxX)),5);
scanVecX = scanVecX1 + shiftFactor;
xAxisMeas = tVec*c;

%Calculate rayleigh length
pos6dB      = scanVecX1(p2pXdB >= -6);
focalLen    = abs(pos6dB(1) -  pos6dB(end));
focalPoint  = scanVecX(idxX);

%Pulse Parameters:
envNorm    = envUp./max(envUp);
envdB      = db(envNorm);
pos6dB     = xAxisMeas(envdB >= -6);
pulseDepth = abs(pos6dB(1) -  pos6dB(end));

% Y:
scanVecY1 = scanYVars.scanVec;
p2pYScan = peak2peak(scanYMat,2);
p2pYNorm = normMatf(p2pYScan,1);
p2pYdB = db(p2pYNorm);
[maxValY, idxY] = max(p2pYScan);
sigYF = scanYMat(idxY,:);
scanVecY = scanVecY1 - scanVecY1(idxY);

% Z:
scanVecZ1 = scanZVars.scanVec;
p2pZScan = peak2peak(scanZMat,2);
p2pZNorm = normMatf(p2pZScan,1);
p2pZdB = db(p2pZNorm);
[maxValZ, idxZ] = max(p2pZScan);
sigZF = scanZMat(idxZ,:);
scanVecZ = scanVecZ1 - scanVecZ1(idxZ);

% Waist Size:
pos6dBY = scanVecY(p2pYdB >= -6);
pos6dBZ = scanVecZ(p2pZdB >= -6);
waistSize = mean ([abs(pos6dBY(1) -  pos6dBY(end)), abs(pos6dBZ(1) -  pos6dBZ(end))]);

% Max Pressure:
maxPressure = mean([maxValX, maxValY, maxValZ])/837*1e3; %[KPa]

Parameters = ["Focal Point"; "Focal Len"; "Waist Size"; "Pulse Depth"; "Max Pressure"];
Units = ["mm"; "mm"; "mm"; "mm"; "KPa"];
Values = [focalPoint; focalLen; waistSize; pulseDepth; maxPressure];
TransducerParams = table(Parameters, Values, Units)


% Temporal Response along focal axis
x = 60:10:140;
for i=1:length(x)
    [~, idxs(i)] = min(abs(scanVecX - x(i)));
end
sampleSig = scanXMat(idxs,:);

x2 = 54:5:94;
for i=1:length(x)
    [~, idxs2(i)] = min(abs(scanVecX - x2(i)));
end
sampleSig2 = scanXMat(idxs2,:);

%% Convolve Simulations:
dxt        = 0.005;
alignPoint = focalPoint+20;

xtVec   = xAxisMeas(1):dxt:xAxisMeas(end);
sigFInt = interp1(xAxisMeas, sigXFNorm, xtVec);

[envupInt, ~] = envelope(sigFInt, 10, 'peaks');

[~, pulsePivotIdx]    = max(sigFInt);
pulsePivotPos         = xtVec(pulsePivotIdx);
pulseNegDist          = 2;
pulsePosDist          = 15;
pulseNegPos           = pulsePivotPos - pulseNegDist;
pulsePosPos           = pulsePivotPos + pulsePosDist;
[~, pulseNegWidthIdx] = min(abs(xtVec - pulseNegPos));
[~, pulsePosWidthIdx] = min(abs(xtVec - pulsePosPos));
pulseFullWidthIdx     = pulsePosWidthIdx - pulseNegWidthIdx + 1;
pulseFullLen          = pulseFullWidthIdx*dxt;

figure();
subplot(1,2,1)
plot(xtVec, sigFInt);
xlim([pulseNegPos, pulsePosPos])
ylim([-1.1, 1.1])
xlabel('X [mm]');
ylabel('Pressure [MPa]');
subplot(1,2,2)
plot(sigFInt); hold on,
plot(envupInt);
xlim([pulseNegWidthIdx, pulsePosWidthIdx])
ylim([-1.1, 1.1])
xlabel('X [mm]');
ylabel('Pressure [MPa]');

%% -----------------------------------------------------------
% Interpulate Simulations to match temporal Resolution:
%------------------------------------------------------------
% Interpolate simulation for 0.01 mm res
simDepth =  abs(depthVecSim(end)- depthVecSim(1));
depthVecSimInt = depthVecSim(1) : dxt: depthVecSim(end);
phIdxVec = 1:1:numOfPhantoms;

[X, Y]    = meshgrid(depthVecSim, phIdxVec);
[Xq, Yq]  = meshgrid(depthVecSimInt, phIdxVec);
phiSimInt = interp2(X,  Y, phiSimNorm, Xq, Yq);

% Align to transducer direction
depthVecSimAlign = flip(abs(depthVecSimInt - alignPoint));
phiSimAlign      = flip(phiSimInt,2);

%% ------------------------------------------------------------
% Pad with zeros for convolution
%------------------------------------------------------------
%simulation
phiSimIntPad   = [ zeros(numOfPhantoms,pulseFullWidthIdx), phiSimAlign,  zeros(numOfPhantoms,pulseFullWidthIdx)]';
depthVecSimPad = [ dxt*(-pulseFullWidthIdx:1:-1)+depthVecSimAlign(1), depthVecSimAlign,  dxt*(1:pulseFullWidthIdx)+depthVecSimAlign(end)];

% if depthVecSimPad(1) < scanVecX(1)
%     
% end

%% -------------------------------------------------------------
% Interpolate field to match spatial resolution:
%-------------------------------------------------------------
% Interpolate1: interpolate to a dxt res on time axis:
[X,Y]       = meshgrid(xAxisMeas, scanVecX);
[Xq, Yq]    = meshgrid(xtVec,scanVecX);
scanXMatInt = interp2(X,Y, scanXMat, Xq, Yq);

% Cut Relevant Indices On Time Axis:
% [xtVec(1), xtVec(end)]
% [depthVecSimPad(1), depthVecSimPad(end)]
mask      = (xtVec >= depthVecSimPad(1)-dxt) & (xtVec <= (depthVecSimPad(end)+dxt));
xtRed     = xtVec(mask);
scanXMatRed1 = scanXMatInt(:, mask);

% Interpolate2: interpolate to a dxt res on scan axis:
[X,Y]        = meshgrid(xtRed, scanVecX);
scanVecInt   = scanVecX(1):dxt:scanVecX(end);
[Xq, Yq]     = meshgrid(xtRed, scanVecInt);
scanXMatInt2 = interp2(X,Y, scanXMatRed1, Xq, Yq, 'cubic');

% Cut Relevant Indices On Scan Axis:
% [scanVecInt(1), scanVecInt(end)]
% [depthVecSimPad(1), depthVecSimPad(end)]

mask        = (scanVecInt >= depthVecSimPad (1)) & (scanVecInt <= (depthVecSimPad (end)+dxt));
scanVecRed  = scanVecInt(mask);
scanMat     = scanXMatInt2(mask,:)';

% figure(); imagesc(scanVecRed, xtRed, scanMat); axis tight equal

scanMatNorm = scanMat/max(scanMat(:));
scanMatNorm = scanMatNorm - mean(scanMatNorm(:));
%% convolve: 
[envMat, ~]  = envelope(scanMatNorm', 100, 'peaks');
envMat = envMat';
envMatNorm = envMat/max(envMat(3000*size(envMat,1):((size(envMat,1)^2-size(envMat,1)*1000))));

figure(); 
plot(envMatNorm(1:1000:end,:)'); hold on
plot(phiSimIntPad(:,1))


phiSimConv     = envMatNorm*phiSimIntPad;
phiSimConvNorm = normMatf(phiSimConv,1);
phiSimConvLog  = log(phiSimConvNorm);
% phiSimConvLog = log(phiSimConv);
figure();
subplot(1,2,1)
plot(phiSimConv); hold on
subplot(1,2,2)
plot(phiSimConvLog)

depthVecMeas = -flip(depthAxAlign) + (alignPoint )+6;
%%

figure()
for i=1:numOfPhantoms
    ax(i) = subplot(3,2,i);
    plot(depthVecSimPad, log(phiSimIntPad(:,i))); hold on
    plot(scanVecRed+1, phiSimConvLog(:,i))
    plot(depthVecMeas, flip(phiAvgLog(i,:)))
%     xlim([65,110])
    legend("Clean Sim", "Conv Sim", "Meas.", 'Location', 'northwest')
end
linkaxes(ax);
