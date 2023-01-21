close all
clear all
clc;

%% Load Data
dir = "C:/Users\sahiadl.EED/OneDrive - Technion/Graduate/AcoustoOpticSystem/Measurements/Transducer Pressure Field/Calibrated/Focused";
resDepth = load(sprintf("%s/%s", dir, "1.25MHz-Depth-Axis"));
resTrans = load(sprintf("%s/%s", dir, "1.25MHz-Transversal-Axis"));
resPulse = load(sprintf("%s/%s", dir, "1.25MHz-FocalPoint"));
resHad   = load(sprintf("%s/%s", dir, "1.25MHz-FocalPoint-HadamardSqnc"));

% %% Load Data
% dir = "C:/Users\sahiadl.EED/OneDrive - Technion/Graduate/AcoustoOpticSystem/Measurements/Transducer Pressure Field/Calibrated/UnFocused";
% resDepth = load(sprintf("%s/%s", dir, "1.25MHz-Depth-Axis"));
% resTrans = load(sprintf("%s/%s", dir, "1.25MHz-Transversal-Axis"));
% resPulse = load(sprintf("%s/%s", dir, "1.25MHz-FocalPoint"));
% resHad   = load(sprintf("%s/%s", dir, "1.25MHz-FocalPoint-HadamardSqnc"));

%% Create Variables for 3D Scan:
csVars = resTrans.csVars;
csVars.axDisc2Span   = resDepth.csVars.axScanSpan;
csVars.axDisc2Span   = resDepth.csVars.axScanSpan;
csVars.axDisc2Stride = resDepth.csVars.axScanStride;
csVars.axDisc2Ref    = resDepth.csVars.axDisc2Ref;
csVars.disc2Vec      = resDepth.csVars.scanVec;
csVars.disc2VecBin   = resDepth.csVars.scanVecBin;

close all
%% Create Polar Grid:
close all;
transMeas = squeeze(resTrans.resCs);
depthMeas = squeeze(resDepth.resCs);

[tr1Size, tr2Size, depthSize] = size(transMeas);
tr1FocalIdx = floor(tr1Size/2) +1;
tr2FocalIdx = floor(tr2Size/2) +1;

depthMeasAC = depthMeas - mean(depthMeas, 2);
depthMeasPressure = (depthMeasAC*1e3/837) *1e3;

p2pDepth = peak2peak(depthMeasPressure, 2);
p2pLpf   = conv(p2pDepth, (1/7)*ones(1,7), 'same');
env      = envelope(p2pLpf, 10, 'peak');

p2pTrans = peak2peak(transMeas,3);
p2pTransNorm = p2pTrans/max(p2pTrans(:));

trans1Profile = p2pTransNorm(:, tr2FocalIdx);
trans2Profile = p2pTransNorm(tr1FocalIdx, :);

figure(); 
subplot(2,2,1)
plot(csVars.disc2Vec, p2pDepth); hold on;
plot(csVars.disc2Vec, p2pLpf);
plot(csVars.disc2Vec, env)
subplot(2,2,2)
imagesc(log(p2pTransNorm))
axis equal tight
colorbar
subplot(2,2,3)
plot(trans1Profile)
subplot(2,2,4)
plot(trans2Profile)

[maxPressure, focalIdx] = max(env);

%% Speed of sound
close all
clc 

dX = csVars.axDisc2Stride;
tVec = csVars.tVec;
for i=1:20
    [~, idxT(i)] = max(abs(depthMeasAC(focalIdx -(i-1),:)));
end
dTpeak =  abs(mean(tVec(idxT(2:end)) - tVec(idxT(1:end-1))));
c = dX*1e-3/dTpeak;

depthVec = c * tVec *1e3; %[mm]
dDepth = c * (tVec(2) - tVec(1)) *1e3;

sig = depthMeasAC(focalIdx, :);
[~, focalMaxDelayIdx] = max(sig);

gSig = gradient(sig);
envSize = 40;
idxVec = focalMaxDelayIdx-envSize:focalMaxDelayIdx+envSize;
gSigFoci = gSig(idxVec);
gSigFociNorm = gSigFoci - gSigFoci(1);
i = find(abs(gSigFociNorm)>(0.01*max(gSigFociNorm)) ,1);
focalLen = depthVec(idxVec(i));

scanVecFix = csVars.disc2VecBin;


figure(); 
subplot(1,2,1)
plot(sig); hold on
yyaxis right
% plot(gradient(sig)); hold on
plot(idxVec, gSigFociNorm);



%% Collect Vars
res.focalProfile = p2pDepth;
res.focalSig     = 

%%
dir = "C:/Users/sahiadl.EED/OneDrive - Technion/Graduate/AcoustoOpticSystem/Measurements/Transducer Pressure Field/Calibrated/Focused";
filename = "FocusedTransducerAnalysed.mat";
path = sprintf("%s/%s", dir, filename);
save(path,);
