close all
clear all
clc

%%

% resPulseObj = load("../../Other/contScan/12-Feb-2022 19-43-44-atFocus-1.25MHz-sin.mat");
resPulseObj = load("../../Other/contScan/14-Feb-2022 11-58-37-atFocus-73.6-1.25MHz-sin.mat");
% resPulseObj = load("../../Other/contScan/14-Feb-2022 12-04-37-atFocus-95.6-1.25MHz.mat");
varsPulse = resPulseObj.csVars;
resPulse = squeeze(resPulseObj.resCs);

N = length(resPulse);
fs = varsPulse.fs;

% AC Coupling
acCoupling = resPulse - mean(resPulse);

% FFT
fftPulseAC     = fftshift(fft(fftshift(acCoupling)));
phaseAC        = angle(fftPulseAC) / pi;
fftPulseACNorm = abs(fftPulseAC).^2;

% FFT Normalization (AC)
minVal = min(fftPulseACNorm);
maxValPulse = max(fftPulseACNorm);
span = maxValPulse - minVal;

fftPulseACNorm = (fftPulseACNorm - minVal)/span;

tVec = (0:1:N-1)/fs;
fVec = (fs/N) *  ( (-N/2) : 1 : (N/2)-1 );

figure();
subplot(1,2,1)
plot(tVec*1e6, acCoupling);
subplot(1,2,2)
plot(fVec/1e6, fftPulseACNorm); hold on
xlim([-5,5])

figure();
subplot(1,2,1)
plot(acCoupling);
subplot(1,2,2)
plot(fVec/1e6, fftPulseACNorm); hold on
xlim([-5,5])

%% envelope
idxLow = 4778;
idxHig = 7529;
noi = idxHig-idxLow+1;
pulseSig = acCoupling(idxLow:idxHig);
depthSig = 1480*tVec(1:noi)*1e3;
[envUp, envDo] = envelope(pulseSig, 10, 'peaks');

figure();
plot(depthSig, envUp)
