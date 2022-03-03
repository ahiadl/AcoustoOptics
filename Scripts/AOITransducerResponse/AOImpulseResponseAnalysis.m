% close all
% clear all
% clc

%%
% resPulse = load("14-Feb-2022 12-02-03-atFocus-73.6-pulse.mat");
% resPulse = load("14-Feb-2022 12-02-58-atFocus-95.6-pulse.mat");
resPulse = load("15-Feb-2022 12-09-38-atSinFocus-83.6-Pulse.mat");
varsPulse = resPulse.csVars;
resPulse = squeeze(resPulse.resCs);
N = length(resPulse);
fs = varsPulse.fs;

acCoupling = resPulse - mean(resPulse);
fftPulse = fftshift(fft(fftshift(resPulse)));
fftPulseNorm = abs(fftPulse).^2;

minVal = min(fftPulseNorm);
maxValPulse = max(fftPulseNorm);
span = maxValPulse - minVal;

fftPulseNorm = (fftPulseNorm - minVal)/span;

fftPulseAC = fftshift(fft(fftshift(acCoupling)));
phaseAC = angle(fftPulseAC) / pi;
fftPulseACNorm = abs(fftPulseAC).^2;

minVal = min(fftPulseACNorm);
maxValPulse = max(fftPulseACNorm);
span = maxValPulse - minVal;

fftPulseACNorm = (fftPulseACNorm - minVal)/span;

tVec = (0:1:N-1)/fs;
fVec = (fs/N) *  ( (-N/2) : 1 : (N/2)-1 );

% figure();
% subplot(1,2,1)
% plot(tVec*1e6, resPulse);
% subplot(1,2,2)
% plot(fVec/1e6, fftPulseNorm);
% xlim([-5,5])

figure();
subplot(1,2,1)
plot(tVec*1e6, acCoupling);
subplot(1,2,2)
plot(fVec/1e6, fftPulseACNorm); hold on
% plot(fVec/1e6, phaseAC)
xlim([-5,5])
%%
sClk = 100e6;
N = 65536;
sigTime = N/sClk;
% delay = 51.28e-6;
delay = 0;
% AO Signal
fsin = 1.25e6;

% Create a sin
tVecSin  = (0:1:((sClk/fsin)-1))*(1/sClk);
Nsin = length(tVecSin);

tVecSinPad = (0:1:N-1)./sClk;
fVecSig = (sClk/N) *  ( (-N/2) : 1 : (N/2)-1 ); 
sig = sin(2*pi*fsin*tVecSin);

delaySamples = floor(delay * sClk);
padLen = N-length(sig)-delaySamples;
sig = [zeros(1,delaySamples), sig, zeros(1,padLen)];

fVecSig = (sClk/N) *  ( (-N/2) : 1 : (N/2)-1 );

sigFFT = fftshift(fft(fftshift(sig)));
sigFFTNorm = abs(sigFFT).^2;

minVal = min(sigFFTNorm);
maxValSig = max(sigFFTNorm);
span = maxValSig - minVal;

sigFFTNorm = (sigFFTNorm - minVal)./span;

figure();
subplot(1,2,1)
plot(tVecSinPad*1e6, sig);
subplot(1,2,2)
plot(fVecSig/1e6, sigFFTNorm); hold on
plot(fVec/1e6, fftPulseACNorm);
xlim([-5,5])

%% Creat Mask
zeroIdx = find(fVecSig == 0);
sigma     =  0.25e6;
mu        = 1.059e6;
power = 4;
fVecGauss = flip(abs(fVecSig(1 : zeroIdx-1)));
gaussian =  (1/(sigma * sqrt(2*pi))) * exp( -0.5 * ( ((fVecGauss - mu).^power)./ (sigma^power) ) );

mask = [flip(gaussian), gaussian(1), gaussian(1:end-1)];

minVal = min(mask);
maxVal = max(mask);
span = maxVal - minVal;

maskNorm = ((mask - minVal)./span );

figure();
plot(fVecSig/1e6, maskNorm); hold on
plot(fVecSig/1e6, fftPulseACNorm)
plot(fVecSig/1e6, sigFFTNorm)
xlim([-5,5])

% Calcualate deformed signal
spanIdx = 1:N;
% spanIdx = (abs(fVecSig)<2.5e6 & abs(fVecSig)>0.1e6 );
% fVecConv = fVecSig(spanIdx);

convSigFFTNP = sigFFT(spanIdx)./fftPulse(spanIdx)'.*mask;
convSigFFTNP(convSigFFTNP == inf) = 0;

convSigFFT          = zeros(1, N);
convSigFFT(spanIdx) = convSigFFTNP;

convSigFFTNorm = abs(convSigFFT).^2;

convSig = (ifft(ifftshift(convSigFFT)));

Nconv = length(convSigFFT);
tVecConv = (0:1:Nconv-1)./sClk; 

minVal = min(convSigFFTNorm);
maxVal = max(convSigFFTNorm);
span = maxVal - minVal;

convSigFFTNorm = (convSigFFTNorm - minVal)./span;

figure();
plot(fVecSig/1e6, sigFFTNorm); hold on
plot(fVec/1e6,    fftPulseNorm);
plot(fVec/1e6,    convSigFFTNorm);
xlim([-5,5])


figure();
plot(convSig)

convSig = convSig - mean(convSig);

minVal = min(convSig);
maxVal = max(convSig);
span = maxVal - minVal;

convSigNorm = ((convSig - minVal)./span )*2 -1;
% convSigNorm = convSigNorm - mean(convSigNorm);

figure();
subplot(1,2,1)
plot(tVecConv*1e6, convSigNorm);
xlim([30,70])
subplot(1,2,2)
plot(fVec/1e6, convSigFFTNorm);
xlim([-5,5])

figure();
subplot(1,2,1)
plot(tVecConv*1e6, convSigNorm);
% xlim([30,70])
subplot(1,2,2)
plot(fVec/1e6, convSigFFTNorm);
xlim([-5,5])


figure();
subplot(1,2,1)
plot(convSigNorm);
% xlim([30,70])
subplot(1,2,2)
plot(fVec/1e6, convSigFFTNorm);
xlim([-5,5])


convPulse   = convSigNorm(5870:6468);
padLen      = N-length(convPulse);
convSigFull = -[convPulse,  mean(convSigNorm)*ones(1,padLen)];

figure(); plot(convSigFull);
