close all
clear all
clc

load('C:\Users\sahiadl.EED\OneDrive - Technion\Graduate\AcoustoOpticSystem\Measurements\Transducer Pressure Field\AO_Transducer.mat');
%%
lenY = csVars.scanSizeBin(1);
lenZ = csVars.scanSizeBin(2);
lenX = csVars.scanSizeBin(3);

dY = csVars.axScanStride; 
dZ = csVars.axDisc1Stride; 
dX = csVars.axDisc2Stride;

yVec = csVars.scanVecBin - csVars.axScanRef;
zVec = csVars.disc1Vec - csVars.axDisc1Ref;
xVecRaw = flip(csVars.disc2Vec);

yDim = 1; zDim = 2; xDim = 3; chDim =4; tDim = 5;

triggerInd = 20; 
sensitivity = 837;

p2pField = peak2peak(resCs(:,:,:,:,triggerInd:end), tDim)/sensitivity;
[maxVal, indMax] = max(p2pField(:));

[maxY, maxZ, maxX] = ind2sub(size(p2pField), indMax);
peakSig = squeeze(resCs(maxY,maxZ,maxX,1,:));
peakSigAC = peakSig - mean(peakSig(60:500));


%% Calculate Speed Of Sound
offPeakSig = squeeze(resCs(maxY,maxZ,maxX+1,1,:));

[~, i] = max(peakSig);
[~, j] = max(offPeakSig);

tVec = csVars.tVec;
deltaT = abs(tVec(i) - tVec(j));
c = (dX/deltaT)/1000;

%% Calculate Frequency
N = length(tVec);
n = 0:N-1;
fs   = genVars.fs;
dfs  = fs/N;
fBar = fs * n /N;
fBarShift2 = ( ((-fs/2)+dfs) : dfs : (fs/2));
sigFFT = fftshift(fft(peakSig));
sigFFTAbs = abs(sigFFT);
sigFFTNorm = (sigFFTAbs - min(sigFFTAbs))/abs(max(sigFFTAbs) - min(sigFFTAbs));

T = tVec(664) - tVec(655);
fSig = 1/T;
fSig = 1.25e6;

%% Calculate Transducer distance from scan span
peakDistance = tVec(i)*c*1000;
spatialOffset =  peakDistance - xVecRaw(maxX); 
xVec = xVecRaw + spatialOffset;

%% Calculate Focuse distance (f), waist size (w0), and Rayleigh range (zR)
maxLine = squeeze(mean(p2pField(maxY-2:maxY+2,maxZ,:),1));
xSectionMax = squeeze(p2pField(:,maxZ, maxX));
% xSectionNorm = xSectionMax - min(xSectionMax);
K = 2*pi*fSig/c;
lambda = c/fSig*1e3; %[mm]
FWHM = dY*sum(xSectionMax>(maxVal*0.5));

f   = xVec(maxX)/10; %[cm]
w0  = FWHM/sqrt(2*log(2)); %[mm] %before interpolation]
Zr  = 0.5*K*((w0/1000)^2)*100; %[cm]

%% Interpolate
ext = 0;
intFactor = 10;
[Z, Y, X] = meshgrid(zVec, yVec, xVec);
yVecHR = yVec(1)-ext : dY/intFactor : yVec(end) + ext;
zVecHR = zVec(1)-ext : dZ/intFactor : zVec(end) + ext;
xVecHR = xVec(1) : -dX/10 : xVec(end);

[ZHR, YHR, XHR] = meshgrid(zVecHR, yVecHR, xVecHR);

p2pFieldInterp = interp3(Z, Y, X, p2pField, YHR, ZHR, XHR);

[~, indMaxHR] = max(p2pFieldInterp(:));
[maxYHR, maxZHR, maxXHR] = ind2sub(size(p2pFieldInterp), indMaxHR);

%% Diaplay
close all;

%peak signal - Pa vs Time and FFT
figure();
subplot(1,2,1)
plot(fBarShift2*1e-6, sigFFTNorm)
xlabel('f [MHz]')
ylabel('Normalized spectral Intensity [AU]')
subplot(1,2,2)
plot(tVec * 1e6, peakSigAC/sensitivity);
xlabel('t [\mus]')
ylabel('MPa')

%y xSection (at max Z)
%z xSection (at max Y)
%x xSection (at focus)

%Before interpolation
figure(); 
subplot(3,1,1)
imagesc(xVec, yVec, squeeze(p2pField(:,maxZ,:)));
xlabel( 'X [mm]')
ylabel( 'Y [mm]')
axis tight equal
h = colorbar;
ylabel(h, 'MPa')
subplot(3,1,2)
imagesc(xVec, zVec, squeeze(p2pField(maxY,:,:)));
xlabel( 'X [mm]')
ylabel( 'Z [mm]')
axis tight equal
h = colorbar;
ylabel(h, 'MPa')
subplot(3,1,3)
imagesc(yVec, zVec, squeeze(p2pField(:,:,maxX)));
xlabel( 'y [mm]')
ylabel( 'Z [mm]')
axis tight equal
h = colorbar;
ylabel(h, 'MPa')

%After interpolation
%y xSection (at max Z)
%z xSection (at max Y)
%x xSection (at focus)
figure(); 
subplot(3,1,1)
imagesc(xVecHR, yVecHR, squeeze(p2pFieldInterp(:,maxZHR,:)));
xlabel( 'X [mm]')
ylabel( 'Y [mm]')
axis tight equal
h = colorbar;
ylabel(h, 'MPa')
subplot(3,1,2)
imagesc(xVecHR, zVecHR, squeeze(p2pFieldInterp(maxYHR,:,:)));
xlabel( 'X [mm]')
ylabel( 'Z [mm]')
axis tight equal
h = colorbar;
ylabel(h, 'MPa')
subplot(3,1,3)
imagesc(yVecHR, zVecHR, squeeze(p2pFieldInterp(:,:,maxXHR)));
xlabel( 'y [mm]')
ylabel( 'Z [mm]')
axis tight equal
h = colorbar;
ylabel(h, 'MPa')

XY = squeeze(p2pField(:,maxZ,:));
XZ = squeeze(p2pField(maxY,:,:));
YZ = squeeze(p2pField(:,:,maxX));

minXY = min(XY(:));
maxXY = max(XY(:));
spanXY = maxXY-minXY;
XYNorm = (XY - minXY)/spanXY;

minXZ = min(XZ(:));
maxXZ = max(XZ(:));
spanXZ = maxXZ-minXZ;
XZNorm = (XZ - minXZ)/spanXZ;

minYZ = min(YZ(:));
maxYZ = max(YZ(:));
spanYZ = maxYZ-minYZ;
YZNorm = (YZ - minYZ)/spanYZ;

figure(); 
subplot(3,1,1)
imagesc(xVec, yVec, db(XYNorm));
xlabel( 'X [mm]')
ylabel( 'Y [mm]')
axis tight equal
h = colorbar;
ylabel(h, 'dB')
subplot(3,1,2)
imagesc(xVec, zVec, db(XZNorm));
xlabel( 'X [mm]')
ylabel( 'Z [mm]')
axis tight equal
h = colorbar;
ylabel(h, 'dB')
subplot(3,1,3)
imagesc(yVec, zVec, db(YZNorm));
xlabel( 'y [mm]')
ylabel( 'Z [mm]')
axis tight equal
h = colorbar;
ylabel(h, 'dB')

% plot(xVec,db(XZNorm (20,:)));
% xlabel( 'X [mm]')
% ylabel(h, 'dB')

figure(); 
subplot(3,1,1)
plot(xVec, db(XYNorm(41,:))); hold on
plot(xVec, -6*ones(1,length(xVec)))
xlabel( 'X [mm]')
ylabel('dB')
title("Pressure profile - Axial axis Z")
subplot(3,1,2)
plot(zVec, db(YZNorm(41,:))); hold on
plot(zVec, -6*ones(1,length(zVec)))
xlabel( 'Z [mm]')
ylabel('dB')
xlim([-10, 10])
title("Pressure profile - Transversal axis Z")
subplot(3,1,3)
plot(yVec, db(YZNorm(:,20))); hold on
plot(yVec, -6*ones(1,length(yVec)))
xlabel( 'Y [mm]')
ylabel('dB')
xlim([-10, 10])
title("Pressure profile - Transversal axis Y")