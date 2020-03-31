close all;
clear all;
clc;

%----------------------------------------------------------%
% In this file: changing the optical intensity and check SNR behavior
% 1. init Acousto Optics element
% 2. Comparing with ultrasound naive and multiplexed analysis
% 3. Comparing without ultrasound naive and multiplexed analysis
%----------------------------------------------------------%

ao = acoustoOptics();
ao.init();

uVars = acoustoOptics.uVarsCreate();

uVars.c                 = 1445;
uVars.fSin              = 1.25e6;              
uVars.fTrain            = 16.5e3;
uVars.cycInPulse        = 4; 
uVars.channels          = 4; %update in digitizer
uVars.phantomDepth      = 4.6e-2;
uVars.distFromPhantom   = 7.15e-2;
uVars.fExtClk           = 5e6; %fs
uVars.fSclk             = 100e6;     %update in fGen
uVars.timeToSample      = 2;
uVars.extClkDcyc        = 50; % [%]
uVars.IOPort            = 1;
uVars.IOLine            = 4;
uVars.useGPU            = true; %algo, digitizer;
uVars.fastAnalysis      = false;
uVars.useQuant          = true;
uVars.quantTime         = 0.002;
uVars.useHadamard       = true;

uVars.gReq.ch  = 4;
uVars.gReq.pos = 1;
uVars.gReq.intExt = 'int';
names = fieldnames(uVars.gReq.validStruct);
for i=1:length(names)
    uVars.gReq.validStruct.(names{i}) = false;
end

uVars.gReq.validStruct.usSignal = true;
uVars.gReq.validStruct.phi = true;


uVars.exportRawData.rawData       = false;
uVars.exportRawData.netSignal     = false;
uVars.exportRawData.deMultiplexed = false;
uVars.exportRawData.reshape       = false;
uVars.exportRawData.FFT           = false;

uVars.fs.scanName        = [];
uVars.fs.resDirPath      = 'D:\Results';
uVars.fs.saveFullData    = false;
uVars.fs.saveFigs        = false;
uVars.fs.saveResults     = false;

repPerConf = 30;

powerVec = [25,50,100,150,200];

%%

% for i =1:length(powerVec)
%     beep
%     f = msgbox(sprintf('Set Optical power to %d',powerVec(i)) ,'Done');
%     waitfor(f)
%     fprintf("Runnins no ultrasound measurement");
% 
%     %no Hadamard:
%     uVars.useHadamard = false;
%     uVars.useGPU = true;
%     name = sprintf("Comparison-woHad-withUltrasound-4Cyc-17n");
%     uVars.fs.scanName = name;
% 
%     ao.setMeasVars(uVars);
%     numOfPosNaive = ao.algo.samples.numOfPos;
%     ao.configPeripherals();
% 
%     fprintf(name);
%     for k = 1:repPerConf
%         resNaive{i,k} = ao.measureAndAnlayse();
%         phiNaive(i,k,:) = resNaive{i,k}.phi;
%         disp(k);
%         pause(0.01);
%     end
% 
%     %with Hadamard:
%     uVars.useHadamard = true;
%     uVars.useGPU = false;
%     name = sprintf("Comparison-wHad-withUltrasound-4Cyc-17n");
%     uVars.fs.scanName = name;
%     fprintf(name);
% 
%     ao.setMeasVars(uVars);
%     numOfPosHad = ao.algo.samples.numOfPos;
%     ao.configPeripherals();
% 
%     fprintf(name);
%     for k = 1:repPerConf
%         resHad{i,k} = ao.measureAndAnlayse();
%         phiHad(i,k,:) = resHad{i,k}.phi;
%         disp(k);
%         pause(0.01);
%     end
% end
% 
% save("./matFiles/FellgettsAdv3.mat", '-v7.3')

%% Analysis
load("./matFiles/FellgettsAdv3.mat");

idx = 8;
env = 3;
vars = ao.getMeasVars();
zVec = vars.algo.len.zVecUSRes;

meanSigHad   = squeeze(mean(phiHad, 2));
meanSigNaive = squeeze(mean(phiNaive,2));

stdSigHad = squeeze(std(phiHad, 0, 2));
stdSigNaive = squeeze(std(phiNaive, 0, 2));

figure();
for i = 1:length(powerVec)
    subplot(2,3,i)
    errorbar(zVec, meanSigHad(i,:), stdSigHad(i,:)); hold on
    errorbar(zVec, meanSigNaive(i,:), stdSigNaive(i,:));
    title(sprintf("%d[W]", powerVec(i)))
    legend("Hadamard", "Naive")
end

meanSigHadPeak   = meanSigHad(:,idx);
meanSigNaivePeak = meanSigNaive(:,idx);
meanSigHadTail   = mean(meanSigHad(:,idx+env:end), 2);
meanSigNaiveTail = mean(meanSigNaive(:,idx+env:end), 2);

stdSigHadPeak   = stdSigHad(:,idx);
stdSigNaivePeak = stdSigNaive(:,idx);

hadTail = phiHad(:,:,idx+env:end);
hadTail = reshape(hadTail, 5, []);
stdSigHadTail   = std(hadTail, 0 ,2);

naiveTail = phiNaive(:,:,idx+env:end);
naiveTail = reshape(naiveTail, 5, []);
stdSigNaiveTail = std(naiveTail, 0 ,2);

%% Peak Results
% Only Exp Ration at Peak
figure();
subplot(1,2,1)
plot(powerVec, meanSigHadPeak, '-+'); hold on
plot(powerVec, meanSigNaivePeak, '-+')
legend("Hadamard", "Naive")
xlabel("Laser Power[W]")
ylabel("E(sig) at Peak")
subplot(1,2,2)
plot(powerVec, meanSigHadPeak./meanSigNaivePeak, '-+')
xlabel("Laser Power[W]")
ylabel("E_{had}/E_{naive}")

% Only Std Ration at Peak
figure();
subplot(1,2,1)
plot(powerVec, stdSigHadPeak, '-+'); hold on
plot(powerVec, stdSigNaivePeak, '-+')
legend("Hadamard", "Naive")
xlabel("Laser Power[W]")
ylabel("STD at Peak Power")

subplot(1,2,2)
plot(powerVec,stdSigNaivePeak./stdSigHadPeak, '-+')
xlabel("Laser Power[W]")
ylabel("Std_{Naive} / Std_{had}")


%Complete Peak SNR & multiplexing advantage
SNRHadPeak   = meanSigHadPeak./stdSigHadPeak;
SNRNaivePeak = meanSigNaivePeak./stdSigNaivePeak;

figure();
subplot(1,2,1)
plot(powerVec, SNRHadPeak); hold on
plot(powerVec, SNRNaivePeak)
legend("Hadamard", "Naive")
xlabel("Laser Power[W]")
ylabel("SNR")

subplot(1,2,2)
plot(powerVec, SNRHadPeak./SNRNaivePeak)
xlabel("Laser Power[W]")
ylabel("Multiplexing Advantage at Peak Power")


%% Tail Results
% Only Exp Ration at Peak
figure();
subplot(1,2,1)
plot(powerVec, meanSigHadTail, '-+'); hold on
plot(powerVec, meanSigNaiveTail, '-+')
legend("Hadamard", "Naive")
xlabel("Laser Power[W]")
ylabel("E(sig) at Tail")
subplot(1,2,2)
plot(powerVec, meanSigHadTail./meanSigNaiveTail, '-+')
xlabel("Laser Power[W]")
ylabel("E_{had}/E_{naive}")

% Only Std Ration at Tail
figure();
subplot(1,2,1)
plot(powerVec, stdSigHadTail, '-+'); hold on
plot(powerVec, stdSigNaiveTail, '-+')
legend("Hadamard", "Naive")
xlabel("Laser Power[W]")
ylabel("STD at Tail")

subplot(1,2,2)
plot(powerVec,stdSigNaiveTail./stdSigHadTail, '-+')
xlabel("Laser Power[W]")
ylabel("Std_{Naive} / Std_{had}")


%Complete Tail SNR & multiplexing advantage
SNRHadTail   = meanSigHadTail./stdSigHadTail;
SNRNaiveTail = meanSigNaiveTail./stdSigNaiveTail;

figure();
subplot(1,2,1)
plot(powerVec, SNRHadTail); hold on
plot(powerVec, SNRNaiveTail)
legend("Hadamard", "Naive")
xlabel("Laser Power[W]")
ylabel("SNR")

subplot(1,2,2)
plot(powerVec, SNRHadTail./SNRNaiveTail)
xlabel("Laser Power[W]")
ylabel("Multiplexing Advantage at Peak Power")

%%

% SNR change, new definition (dynamic range / std)

snrHad   = abs(meanSigHadPeak   - meanSigHadTail)   ./ stdSigHadTail;
snrNaive = abs(meanSigNaivePeak - meanSigNaiveTail) ./ stdSigNaiveTail;

figure();
subplot(1,2,1)
plot(powerVec, snrHad, '-+'); hold on
plot(powerVec, snrNaive, '-x');
legend("Hadamard", "Naive")
xlabel("Laser Power[W]")
ylabel("SNR")
subplot(1,2,2)
plot(powerVec, snrHad./snrNaive, '-+');
xlabel("Laser Power[W]")
ylabel("Fellget Adv")


% SBR change
sbrHad   = abs(meanSigHadPeak   - meanSigHadTail)   ./ meanSigHadTail;
sbrNaive = abs(meanSigNaivePeak - meanSigNaiveTail) ./ meanSigNaiveTail;

figure();
subplot(1,2,1)
plot(powerVec, sbrHad, '-+'); hold on
plot(powerVec, sbrNaive, '-x');
legend("Hadamard", "Naive")
xlabel("Laser Power[W]")
ylabel("SBR")
subplot(1,2,2)
plot(powerVec, sbrHad./sbrNaive, '-+');
xlabel("Laser Power[W]")
ylabel("Fellget Adv")


