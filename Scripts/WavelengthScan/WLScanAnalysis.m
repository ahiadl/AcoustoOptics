%% -------------------
% 19-07-21
% --------------------

% -> SNR vs Quant PMT 0.05Hz 20s

dir = 'D:\WavelengthScan\SNR-Vs.-Quant-PMT-20s-0.05Hz';
load(sprintf('%s/19-Jul-2021 15-00-23-PMT-20s-ExpVars.mat', dir))
expVars.quant = expVars.quant(1:9);

figure()
subplot(1,2,1)
errorbar(expVars.quant, expVars.SNRAvg, expVars.SNRstd)
xlabel('Quant Time[s]')
ylabel('Avg. SNR')
set(gca,'Xscale','log')
subplot(1,2,2)
stem(expVars.quant, expVars.SNRstd./expVars.SNRAvg)
xlabel('Quant Time[s]')
ylabel('(SNR Std)/(SNR Avg)')
set(gca,'Xscale','log')

SNRVec(1) = expVars.SNRAvg(5); 
STDVec(1) = expVars.SNRstd(5);
timeVec(1) = 20;

% -> SNR vs Quant PMT 0.5Hz 2s
dir = 'D:\WavelengthScan\SNR-Vs.-Quant-PMT-2s-0.5Hz';
load(sprintf('%s/19-Jul-2021 11-59-34-PMT-2s-ExpVars.mat', dir))

figure()
subplot(1,2,1)
errorbar(expVars.quant, expVars.SNRAvg, expVars.SNRstd)
xlabel('Quant Time[s]')
ylabel('Avg. SNR')
set(gca,'Xscale','log')
subplot(1,2,2)
stem(expVars.quant, expVars.SNRstd./expVars.SNRAvg)
xlabel('Quant Time[s]')
ylabel('(SNR Std)/(SNR Avg)')
set(gca,'Xscale','log')

SNRVec(2) = expVars.SNRAvg(5); 
STDVec(2) = expVars.SNRstd(5);
timeVec(2) = 2;

% -> SNR vs Quant PMT 5Hz 2s
dir = 'D:\WavelengthScan\SNR-Vs.-Quant-PMT-2s-5Hz';
load(sprintf('%s/19-Jul-2021 16-32-25-PMT-5Hz-ExpVars.mat', dir))

figure()
subplot(1,2,1)
errorbar(expVars.quant, expVars.SNRAvg, expVars.SNRstd)
xlabel('Quant Time[s]')
ylabel('Avg. SNR')
set(gca,'Xscale','log')
subplot(1,2,2)
stem(expVars.quant, expVars.SNRstd./expVars.SNRAvg)
xlabel('Quant Time[s]')
ylabel('(SNR Std)/(SNR Avg)')
set(gca,'Xscale','log')

SNRVec(3) = expVars.SNRAvg(5); 
STDVec(3) = expVars.SNRstd(5);
timeVec(3) = 2;

% -> SNR vs Quant PMT 10Hz 2s
dir = 'D:\WavelengthScan\SNR-Vs.-Quant-PMT-2s-10Hz';
load(sprintf('%s/19-Jul-2021 19-20-54-PMT-2s-10Hz-ExpVars.mat', dir))

figure()
subplot(1,2,1)
errorbar(expVars.quant, expVars.SNRAvg, expVars.SNRstd)
xlabel('Quant Time[s]')
ylabel('Avg. SNR')
set(gca,'Xscale','log')
subplot(1,2,2)
stem(expVars.quant, expVars.SNRstd./expVars.SNRAvg)
xlabel('Quant Time[s]')
ylabel('(SNR Std)/(SNR Avg)')
set(gca,'Xscale','log')

SNRVec(4) = expVars.SNRAvg(5); 
STDVec(4) = expVars.SNRstd(5);
timeVec(4) = 2;

% -> SNR vs Quant PMT 25Hz 2s
dir = 'D:\WavelengthScan\SNR-Vs.-Quant-PMT-2s-25Hz';
load(sprintf('%s/19-Jul-2021 19-48-40-PMT-2s-25Hz-ExpVars.mat', dir))

figure()
subplot(1,2,1)
errorbar(expVars.quant, expVars.SNRAvg, expVars.SNRstd)
xlabel('Quant Time[s]')
ylabel('Avg. SNR')
set(gca,'Xscale','log')
subplot(1,2,2)
stem(expVars.quant, expVars.SNRstd./expVars.SNRAvg)
xlabel('Quant Time[s]')
ylabel('(SNR Std)/(SNR Avg)')
set(gca,'Xscale','log')

SNRVec(5) = expVars.SNRAvg(5); 
STDVec(5) = expVars.SNRstd(5);
timeVec(5) = 2;


% -> SNR vs Quant PMT 50Hz 2s
dir = 'D:\WavelengthScan\SNR-Vs.-Quant-PMT-2s-50Hz-2';
load(sprintf('%s/19-Jul-2021 20-28-53-PMT-2s-50Hz-ExpVars.mat', dir))

figure()
subplot(1,2,1)
errorbar(expVars.quant, expVars.SNRAvg, expVars.SNRstd)
xlabel('Quant Time[s]')
ylabel('Avg. SNR')
set(gca,'Xscale','log')
subplot(1,2,2)
stem(expVars.quant, expVars.SNRstd./expVars.SNRAvg)
xlabel('Quant Time[s]')
ylabel('(SNR Std)/(SNR Avg)')
set(gca,'Xscale','log')

SNRVec(6) = expVars.SNRAvg(5); 
STDVec(6) = expVars.SNRstd(5);
timeVec(6) = 2;


%SNR vs Scan Freq (quant = 0.02s)
normVec = sqrt(timeVec./2);
SNRVecNorm = SNRVec./normVec;
STDVecNorm = STDVec./normVec;
freqVec = [0.05, 0.5, 5, 10, 25, 50];
figure();
errorbar(freqVec, SNRVecNorm, STDVecNorm)
xlabel('Scan Frequency (Hz)')
ylabel('SNR')
set(gca,'Xscale','log')

% -> SNR vs Quant PMT 0Hz 2s
dir = 'D:\WavelengthScan\SNR-Vs.-Quant-PMT-2s-0Hz';
load(sprintf('%s/19-Jul-2021 17-55-56-PMT-2s-0Hz-ExpVars.mat', dir))

figure()
subplot(1,2,1)
errorbar(expVars.quant, expVars.SNRAvg, expVars.SNRstd)
xlabel('Quant Time[s]')
ylabel('Avg. SNR')
set(gca,'Xscale','log')
subplot(1,2,2)
stem(expVars.quant, expVars.SNRstd./expVars.SNRAvg)
xlabel('Quant Time[s]')
ylabel('(SNR Std)/(SNR Avg)')
set(gca,'Xscale','log')