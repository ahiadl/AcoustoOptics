close all;
clear all;
clc;

%----------------------------------------------------------%
%In this file:
%1. Measuring for 1, 2, and 4 cycles in pulse Acousto Optics results with
%   and without multiplexing
%2. Calculating the theoretical SNR value and displaying Results (theory vs measurement).
%
% There are already saved results therefore it is loaded
%----------------------------------------------------------%

%% init Acousto Optics Objects and Params

addpath(genpath(pwd));
% instrreset;

% Stages = stages('COM3');
% Stages.connect();
% 
% Stages.moveStageAbs([121, 91]);

ao = acoustoOptics();
% ao.init();

uVars = acoustoOptics.uVarsCreate();

uVars.c                 = 1445;
uVars.fSin              = 1.25e6;              
uVars.fTrain            = 19e3;
uVars.cycInPulse        = 4; 
uVars.channels          = 4; %update in digitizer
uVars.phantomDepth      = 4.2e-2;
uVars.distFromPhantom   = 6.4e-2;
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

uVars.gReq.validStruct.phi = true;
uVars.gReq.validStruct.usSignal = false;

uVars.exportRawData.rawData       = true;
uVars.exportRawData.netSignal     = false;
uVars.exportRawData.deMultiplexed = false;
uVars.exportRawData.reshape       = false;
uVars.exportRawData.FFT           = false;

uVars.fs.scanName        = [];
uVars.fs.resDirPath      = 'D:\Results';
uVars.fs.saveFullData    = false;
uVars.fs.saveFigs        = false;
uVars.fs.saveResults     = true;

%% Calc Multiplexing Orders
cyc = [1, 2, 4];

% fprintf("calculating Hadamards Order")
% i =1; j =1;
% while i <1000
%    S = createSMatrix(i, 'QR');
%    order(j) = size(S,1);
%    i = order(j) +1;
%    j = j+1;
% end
% save('./Scripts/hadamardOrders.mat', 'order');

load('./Scripts/hadamardOrders.mat');

samplesPerSin = uVars.fExtClk / uVars.fSin;
samplesInTrain = [order*cyc(1)*samplesPerSin;... 
                  order*cyc(2)*samplesPerSin;...
                  order*cyc(3)*samplesPerSin];

fTrains = uVars.fExtClk ./ samplesInTrain;

df = 0.1*(abs(fTrains(:, end) - fTrains(:, end-1)));

cyc4Idx = 3:1:12;
cyc2Idx = 4:2:22;
cyc1Idx = 7:2:25;
% 
% cyc4Idx = [3,12];
% cyc2Idx = [4,22];
% cyc1Idx = [25,7];

N = { order(cyc1Idx);...
      order(cyc2Idx);...
      order(cyc4Idx) };
  
f = { df(1) + fTrains(1, cyc1Idx);...
      df(2) + fTrains(2, cyc2Idx);...
      df(3) + fTrains(3, cyc4Idx); };

repPerConf = 30;

orders = N;
%% Measure
% for i =1:length(cyc) 
%     uVars.cycInPulse = cyc(i);
%     for j = 1:length(N{i})
%         uVars.fTrain = f{i}(j);
%         
%         %no Hadamard:
%         uVars.useHadamard = false;
%         uVars.useGPU = true;
%         name = sprintf("AcoustoOptics-%dcyc-woHad-%dn-20mW", cyc(i), N{i}(j));
%         uVars.fs.scanName = name;
%         
%         
%         ao.setMeasVars(uVars);
%         
%         numOfPosNaive = ao.algo.samples.numOfPos;
%        
%         ao.configPeripherals();
%         
%         fprintf(name);
%         for k = 1:repPerConf
%             ao.measureAndAnlayse();
%             disp(k);
%             pause(0.01);
%         end
%         
%         %with Hadamard:
%         uVars.useHadamard = true;
%         uVars.useGPU = false;
%         name = sprintf("AcoustoOptics-%dcyc-wHad-%dn-20mW", cyc(i), N{i}(j));
%         uVars.fs.scanName = name;
%         fprintf(name);
%         
%         ao.setMeasVars(uVars);
%         numOfPosHad = ao.algo.samples.numOfPos;
%         if (numOfPosHad ~= numOfPosNaive)
%             numOfPosHad = ao.algo.samples.numOfPos;
%             trueFreq = ao.algo.usSignal.fTrain;
%             
%             setpref('Internet','E_mail','ahiadlevi@gmail.com');
%             setpref('Internet','SMTP_Server','smtp.gmail.com');
%             
%             firstLine  = sprintf("Mismatch found in: %d cycles, Hadamard Order: %d.\n", cyc(i),  N{i}(j)); 
%             secondLine = sprintf("The desired frequency is: %.2f[Hz], While the calculated frequency is: %.2f[Hz].\n", f{i}(j), trueFreq); 
%             thirdLine  = sprintf("Hence, Hadamrd Num Of Pos is: %d, while Naive num of pos is: %d.\n", numOfPosHad, numOfPosNaive);
%             text = sprintf("%s%s%s", firstLine, secondLine, thirdLine);
%             
%             sendmail('sahiadl@campus.technion.ac.il','Positions Mismatch In Scan', text);
%             
%             beep
%             fprintf("Notice: Hadamard Positions are different from the naive Positions\n");
%         end
%         
%         ao.configPeripherals();
%         
%         for k = 1:repPerConf
%             ao.measureAndAnlayse();
%             pause(0.01);
%         end
%  
%     end
% end

%% Rename Folders
% moveFilesToDir('FellgettsAdvantage', '24-Jul-2019')

%% Save Phi structure
% for c = 1:length(cyc)
%     for i = 1:length(N{c})
%         for k=1:2
%             if k == 1        
%                 dirName = sprintf("D:/Results/FellgettsAdvantage/AcoustoOptics-%dcyc-woHad-%dn-20mW/", cyc(c), N{c}(i));
% %                 cyc4VarNaive(c,i) = load(sprintf("%sVars.mat", dirName));
% %                 fprintf("Loading: Cycles: %d, Order: %d, Hadamard: %d\n", cyc(c), N{c}(i), cyc4VarNaive(c,i).vars.algo.uVars.useHadamard);
%                 fprintf("Loading: Cycles: %d, Order: %d, Hadamard: %d\n", cyc(c), N{c}(i), 0);
%                 for j =1:repPerConf
%                     filename = sprintf("%s%d-Results",dirName, j);
% %                     cyc4Naive(c,i,j) = load(filename);
%                     cyc4Naive = load(filename);
% %                     cyc4PhiNaive{c,i}(j,:) = gather(cyc4Naive(c,i,j).res.phi);
%                     cyc4PhiNaive{c,i}(j,:) = cyc4Naive.res.phi;
%                     
%                     if j>1
%                       for m=0:log10(j-1)
%                           fprintf('\b'); % delete previous counter display
%                       end
%                     end
%                     fprintf('%d', j);
%                 end
%             else
%                 dirName = sprintf("D:/Results/FellgettsAdvantage/AcoustoOptics-%dcyc-wHad-%dn-20mW/", cyc(c),  N{c}(i));
% %                 cyc4VarHad(c,i) = load(sprintf("%sVars.mat", dirName));
% %                 cyc4VarHad(c,i) = load(sprintf("%sVars.mat", dirName));
% %                 fprintf("Loading: Cycles: %d, Order: %d, Hadamard: %d\n", cyc(c), N{c}(i), cyc4VarHad(c,i).vars.algo.uVars.useHadamard);
%                 fprintf("Loading: Cycles: %d, Order: %d, Hadamard: %d\n", cyc(c), N{c}(i), 1);
%                 for j =1:repPerConf
%                     filename = sprintf("%s%d-Results",dirName, j);
% %                      cyc4Had(c,i,j)= load(filename);
%                     cyc4Had= load(filename);
% %                      cyc4PhiHad{c,i}(j,:) = gather(cyc4Had(c,i,j).res.phi);
%                     cyc4PhiHad{c,i}(j,:) = cyc4Had.res.phi;
%                     
%                     if j>1
%                       for m=0:log10(j-1)
%                           fprintf('\b'); % delete previous counter display
%                       end
%                     end
%                     fprintf('%d', j);
%                     
%                 end
%             end
%         end
%     end
% end

% save('D:/Results/FellgettsAdvantage/phiHad.mat', 'cyc4PhiHad', '-v7.3');
% save('D:/Results/FellgettsAdvantage/phiNaive.mat', 'cyc4PhiNaive', '-v7.3');

%% save rawData structures

% for c = 1:length(cyc)
%     for i = 1:length(N{c})
%         dirName = sprintf("D:/ResultsToKeep/FellgettsAdvantage/AcoustoOptics-%dcyc-woHad-%dn-20mW/", cyc(c), N{c}(i));
%         fprintf("Loading: Cycles: %d, Order: %d, Hadamard: %d\n", cyc(c), N{c}(i), 0);
%         Vars = load(sprintf("%sVars.mat", dirName));
%         uVars.fTrain = Vars.vars.algo.uVars.fTrain;
%         uVars.cycInPulse = Vars.vars.algo.uVars.cycInPulse;
%         uVars.useHadamard = false;
%         uVars.useGPU = true;
%         ao.setMeasVars(uVars);
%         zVecUSRes{c,i} = Vars.vars.algo.len.zVecUSRes;
        
%         for j =1:repPerConf
%             filename = sprintf("%s%d-Results",dirName, j);
%             newData = load(filename);
%             ao.algo.setRawData(newData.res.rawData);
%             res = ao.algo.analyse();
%             dataNaive{c,i}(j,:) = res.phi;
%             fftNaive {c,i}(j,:,:,:) = res.phiChCmplx; %[quant, ch, pos]
%             if j>1
%               for m=0:log10(j-1)
%                   fprintf('\b'); % delete previous counter display
%               end
%             end
%             fprintf('%d', j);
%         end
%         dirName = sprintf("D:/ResultsToKeep/FellgettsAdvantage/AcoustoOptics-%dcyc-wHad-%dn-20mW/", cyc(c), N{c}(i));
%         fprintf("Loading: Cycles: %d, Order: %d, Hadamard: %d\n", cyc(c), N{c}(i), 1);
%         uVars.useHadamard = true;
%         uVars.useGPU = false;
%         ao.algo.setRawData(newData.res.rawData);
%         ao.setMeasVars(uVars);
%         for j =1:repPerConf
%             filename = sprintf("%s%d-Results",dirName, j);
%             newData = load(filename);
%             ao.algo.setRawData(newData.res.rawData);
%             res = ao.algo.analyse();
%             dataHad{c,i}(j,:) = res.phi;
%             fftHad {c,i}(j,:,:,:) = res.phiChCmplx;
%             if j>1
%               for m=0:log10(j-1)
%                   fprintf('\b'); % delete previous counter display
%               end
%             end
%             fprintf('%d', j);
%         end
%     end
% end
% save('D:/ResultsToKeep/FellgettsAdvantage/zVecUSRes.mat', 'zVecUSRes', '-v7.3');
% save('D:/ResultsToKeep/FellgettsAdvantage/rawDataHad.mat', 'dataHad', '-v7.3');
% save('D:/ResultsToKeep/FellgettsAdvantage/rawDataNaive.mat', 'dataNaive', '-v7.3');
% save('D:/ResultsToKeep/FellgettsAdvantage/fftHad.mat', 'fftHad', '-v7.3');
% save('D:/ResultsToKeep/FellgettsAdvantage/fftNaive.mat', 'fftNaive', '-v7.3');

%% Analysis
load('D:/ResultsToKeep/FellgettsAdvantage/zVecUSRes.mat')
load('D:/ResultsToKeep/FellgettsAdvantage/phiHad.mat')
load('D:/ResultsToKeep/FellgettsAdvantage/phiNaive.mat')
load('D:/ResultsToKeep/FellgettsAdvantage/fftHad.mat')
load('D:/ResultsToKeep/FellgettsAdvantage/fftNaive.mat')
load('D:/ResultsToKeep/FellgettsAdvantage/rawDataHad.mat')
load('D:/ResultsToKeep/FellgettsAdvantage/rawDataNaive.mat')
%%
% dims

idx1 = 29;
idx2 = 15;
idx4 = 8;

maxIdx = [idx1, idx2, idx4];
postMaxPts = [8, 11, 7];

samp = 1;
ch = 3;
quant = 2;
pos = 4;

for i = 1:length(cyc)
    for j = 1:length(N{i})

        % Had
        fftMat = fftHad{i,j}(:,:,:,:);
        absMeanFFTtmp   = mean(abs(mean(fftMat,     ch)), quant);
        AbsFFTtmp       = mean(mean(abs(fftMat),    ch),  quant);
        AbsSqFFTtmp     = mean(mean(abs(fftMat).^2, ch),  quant);

        hadAbsMeanFFT{i,j}   = absMeanFFTtmp;
        hadMeanAbsFFT{i,j}   = AbsFFTtmp;
        hadMeanAbsSqFFT{i,j} = AbsSqFFTtmp;

        %Average over different samples
        expHadAbsMeanFFT{i,j}   = squeeze(mean(absMeanFFTtmp, samp));
        stdHadAbsMeanFFT{i,j}   = squeeze(std(absMeanFFTtmp, 0, samp));
        
        expHadMeanAbsFFT{i,j}   = squeeze(mean(AbsFFTtmp, samp));
        stdHadMeanAbsFFT{i,j}   = squeeze(std(AbsFFTtmp, 0, samp));
        
        expHadMeanAbsSqFFT{i,j} = squeeze(mean(AbsSqFFTtmp, samp));
        stdHadMeanAbsSqFFT{i,j} = squeeze(std(AbsSqFFTtmp, 0, samp));

        if j == 1
            sigStdHadAbsMean(i,j)  = stdHadAbsMeanFFT{i,j}(maxIdx(i));
            bkgStdHadAbsMean(i,j)  = stdHadAbsMeanFFT{i,j}(1);
            sigMeanHadAbsMean(i,j) = expHadAbsMeanFFT{i,j}(maxIdx(i));
            bkgMeanHadAbsMean(i,j) = expHadAbsMeanFFT{i,j}(1);
            
            sigStdHadAbs(i,j)  = stdHadMeanAbsFFT{i,j}(maxIdx(i));
            bkgStdHadAbs(i,j)  = stdHadMeanAbsFFT{i,j}(1);
            sigMeanHadAbs(i,j) = expHadMeanAbsFFT{i,j}(maxIdx(i));
            bkgMeanHadAbs(i,j) = expHadMeanAbsFFT{i,j}(1);
            
            sigStdHadAbsSq(i,j)  = stdHadMeanAbsSqFFT{i,j}(maxIdx(i));
            bkgStdHadAbsSq(i,j)  = stdHadMeanAbsSqFFT{i,j}(1);
            sigMeanHadAbsSq(i,j) = expHadMeanAbsSqFFT{i,j}(maxIdx(i));
            bkgMeanHadAbsSq(i,j) = expHadMeanAbsSqFFT{i,j}(1);

            SNRHadAbsMean(i,j)   = abs( expHadAbsMeanFFT{i,j}(maxIdx(i)) - expHadAbsMeanFFT{i,j}(1)) / stdHadAbsMeanFFT{i,j}(1);
            SNRHadMeanAbs(i,j)   = abs( expHadMeanAbsFFT{i,j}(maxIdx(i)) - expHadMeanAbsFFT{i,j}(1)) / stdHadMeanAbsFFT{i,j}(1);
            SNRHadMeanAbsSq(i,j) = abs( expHadMeanAbsSqFFT{i,j}(maxIdx(i)) - expHadMeanAbsSqFFT{i,j}(1)) / stdHadMeanAbsSqFFT{i,j}(1);

            SBRHadAbsMean(i,j)   = abs( expHadAbsMeanFFT{i,j}(maxIdx(i)) - expHadAbsMeanFFT{i,j}(1)) / expHadAbsMeanFFT{i,j}(1);
            SBRHadMeanAbs(i,j)   = abs( expHadMeanAbsFFT{i,j}(maxIdx(i)) - expHadMeanAbsFFT{i,j}(1)) / expHadMeanAbsFFT{i,j}(1);
            SBRHadMeanAbsSq(i,j) = abs( expHadMeanAbsSqFFT{i,j}(maxIdx(i)) - expHadMeanAbsSqFFT{i,j}(1)) / expHadMeanAbsSqFFT{i,j}(1);

        else
            
            tmpBkgStdVec = absMeanFFTtmp(:,:,:, (maxIdx(i)+postMaxPts(i)):end);
            tmpBkgStdVec = tmpBkgStdVec(:);
            bkgStdAbsMean  = std(tmpBkgStdVec);
            bkgMeanAbsMean = mean(tmpBkgStdVec);
            
            tmpBkgStdVec = AbsFFTtmp(:,:,:, (maxIdx(i)+postMaxPts(i)):end);
            tmpBkgStdVec = tmpBkgStdVec(:);
            bkgStdAbs  = std(tmpBkgStdVec);
            bkgMeanAbs = mean(tmpBkgStdVec);

            tmpBkgStdVec = AbsSqFFTtmp(:,:,:, (maxIdx(i)+postMaxPts(i)):end);
            tmpBkgStdVec = tmpBkgStdVec(:);
            bkgStdAbsSq = std(tmpBkgStdVec);
            bkgMeanAbsSq = mean(tmpBkgStdVec);

            
            sigStdHadMean1(i,j) = stdHadAbsMeanFFT{i,j}(maxIdx(i));
            bkgStdHadMean(i,j)  = bkgStdAbsMean;
            sigMeanHadMean(i,j) = expHadAbsMeanFFT{i,j}(maxIdx(i));
            bkgMeanHadMean(i,j) = bkgMeanAbsMean;
            
            sigStdHadAbs(i,j)  = stdHadMeanAbsFFT{i,j}(maxIdx(i));
            bkgStdHadAbs(i,j)  = bkgStdAbs;
            sigMeanHadAbs(i,j) = expHadMeanAbsFFT{i,j}(maxIdx(i));
            bkgMeanHadAbs(i,j) = bkgMeanAbs;
            
            sigStdHadAbsSq(i,j)  = stdHadMeanAbsSqFFT{i,j}(maxIdx(i));
            bkgStdHadAbsSq(i,j)  = bkgStdAbsSq;
            sigMeanHadAbsSq(i,j) = expHadMeanAbsSqFFT{i,j}(maxIdx(i));
            bkgMeanHadAbsSq(i,j) = bkgMeanAbsSq;

            SNRHadAbsMean(i,j)   = abs( expHadAbsMeanFFT{i,j}(maxIdx(i)) - expHadAbsMeanFFT{i,j}(end)) / stdHadAbsMeanFFT{i,j}(end);
            SNRHadMeanAbs(i,j)   = abs( expHadMeanAbsFFT{i,j}(maxIdx(i)) - expHadMeanAbsFFT{i,j}(end)) / stdHadMeanAbsFFT{i,j}(end);
            SNRHadMeanAbsSq(i,j) = abs( expHadMeanAbsSqFFT{i,j}(maxIdx(i)) - expHadMeanAbsSqFFT{i,j}(end)) / stdHadMeanAbsSqFFT{i,j}(end);

            SBRHadAbsMean(i,j)   = abs( expHadAbsMeanFFT{i,j}(maxIdx(i)) - expHadAbsMeanFFT{i,j}(end)) / expHadAbsMeanFFT{i,j}(end);
            SBRHadMeanAbs(i,j)   = abs( expHadMeanAbsFFT{i,j}(maxIdx(i)) - expHadMeanAbsFFT{i,j}(end)) / expHadMeanAbsFFT{i,j}(end);
            SBRHadMeanAbsSq(i,j) = abs( expHadMeanAbsSqFFT{i,j}(maxIdx(i)) - expHadMeanAbsSqFFT{i,j}(end)) / expHadMeanAbsSqFFT{i,j}(end);

        end

        %Naive
        fftMat = fftNaive{i,j}(:,:,:,:);
        absMeanFFTtmp   = mean(abs(mean(fftMat, ch)),quant);
        AbsFFTtmp       = mean(mean(abs(fftMat),ch), quant);
        AbsSqFFTtmp     = mean(mean(abs(fftMat).^2, ch), quant);

        naiveAbsMeanFFT{i,j}   = absMeanFFTtmp;
        naiveMeanAbsFFT{i,j}   = AbsFFTtmp;
        naiveMeanAbsSqFFT{i,j} = AbsSqFFTtmp;
        
        % Average over samples
        expNaiveAbsMeanFFT{i,j}   = squeeze(mean(absMeanFFTtmp, samp));
        stdNaiveAbsMeanFFT{i,j}   = squeeze(std(absMeanFFTtmp, 0, samp));

        expNaiveMeanAbsFFT{i,j}   = squeeze(mean(AbsFFTtmp, samp));
        stdNaiveMeanAbsFFT{i,j}   = squeeze(std(AbsFFTtmp, 0, samp));

        expNaiveMeanAbsSqFFT{i,j} = squeeze(mean(AbsSqFFTtmp, samp));
        stdNaiveMeanAbsSqFFT{i,j} = squeeze(std(AbsSqFFTtmp, 0, samp));

        if j == 1
            sigStdNaiveAbsMean(i,j)  = stdNaiveAbsMeanFFT{i,j}(maxIdx(i));
            bkgStdNaiveAbsMean(i,j)  = stdNaiveAbsMeanFFT{i,j}(1);
            sigMeanNaiveAbsMean(i,j) = expNaiveAbsMeanFFT{i,j}(maxIdx(i));
            bkgMeanNaiveAbsMean(i,j) = expNaiveAbsMeanFFT{i,j}(1);
            
            sigStdNaiveAbs(i,j)  = stdNaiveMeanAbsFFT{i,j}(maxIdx(i));
            bkgStdNaiveAbs(i,j)  = stdNaiveMeanAbsFFT{i,j}(1);
            sigMeanNaiveAbs(i,j) = expNaiveMeanAbsFFT{i,j}(maxIdx(i));
            bkgMeanNaiveAbs(i,j) = expNaiveMeanAbsFFT{i,j}(1);
            
            sigStdNaiveAbsSq(i,j)  = stdNaiveMeanAbsSqFFT{i,j}(maxIdx(i));
            bkgStdNaiveAbsSq(i,j)  = stdNaiveMeanAbsSqFFT{i,j}(1);
            sigMeanNaiveAbsSq(i,j) = expNaiveMeanAbsSqFFT{i,j}(maxIdx(i));
            bkgMeanNaiveAbsSq(i,j) = expNaiveMeanAbsSqFFT{i,j}(1);

            SNRNaiveAbsMean(i,j)   = abs( expNaiveAbsMeanFFT{i,j}(maxIdx(i)) - expNaiveAbsMeanFFT{i,j}(1)) / stdNaiveAbsMeanFFT{i,j}(1);
            SNRNaiveMeanAbs(i,j)   = abs( expNaiveMeanAbsFFT{i,j}(maxIdx(i)) - expNaiveMeanAbsFFT{i,j}(1)) / stdNaiveMeanAbsFFT{i,j}(1);
            SNRNaiveMeanAbsSq(i,j) = abs( expNaiveMeanAbsSqFFT{i,j}(maxIdx(i)) - expNaiveMeanAbsSqFFT{i,j}(1)) / stdNaiveMeanAbsSqFFT{i,j}(1);

            SBRNaiveAbsMean(i,j)   = abs( expNaiveAbsMeanFFT{i,j}(maxIdx(i)) - expNaiveAbsMeanFFT{i,j}(1)) / expNaiveAbsMeanFFT{i,j}(1);
            SBRNaiveMeanAbs(i,j)   = abs( expNaiveMeanAbsFFT{i,j}(maxIdx(i)) - expNaiveMeanAbsFFT{i,j}(1)) / expNaiveMeanAbsFFT{i,j}(1);
            SBRNaiveMeanAbsSq(i,j) = abs( expNaiveMeanAbsSqFFT{i,j}(maxIdx(i)) - expNaiveMeanAbsSqFFT{i,j}(1)) / expNaiveMeanAbsSqFFT{i,j}(1);

        else
            tmpBkgStdVec = absMeanFFTtmp(:,:,:, (maxIdx(i)+postMaxPts(i)):end);
            tmpBkgStdVec = tmpBkgStdVec(:);
            bkgStdAbsMean  = std(tmpBkgStdVec);
            bkgMeanAbsMean = mean(tmpBkgStdVec);
            
            tmpBkgStdVec = AbsFFTtmp(:,:,:, (maxIdx(i)+postMaxPts(i)):end);
            tmpBkgStdVec = tmpBkgStdVec(:);
            bkgStdAbs  = std(tmpBkgStdVec);
            bkgMeanAbs = mean(tmpBkgStdVec);

            tmpBkgStdVec = AbsSqFFTtmp(:,:,:, (maxIdx(i)+postMaxPts(i)):end);
            tmpBkgStdVec = tmpBkgStdVec(:);
            bkgStdAbsSq = std(tmpBkgStdVec);
            bkgMeanAbsSq = mean(tmpBkgStdVec);
            
            sigStdNaiveAbsMean(i,j)  = stdNaiveAbsMeanFFT{i,j}(maxIdx(i));
            bkgStdNaiveAbsMean(i,j)  = bkgStdAbsMean;
            sigMeanNaiveAbsMean(i,j) = expNaiveAbsMeanFFT{i,j}(maxIdx(i));
            bkgMeanNaiveAbsMean(i,j) = bkgMeanAbsMean;
            
            sigStdNaiveAbs(i,j)  = stdNaiveMeanAbsFFT{i,j}(maxIdx(i));
            bkgStdNaiveAbs(i,j)  = bkgStdAbs;
            sigMeanNaiveAbs(i,j) = expNaiveMeanAbsFFT{i,j}(maxIdx(i));
            bkgMeanNaiveAbs(i,j) = bkgMeanAbs;

            sigStdNaiveAbsSq(i,j)  = stdNaiveMeanAbsSqFFT{i,j}(maxIdx(i));
            bkgStdNaiveAbsSq(i,j)  = bkgStdAbsSq;
            sigMeanNaiveAbsSq(i,j) = expNaiveMeanAbsSqFFT{i,j}(maxIdx(i));
            bkgMeanNaiveAbsSq(i,j) = bkgMeanAbsSq;

            SNRNaiveAbsMean(i,j)   = abs( expNaiveAbsMeanFFT{i,j}(maxIdx(i)) - expNaiveAbsMeanFFT{i,j}(end)) / stdNaiveAbsMeanFFT{i,j}(end);
            SNRNaiveMeanAbs(i,j)   = abs( expNaiveMeanAbsFFT{i,j}(maxIdx(i)) - expNaiveMeanAbsFFT{i,j}(end)) / stdNaiveMeanAbsFFT{i,j}(end);
            SNRNaiveMeanAbsSq(i,j) = abs( expNaiveMeanAbsSqFFT{i,j}(maxIdx(i)) - expNaiveMeanAbsSqFFT{i,j}(end)) / stdNaiveMeanAbsSqFFT{i,j}(end);

            SBRNaiveAbsMean(i,j)   = abs( expNaiveAbsMeanFFT{i,j}(maxIdx(i)) - expNaiveAbsMeanFFT{i,j}(end)) / expNaiveAbsMeanFFT{i,j}(end);
            SBRNaiveMeanAbs(i,j)   = abs( expNaiveMeanAbsFFT{i,j}(maxIdx(i)) - expNaiveMeanAbsFFT{i,j}(end)) / expNaiveMeanAbsFFT{i,j}(end);
            SBRNaiveMeanAbsSq(i,j) = abs( expNaiveMeanAbsSqFFT{i,j}(maxIdx(i)) - expNaiveMeanAbsSqFFT{i,j}(end)) / expNaiveMeanAbsSqFFT{i,j}(end);

        end
    end
end

%%
FASNRAbsMean   = SNRHadAbsMean./SNRNaiveAbsMean;
FASNRMeanAbs   = SNRHadMeanAbs./SNRNaiveMeanAbs;
FASNRMeanAbsSq = SNRHadMeanAbsSq./SNRNaiveMeanAbsSq;

FASBRAbsMean = SBRHadAbsMean./SBRNaiveAbsMean;
FASBRMeanAbs = SBRHadMeanAbs./SBRNaiveMeanAbs;
FASBRMeanAbsSq = SBRHadMeanAbsSq./SBRNaiveMeanAbsSq;

FASNRStdSigAbs   = sigStdNaiveAbs ./ sigStdHadAbs;
FASNRStdSigAbsSq = sigStdNaiveAbsSq ./ sigStdHadAbsSq;
FASNRStdBkgAbs   = bkgStdNaiveAbs ./ bkgStdHadAbs;
FASNRStdBkgAbsSq = bkgStdNaiveAbsSq ./ bkgStdHadAbsSq;

%% save Post Analysis Results

save('D:/ResultsToKeep/FellgettsAdvantage/analysisResultsHad.mat', ...
    'hadMeanAbsFFT', 'expHadMeanAbsFFT', 'stdHadMeanAbsFFT',...
    'sigStdHadAbs', 'bkgStdHadAbs', 'sigMeanHadAbs', 'bkgMeanHadAbs',...
    'SNRHadMeanAbs', '-v7.3');

save('D:/ResultsToKeep/FellgettsAdvantage/analysisResultsNaive.mat', ...
    'naiveMeanAbsFFT', 'expNaiveMeanAbsFFT', 'stdNaiveMeanAbsFFT',...
    'sigStdNaiveAbs', 'bkgStdNaiveAbs', 'sigMeanNaiveAbs', 'bkgMeanNaiveAbs',...
    'SNRNaiveMeanAbs', '-v7.3');

save('D:/ResultsToKeep/FellgettsAdvantage/postAnalysisAdvantageResults.mat', ...
    'FASNRAbsMean', 'FASNRMeanAbs', 'FASNRMeanAbsSq',...
    'FASBRAbsMean', 'FASBRMeanAbs', 'FASBRMeanAbsSq',...
    'FASNRStdSigAbs', 'FASNRStdSigAbsSq', 'FASNRStdBkgAbs', 'FASNRStdBkgAbsSq',...
    'N','-v7.3');


%% Phd Figures
figure()
plot(zVecUSRes{3,1}, squeeze(naiveMeanAbsFFT{3,1}(1,:,:,:)), '-^'); hold on;
plot(zVecUSRes{1,1},squeeze(naiveMeanAbsFFT{1,1}(1,:,:,:))-3e-3,'-+');
plot(zVecUSRes{3,1},squeeze(hadMeanAbsFFT{3,1}(1,:,:,:)), '-^'); hold on;
plot(zVecUSRes{1,1},squeeze(hadMeanAbsFFT{1,1}(1,:,:,:)),'-+');
set(gca, 'FontSize', 14)
legend("Trad 4", "Trad 1", "Coded 4", "Coded 1", 'FontSize', 14, 'Location', 'northwest')
xlabel("Z[mm]", 'FontSize', 14)
ylabel("Relative Fluence", 'FontSize', 14)


% figure()
% for i = 1:3
%     subplot(1,3,i)
%     plot(N{i}, (N{i}+1)./(2*sqrt(N{i})), '-+'); hold on;
%     plot(N{i}, FASNRStdBkgAbs(i,:), '-v');
% %     plot(N{i}, FASNRStdBkgAbsSq(i,:), '-x');
%     legend("theory", "Measurement");
%     title(sprintf("FA STD_{bkg} for %d Cycles", cyc(i)))
% end

figure()
i=3;
plot(N{i}, (N{i}+1)./(2*sqrt(N{i})), '-+'); hold on;
plot(N{i}, FASNRStdBkgAbs(i,:), '-v');
set(gca, 'FontSize', 14)
legend("Theory", "Measurement", 'FontSize', 14, 'Location', 'northwest');
xlabel("Code Order", 'FontSize', 14)
ylabel("Multiplexing Advantage", 'FontSize', 14)


%%
figure()
% plot(squeeze(hadAbsMeanFFT{3,1}(1,:,:,:))); hold on
subplot(1,2,1)
plot(squeeze(hadMeanAbsFFT{3,1}(1,:,:,:)));
title("Abs");
subplot(1,2,2)
plot(squeeze(hadMeanAbsSqFFT{3,1}(1,:,:,:)));
title("Abs^{2}");
% legend("abs(mean)", "mean(abs)", "mean(abs.^2))")

figure()
plot(N{3}, (N{3}+1)./(2*sqrt(N{3})), '-+'); hold on;
% plot(N{3},FASNRAbsMean(3,:), '-v');
plot(N{3},FASNRMeanAbs(3,:), '-o');
plot(N{3},FASNRMeanAbsSq(3,:), '-x');
title("SNR")
legend("theory", "Abs", "Abs^2")

figure()
plot(N{3}, (N{3}+1)./(2*sqrt(N{3})), '-+'); hold on;
% plot(N{3},FASBRAbsMean(3,:), '-v');
plot(N{3},FASBRMeanAbs(3,:), '-o');
plot(N{3},FASBRMeanAbsSq(3,:), '-x');
title("SBR")
legend("theory", "Abs", "Abs^2")

figure()
for i = 1:3
    subplot(1,3,i)
    plot(N{i}, (N{i}+1)./(2*sqrt(N{i})), '-+'); hold on;
    plot(N{i},FASBRMeanAbs(i,:), '-v');
    plot(N{i},FASBRMeanAbsSq(i,:), '-v');
    legend("theory", "Abs", "Abs^{2}");
    title(sprintf("SBR for %d Cycles in pulse", cyc(i)))
end


figure()
for i = 1:3
    subplot(1,3,i)
    plot(N{i}, (N{i}+1)./(2*sqrt(N{i})), '-+'); hold on;
    plot(N{i},FASNRMeanAbs(i,:), '-v');
    plot(N{i},FASNRMeanAbsSq(i,:), '-v');
    legend("theory", "measurement");
    title(sprintf("SNR for %d Cycles in pulse", cyc(i)))
end

FASNRStdSigAbs   = sigStdNaiveAbs ./ sigStdHadAbs;
FASNRStdSigAbsSq = sigStdNaiveAbsSq ./ sigStdHadAbsSq;
FASNRStdBkgAbs   = bkgStdNaiveAbs ./ bkgStdHadAbs;
FASNRStdBkgAbsSq = bkgStdNaiveAbsSq ./ bkgStdHadAbsSq;

figure()
for i = 1:3
    subplot(1,3,i)
    plot(N{i}, (N{i}+1)./(2*sqrt(N{i})), '-+'); hold on;
    plot(N{i}, FASNRStdBkgAbs(i,:), '-v');
    plot(N{i}, FASNRStdBkgAbsSq(i,:), '-x');
    legend("theory", "Abs" , "Abs^{2}");
    title(sprintf("FA STD_{bkg} for %d Cycles in pulse", cyc(i)))
end

figure()
for i = 1:3
    subplot(1,3,i)
    plot(N{i}, (N{i}+1)./(2*sqrt(N{i})), '-+'); hold on;
    plot(N{i}, FASNRStdSigAbs(i,:), '-v');
    plot(N{i}, FASNRStdSigAbsSq(i,:), '-x');
    legend("theory", "Abs" , "Abs^{2}");
    title(sprintf("FA STD_{sig} for %d Cycles in pulse", cyc(i)))
end


figure()
for i = 1:3
    subplot(1,3,i)
    plot(N{i}, (N{i}+1)./(2*sqrt(N{i})), '-+'); hold on;
    plot(N{i}, FASNRStdBkgAbs(i,:), '-v');
    plot(N{i}, FASNRStdBkgAbsSq(i,:), '-x');
    legend("theory", "Abs", "Abs^2");
    title(sprintf("FA STD_{bkg} for %d Cycles", cyc(i)))
end

figure()
for i = 1:3
    subplot(1,3,i)
    plot(N{i}, (N{i}+1)./(2*sqrt(N{i})), '-+'); hold on;
    plot(N{i}, FASNRStdBkgAbs(i,:), '-v');
%     plot(N{i}, FASNRStdBkgAbsSq(i,:), '-x');
    legend("theory", "Measurement");
    title(sprintf("FA STD_{bkg} for %d Cycles", cyc(i)))
end
%%
% figure();
% subplot(1,3,1)
% plot(cyc4PhiHad{1,1}(1,:), '-x'); hold on;
% plot(cyc4PhiNaive{1,1}(1,:), '-x');
% title("1 Cycles Average Phi");
% legend('Hadamard', 'Naive')
% subplot(1,3,2)
% plot(cyc4PhiHad{2,1}(1,:), '-x'); hold on;
% plot(cyc4PhiNaive{2,1}(1,:), '-x');
% title("2 Cycles Average Phi");
% legend('Hadamard', 'Naive');
% subplot(1,3,3)
% plot(cyc4PhiHad{3,1}(1,:), '-x'); hold on;
% plot(cyc4PhiNaive{3,1}(1,:), '-x');
% title("4 Cycle Average Phi");
% legend('Hadamard', 'Naive');

idx1 = 29;
idx2 = 15;
idx4 = 8;

maxIdx = [idx1, idx2, idx4];

prePoints = [20, 9, 5];
postMaxPts = [8, 11, 7];

hFig = figure();
ax = axes();
for c = 1:length(cyc)
    for i=1:length(N{c})
        meanSigHad = mean(cyc4PhiHad{c,i}, 1);
        meanSigNaive = mean(cyc4PhiNaive{c,i}, 1);
        
%         figure()
%         plot(meanSigHad); hold on;
%         plot(meanSigNaive);
%         title(sprintf("Cyc: %d, Order: %d", c, N{c}(i)));

        hadMaxAvg(c,i) = mean(cyc4PhiHad{c,i}(:, maxIdx(c)));
        hadMaxStd(c,i) = std(cyc4PhiHad{c,i}(:, maxIdx(c)));
        if i < 2 
            hadMinAvg(c,i) = mean(mean(cyc4PhiHad{c,i}(:, 1:prePoints(c)) ,2));
            hadMinTail     = cyc4PhiHad{c,i}(:, 1:prePoints(c));
            hadMinStd(c,i) = std(hadMinTail(:));
%             hadMinStd(c,i) = mean(std(cyc4PhiHad{c,i}(:, 1:prePoints(c)), 0, 2));
        else
            hadMinAvg(c,i) = mean(mean(cyc4PhiHad{c,i}(:, (maxIdx(c)+postMaxPts(c)):end) ,2));
            hadMinTail     = cyc4PhiHad{c,i}(:, (maxIdx(c)+postMaxPts(c)):end);
            hadMinStd(c,i) = std(hadMinTail(:));
%             hadMinStd(c,i) = mean(std(cyc4PhiHad{c,i}(:, (maxIdx(c)+postMaxPts(c)):end), 0, 2));
        end
        hadAVGDynamicRange =  abs(hadMaxAvg(c,i) - hadMinAvg(c,i))/hadMinAvg(c,i);
        hadSNRMax  = hadMaxAvg(c,i)/hadMaxStd(c,i);
        hadSNRMMin = hadMinAvg(c,i)/hadMinStd(c,i);
        hadSNR = (hadMaxAvg(c,i)-hadMinAvg(c,i))/hadMinStd(c,i);
        
        naiveMaxAvg(c,i) = mean(cyc4PhiNaive{c,i}(:, maxIdx(c)));
        naiveMaxStd(c,i) = std(cyc4PhiNaive{c,i}(:, maxIdx(c)));
        if i < 2
            naiveMinAvg(c,i) = mean(mean(cyc4PhiNaive{c,i}(:, 1:prePoints(c)), 2));
            naiveMinTail     = cyc4PhiNaive{c,i}(:, 1:prePoints(c));
            naiveMinStd(c,i) = std(naiveMinTail(:));
%             naiveMinStd(c,i) = std(cyc4PhiNaive{c,i}(:, 1:prePoints(c)));
        else
            naiveMinAvg(c,i) = mean(mean(cyc4PhiNaive{c,i}(:, (maxIdx(c)+postMaxPts(c)):end), 2));
            naiveMinTail     = cyc4PhiNaive{c,i}(:, (maxIdx(c)+postMaxPts(c)):end);
            naiveMinStd(c,i) = std(naiveMinTail(:));
%             naiveMinStd(c,i) = mean(std(cyc4PhiNaive{c,i}(:, (maxIdx(c)+postMaxPts(c)):end), 0, 2));
        end    
        naiveAVGDynamicRange =  abs(naiveMaxAvg(c,i) - naiveMinAvg(c,i))/naiveMinAvg(c,i);
        naiveSNRMax = naiveMaxAvg(c,i) / naiveMinAvg(c,i);
        naiveSNRMin = naiveMinAvg(c,i)/ naiveMinStd(c,i);
        naiveSNR = (naiveMaxAvg(c,i) - naiveMinAvg(c,i))/naiveMinStd(c,i);
        

%         dynamicRangeFactor{c}(i) = hadAVGDynamicRange/naiveAVGDynamicRange;

%         dynamicRangeFactor{c}(i) = naiveMinAvg(c,i) / hadMinAvg(c,i);
        dynamicRangeFactor{c}(i) = naiveMinStd(c,i) / hadMinStd(c,i); 
%         dynamicRangeFactor{c}(i) = naiveMaxAvg(c,i) / hadMaxAvg(c,i);
%         dynamicRangeFactor{c}(i) = naiveMaxStd(c,i) / hadMaxStd(c,i);
% 
%         dynamicRangeFactor{c}(i) = naiveMinStd(c,i) / hadMinStd(c,i);
% 
%         dynamicRangeFactor{c}(i) = hadSNRMax/naiveSNRMax;
%         dynamicRangeFactor{c}(i) = hadSNRMMin / naiveSNRMin;
%         dynamicRangeFactor{c}(i) = hadSNR / naiveSNR;

    end
    set(0, 'currentfigure',  hFig);
    subplot(1,length(cyc), c);
    plot(N{c}, (N{c}+1)./(2*sqrt(N{c})), '-+'); hold on; % N{c}/6
    plot(N{c}, dynamicRangeFactor{c}, '-v'); hold off;
    legend('Theory', 'Measurement');
    title(sprintf("%d Cycles In Pulse", cyc(c)));
    xlabel('Code Order')
    if(c==1)
%         ylabel("E(x)_{naive} / E(x)_{had}")
    end
end

%% Analyse Speckle Decorrelation
ylims = [2.5, 14]*1e-3;

figure();
subplot(1,3,1)
errorbar(N{1}, hadMaxAvg(1,:), hadMaxStd(1,:)); hold on;
errorbar(N{1}, naiveMaxAvg(1,:), naiveMaxStd(1,:));
ylim(ylims)
title('1 Cycle In Pulse');
xlabel('Hadamard Order \ Pulse in Z Axes')
legend('Hadamard', 'Naive');

subplot(1,3,2)
errorbar(N{2}, hadMaxAvg(2,:), hadMaxStd(2,:)); hold on;
errorbar(N{2}, naiveMaxAvg(2,:), naiveMaxStd(2,:));
ylim(ylims)
title('2 Cycle In Pulse');
xlabel('Hadamard Order \ Pulse in Z Axes')
legend('Hadamard', 'Naive');

subplot(1,3,3)
errorbar(N{3}, hadMaxAvg(3,:), hadMaxStd(3,:)); hold on;
errorbar(N{3}, naiveMaxAvg(3,:), naiveMaxStd(3,:));
title('4 Cycle In Pulse');
ylim(ylims)
legend('Hadamard', 'Naive');
xlabel('Hadamard Order \ Pulse in Z Axes')
