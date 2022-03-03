projPath = "D:\Results\20-May-2020 21-01-27-DV 2-200";
vars = load(sprintf("%s/DeepView-Vars.mat", projPath));
dvRes = load(sprintf("%s/DeepView-Results.mat", projPath));

loadRes(vars.time.numOfTimeFrames,vars.s2D.grid.firstIdxLen) = struct('qAvgChFFT', [], 'unFittedFFT', [], 'unFittedFFTShift', [], 'fitModel', [], 'fittedFFT', [], 'phi', []);
for i =1:vars.time.numOfTimeFrames
    aoVars(i)    = load(sprintf("%s/TimeResults/T=%d/AO-Vars.mat", projPath, vars.time.timeVec(i)));
    s2DVars(i)   = load(sprintf("%s/TimeResults/T=%d/2DScan-T=%d-Vars.mat", projPath, vars.time.timeVec(i), vars.time.timeVec(i)));
    for j = 1:vars.s2D.grid.firstIdxLen
        loadRes(i,j) = load(sprintf("%s/TimeResults/T=%d/AOResults/AO-R-1-X-%.2f-Results.mat", projPath, vars.time.timeVec(i), vars.s2D.grid.firstVec(j)));
    end
end
idx = 20;

res2Norm = dvRes.phi(:,1,idx) - min(dvRes.phi(:,1,idx));
res2Norm = res2Norm / max(res2Norm);

res200Norm = dvRes.phi(:,2,idx) - min(dvRes.phi(:,2,idx));
res200Norm = res200Norm / max(res200Norm);

figure();
subplot(1,2,1)
plot(vars.s2D.grid.firstVec, dvRes.phi(:,1,idx)); hold on
plot(vars.s2D.grid.firstVec, dvRes.phi(:,2,idx)); hold on
legend("T=2", "T=200")
xlabel("X[mm]")
ylabel("\phi")
title("fitting after channel mean")
subplot(1,2,2)
plot(vars.s2D.grid.firstVec, db(res2Norm)); hold on
plot(vars.s2D.grid.firstVec, db(res200Norm)); hold on
legend("T=2", "T=200")
xlabel("X[mm]")
ylabel("\phi [dB]")
title("fitting after channel mean (dB)")

figure()
subplot(1,2,1)
imagesc(aoVars(1).measVars.algo.len.zVecUSRes*1e3, vars.s2D.grid.firstVec, squeeze(dvRes.phi(:,1,:)))
title("2s")
xlabel("Z[mm]")
ylabel("X[mm]")
axis equal tight
colorbar
subplot(1,2,2)
imagesc(aoVars(1).measVars.algo.len.zVecUSRes*1e3, vars.s2D.grid.firstVec, squeeze(dvRes.phi(:,2,:)))
title("200s")
xlabel("Z[mm]")
ylabel("X[mm]")
axis equal tight
colorbar

figure()
subplot(1,2,1)
imagesc(aoVars(1).measVars.algo.len.zVecUSRes*1e3, vars.s2D.grid.firstVec, db(squeeze(dvRes.phi(:,1,:))))
title("2s")
xlabel("Z[mm]")
ylabel("X[mm]")
axis equal tight
colorbar
subplot(1,2,2)
imagesc(aoVars(1).measVars.algo.len.zVecUSRes*1e3, vars.s2D.grid.firstVec, db(squeeze(dvRes.phi(:,2,:))))
title("200s")
xlabel("Z[mm]")
ylabel("X[mm]")
axis equal tight
colorbar
%%
im2   = db(sqrt(squeeze(dvRes.phi(:,1,:))));
im200 = db(sqrt(squeeze(dvRes.phi(:,2,:))));

lpf = (1/9)*[ 1 1 1 ; ...
              1 1 1 ; ...
              1 1 1];
    
filtIm2   = conv2(im2, lpf, 'same');
filtIm200 = conv2(im200, lpf, 'same');

dx = abs(s2DVars(1).grid.firstVec(1) - s2DVars(1).grid.firstVec(2));
dz = abs(aoVars(1).measVars.algo.len.zVecUSRes(1) -  aoVars(1).measVars.algo.len.zVecUSRes(2));

[gx2, gz2] = gradient(filtIm2);
gx2Norm = gx2*dx;
gz2Norm = gz2*dz;

gTot2 = sqrt(gx2Norm.^2 + gz2Norm.^2);

[gx200, gz200] = gradient(filtIm200);
gx2Norm = gx200*dx;
gz2Norm = gz200*dz;

gTot200 = sqrt(gx2Norm.^2 + gz2Norm.^2);

figure()
subplot(1,2,1)
imagesc(gTot2, [0,2])
colorbar
title("2s")
subplot(1,2,2)
imagesc(gTot200, [0,2])
colorbar
title("200s")

figure()
subplot(1,2,1)
imagesc(aoVars(1).measVars.algo.len.zVecUSRes*1e3, vars.s2D.grid.firstVec, gTot2, [0,2])
title("2s")
xlabel("Z[mm]")
ylabel("X[mm]")
axis equal tight
colorbar
subplot(1,2,2)
imagesc(aoVars(1).measVars.algo.len.zVecUSRes*1e3, vars.s2D.grid.firstVec, gTot200, [0,2])
title("200s")
xlabel("Z[mm]")
ylabel("X[mm]")
axis equal tight
colorbar

%%
algo = Algo();

%% Check if fitting before channel mean has any signnificant impact
idx = 20;
for i=1:2
    aoVars(1).extVars.algo.envHar = 200e3;
    aoVars(1).extVars.algo.envUS = 78e3;
    aoVars(i).extVars.algo.export.netSignal      = false;
    aoVars(i).extVars.algo.export.deMultiplexed  = false;
    aoVars(i).extVars.algo.export.reshapedSignal = false;
    aoVars(i).extVars.algo.export.fftRes         = false;
    aoVars(i).extVars.algo.export.usCompCmplx    = false;
end

clear resNew
% resNew(vars.time.numOfTimeFrames,vars.s2D.grid.firstIdxLen) = struct('qAvgChFFT', [], 'unFittedFFT', [], 'unFittedFFTShift', [], 'fitModel', [], 'fittedFFT', [], 'phi', []);
resNew(vars.time.numOfTimeFrames,vars.s2D.grid.firstIdxLen) = struct('qAvgChFFT', [], 'fitModel', [], 'fittedFFT', [], 'fittedMeanFFT', [], 'phi', []);
dvResNew = zeros(vars.time.numOfTimeFrames, vars.s2D.grid.firstIdxLen);
for i =1:vars.time.numOfTimeFrames
    algo.updateAlgoUserVars(aoVars(i).extVars.algo);
    fprintf("analyzing T=%d\n", vars.time.timeVec(i));
    for j = 1:vars.s2D.grid.firstIdxLen
        resNew(i,j) = algo.analyseLoadedData(loadRes(i,j), 1);
        dvResNew(i,j) = resNew(i,j).phi(idx);
    end
end

res2NormNew = dvResNew(1,:) - min(dvResNew(1,:));
res2NormNew = res2NormNew / max(res2NormNew);

res200NormNew = dvResNew(2,:) - min(dvResNew(2,:));
res200NormNew = res200NormNew / max(res200NormNew);

figure();
plot(vars.s2D.grid.firstVec, dvResNew(1,:)); hold on;
plot(vars.s2D.grid.firstVec, dvResNew(2,:));
legend("T=2", "T=200")
title("fitting before channel mean")

figure();
plot(vars.s2D.grid.firstVec, db(res2NormNew)); hold on;
plot(vars.s2D.grid.firstVec, db(res200NormNew));
legend("T=2", "T=200")
title("fitting before channel mean (dB)")

figure();
plot(algo.freq.frequencyBarShifted, resNew(1,2).unFittedFFTShift(:, idx));hold on
plot(algo.freq.frequencyBarShifted(algo.freq.fitIdxShift), 0.5e-3*ones(1,length(algo.freq.fitIdxShift)), '+')
plot(algo.freq.frequencyBarShifted, resNew(2,2).fitModel(:,20));

% Answer is no :( 

%% calculate SNR

SNR2 = dvRes.phi(2,1,idx)/std(dvRes.phi(2,1,37:end));
SNR200 = dvRes.phi(2,2,idx)/std(dvRes.phi(2,2,37:end));

factorSNR = SNR200/SNR2;

%%
uVarsAO = aoVars(1).uVars;
uVarsAO.ao.envHar = 200e3;
uVarsAO.ao.envUS  = 78e3;
uVarsAO.ao.skipParamsCheck = false;

uVarsAO.fileSystem.extProject  = false;
uVarsAO.fileSystem.saveResults = false;
uVarsAO.fileSystem.saveFigs    = false;

uVarsAO.figs.intExt = 'int';
uVarsAO.figs.validStruct.phi = true;
uVarsAO.figs.ch    = 1;
uVarsAO.figs.zIdx  = 20;
uVarsAO.figs.quant = 1;

uVarsAO.figs.reopenFigures = true;
uVarsAO.figs.plotPhiInd = true; 
uVarsAO.figs.displayFullFFT = true;
uVarsAO.figs.FFTenv = 250e3; %Hz

uVarsAO.ao.timeToSample = 200;

ao.setMeasVars(uVarsAO);
ao.configPeripherals();
ao.runAcoustoOptics();
for j=1:5
    resAO{j} = ao.runAcoustoOptics();
end

%%
for j=1:5
    for i=1:ao.measVars.AO.splitNum
        phiMat(i,j,:) = resAO{j}.splitedRes(i).phi(:);
    end
    phiVec(j) =resAO{j}.phi(idx);
end

figure()
for i=1:5
    plot(phiMat(:,i,idx)); hold on
end

figure()
plot(phiVec)

%TODO: perform a consistency test