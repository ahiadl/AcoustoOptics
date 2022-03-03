res = load('D:\ResultsToKeep\FellgettsAdvantage\AcoustoOptics-1cyc-woHad-43n-20mW\1-Results.mat');
vars = load('D:\ResultsToKeep\FellgettsAdvantage\AcoustoOptics-1cyc-woHad-43n-20mW\Vars.mat');

rawData= res.res.rawData;
uVarsAlgo = vars.vars.algo.uVars;
uVarsAlgo.export.deMultiplexed  = true;
uVarsAlgo.export.fftRes         = true;
uVarsAlgo.export.reshapedSignal = true;
uVarsAlgo.export.usCompCmplx    = true;
uVarsAlgo.export.rawData        = true;
uVarsAlgo.export.netSignal      = true;
uVarsAlgo.useGPU = false;
algo = Algo();

algoOld = AlgoOld();
% algoTmp = AlgoTmp();
% algoTmp2 = AlgoOld2();
%% New
uVarsAlgo.useHadamard           = false;
uVarsAlgo.envDC                 = 100e3;
uVarsAlgo.envUS                 = 78e3;
uVarsAlgo.envHar                = 200e3;
uVarsAlgo.longMeas              = false;

algo.updateAlgoUserVars(uVarsAlgo);
algo.setRawData(rawData);
resNew = algo.analyse();

phiNew     = resNew.phi;
maxValNew  = max(phiNew);
minValNew  = min(phiNew);
spanPhiNew = maxValNew - minValNew;
phiNewNorm = (phiNew - minValNew)/spanPhiNew;

SNRNew = 1/std(phiNewNorm(34:end));

figure(); stem(resNew.rawPhi); title("new")
figure(); stem(res.res.phi); title("old")




%%



%% Old
algoOld.updateAlgoUserVars(uVarsAlgo);
algoOld.setRawData(res.rawData);
resOld = algoOld.analyse();

phiOld     = resOld.phi;
maxValOld  = max(phiOld);
minValOld  = min(phiOld);
spanPhiOld = maxValOld - minValOld;
phiOldNorm = (phiOld - minValOld)/spanPhiOld;

SNROld = 1/std(phiOldNorm(34:end));

%% Temporary
% algoTmp.updateAlgoUserVars(vars.measVars.algo.uVars);
% algoTmp.setRawData(res.rawData);
% resTmp = algoTmp.analyse();
% 
% phiTmp     = resTmp.phi;
% maxValTmp  = max(phiTmp);
% minValTmp  = min(phiTmp);
% spanPhiTmp = maxValTmp - minValTmp;
% phiTmpNorm = (phiTmp - minValTmp)/spanPhiTmp;
% 
% SNRTmp = 1/std(phiTmpNorm(34:end));

%% Temporary 2
% algoTmp2.updateAlgoUserVars(vars.measVars.algo.uVars);
% algoTmp2.setRawData(res.rawData);
% resTmp2 = algoTmp2.analyse();
% 
% phiTmp2     = resTmp2.phi;
% maxValTmp2  = max(phiTmp2);
% minValTmp2  = min(phiTmp2);
% spanPhiTmp2 = maxValTmp2 - minValTmp2;
% phiTmp2Norm = (phiTmp2 - minValTmp2)/spanPhiTmp2;
% 
% SNRTmp2 = 1/std(phiTmp2Norm(34:end));

% %% Differences Analysis
% 
% % Demultiplexing is OK
% diffMat = abs((resOld.deMultiplexed - resNew.deMultiplexed));
% maxDiff = max(max(diffMat));
% diffVal = sum(sum(abs(diffMat)));
% 
% % Reshape is OK
% diffMat2 = abs(resOld.reshapedSignal - resNew.reshapedSignal);
% maxDiff2 = max(max(max(max(diffMat2))));
% diffVal2 = sum(sum(sum(sum(abs(diffMat2)))));
% 
% %
% diffMat3 = abs(abs(resOld.phiChCmplx) - abs(resNew.usCompCmplx));
% maxDiff3 = max(max(max(diffMat3)));
% diffVal3 = sum(sum(sum(diffMat3)));


%%
figure()
subplot(1,3,1)
stem(phiNewNorm)
title("New")
subplot(1,3,2)
stem(phiOldNorm)
title("Old")
% subplot(1,3,3)
% stem(phiTmp2Norm)
% title("Tmp")


%%
% rawPhi = squeeze(s2D.res2D.rawPhi);
maxValNew  = max(max(rawPhi));
minValNew  = min(min(rawPhi));
spanPhiNew = maxValNew - minValNew;
normS2DResRaw = (rawPhi - minValNew)/spanPhiNew;

figure(); 
subplot(2,1,1); 
imagesc(rawPhi); 
axis equal tight
subplot(2,1,2); 
imagesc(log(normS2DResRaw))
axis equal tight


rawPhi2 = squeeze(s2D.res2D.rawPhi);
maxValNew  = max(max(rawPhi2));
minValNew  = min(min(rawPhi2));
spanPhiNew = maxValNew - minValNew;
normS2DResRaw2 = (rawPhi2 - minValNew)/spanPhiNew;

figure(); 
subplot(2,1,1); 
imagesc(rawPhi2); 
axis equal tight
subplot(2,1,2); 
imagesc(log(normS2DResRaw2))
axis equal tight




%%
res = load('D:\ResultsToKeep\Calibration\Scan1D-wHad-100mW-2Fibers\2DScan-Results.mat');
vars = load('D:\ResultsToKeep\Calibration\Scan1D-wHad-100mW-2Fibers\2DScan-Vars.mat');

phi = squeeze(res.phi);
phiNorm = phi+min(min(phi));

figure();
subplot(2,1,1)
imagesc(phi)
title("2Fibers")
axis tight equal
subplot(2,1,2)
imagesc(db(phi))
axis tight equal


res = load('D:\ResultsToKeep\Calibration\Scan1D-wHad-100mW-Bundle\2DScan-Results.mat');
vars = load('D:\ResultsToKeep\Calibration\Scan1D-wHad-100mW-Bundle\2DScan-Vars.mat');

phi = squeeze(res.phi);
phiNorm = phi+min(min(phi));

figure();
subplot(2,1,1)
imagesc(phi)
title("1 fiber")
axis tight equal
subplot(2,1,2)
imagesc(db(phi))
axis tight equal










