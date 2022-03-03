close all
clear all
clc;

%%
meshName = "1LayerSlabMesh-60x40-DxW-mm";
resDirName = "E:\AcoustoOptics\Simulations\Results\60x40-DxW-mm";
fullResFilename = sprintf("%s/%s - FullResults.mat", resDirName, meshName);

resTot = load(fullResFilename);
resTot = resTot.resTot;

muEffVec   = resTot.muEffVec;
muEffPivot = resTot.muEffPivot;
depthVec   = resTot.depthVec;
depthVecHR = resTot.depthVecHR;

slopeIdx = floor(length(depthVecHR)/2);
dl = abs(depthVecHR(2) - depthVecHR(1));
numOfmuEff = length(muEffVec);

for i=1:numOfmuEff
    phiHRMidUniMat(i,:)    = resTot.resArr{i}.phiHRMidUni;    %1D
    phiHRMidUniLogMat(i,:) = resTot.resArr{i}.phiHRMidUniLog; %1D
    phiHRMid2InMat(i,:)    = resTot.resArr{i}.phiHRMid2In;    %1D
    phiHRMid2InLogMat(i,:) = resTot.resArr{i}.phiHRMid2InLog; %1D
    phiHRMid1InMat(i,:)    = resTot.resArr{i}.phiHRMid1In;    %1D
    phiHRMid1InLogMat(i,:) = resTot.resArr{i}.phiHRMid1InLog; %1D
end

%%  
phiHRMidUniMatNorm    = normMatf(phiHRMidUniMat);
phiHRMidUniMatNormLog = log(phiHRMidUniMatNorm);
phiHRMid2InMatNorm    = normMatf(phiHRMid2InMat);
phiHRMid2InMatNormLog = log(phiHRMid2InMatNorm);
phiHRMid1InMatNorm    = normMatf(phiHRMid1InMat);
phiHRMid1InMatNormLog = log(phiHRMid1InMatNorm);

gValUni = abs(mean(diff(phiHRMidUniMatNormLog(:,1:slopeIdx),1,2),2))/dl;
gVal2In = abs(mean(diff(phiHRMid2InMatNormLog(:,1:slopeIdx),1,2),2))/dl;
gVal1In = abs(mean(diff(phiHRMid1InMatNormLog(:,1:slopeIdx),1,2),2))/dl;

gValUni(i) = abs(mean(gVecUni(1:slopeIdx)))/dl;
[gVec2In] = gradient(phiHRMid2InMatNormLog(i,:));
gVal2In(i) = abs(mean(gVec2In(1:slopeIdx)))/dl;
[gVec1In] = gradient(phiHRMid1InMatNormLog(i,:));
gVal1In(i) = abs(mean(gVec1In(1:slopeIdx)))/dl;

for i=1:numOfmuEff
    [gVecUni] = gradient(phiHRMidUniMatNormLog(i,:));
    gValUni(i) = abs(mean(gVecUni(1:slopeIdx)))/dl;
    [gVec2In] = gradient(phiHRMid2InMatNormLog(i,:));
    gVal2In(i) = abs(mean(gVec2In(1:slopeIdx)))/dl;
    [gVec1In] = gradient(phiHRMid1InMatNormLog(i,:));
    gVal1In(i) = abs(mean(gVec1In(1:slopeIdx)))/dl;
end

gValUniPivot = gValUni(3);
gVal2InPivot = gVal2In(3);
gVal1InPivot = gVal1In(3);

% Fluence in each illumination size for different phantoms
figure();
subplot(2,3,1)
for i=1:numOfmuEff
    plot(phiHRMidUniMatNorm(i,:)); hold on
    legStr(i) = sprintf("mu_{eff}=%.2f [cm^{-1}]", muEffVec(i));
end
legend(legStr)
title("Full Illumination Over Simulated Area")
subplot(2,3,2)
for i=1:numOfmuEff
    plot(phiHRMid2InMatNorm(i,:)); hold on
    legStr(i) = sprintf("mu_{eff}=%.2f [cm^{-1}]", muEffVec(i));
end
legend(legStr)
title("2in Illumination Over Simulated Area")
subplot(2,3,3)
for i=1:numOfmuEff
    plot(phiHRMid1InMatNorm(i,:)); hold on
    legStr(i) = sprintf("mu_{eff}=%.2f [cm^{-1}]", muEffVec(i));
end
legend(legStr)
title("1in Illumination Over Simulated Area")
subplot(2,3,4)
for i=1:numOfmuEff
    plot(phiHRMidUniMatNormLog(i,:)); hold on
    legStr(i) = sprintf("mu_{eff}=%.2f [cm^{-1}]", muEffVec(i));
end
legend(legStr, 'Location', 'southwest')
title("Full Illumination Over Simulated Area (Log)")
subplot(2,3,5)
for i=1:numOfmuEff
    plot(phiHRMid2InMatNormLog(i,:)); hold on
    legStr(i) = sprintf("mu_{eff}=%.2f [cm^{-1}]", muEffVec(i));
end
legend(legStr, 'Location', 'southwest')
title("2in Illumination Over Simulated Area (Log)")
subplot(2,3,6)
for i=1:numOfmuEff
    plot(phiHRMid1InMatNormLog(i,:)); hold on
    legStr(i) = sprintf("mu_{eff}=%.2f [cm^{-1}]", muEffVec(i));
end
legend(legStr, 'Location', 'southwest')
title("1in Illumination Over Simulated Area (Log)")

% Simulated muEff
figure();
subplot(1,2,1)
plot(muEffVec*10, muEffVec*10, '-+'); hold on
plot(muEffVec*10,gValUni*10, '-+'); hold on
plot(muEffVec*10,gVal2In*10, '-+');
plot(muEffVec*10,gVal1In*10, '-+');
title("Calculated \mu_{eff} vs. Simulated \mu_{eff}");
legend("Reference", "Uni", "2In", "1In");
ylabel("Simulated \mu_{eff} [1/cm]");
xlabel("Calculated \mu_{eff} [1/cm]");
subplot(1,2,2)
plot(muEffVec./muEffPivot, muEffVec./muEffPivot, '-+'); hold on
plot(muEffVec./muEffPivot,gValUni./gValUniPivot, '-+'); hold on
plot(muEffVec./muEffPivot,gVal2In./gVal2InPivot, '-+');
plot(muEffVec./muEffPivot,gVal1In./gVal1InPivot, '-+');
title(["Calculated \mu_{eff} vs. Simulated \mu_{eff}"; "Normalized by Pivot Value"]);
legend("Reference", "Uni", "2In", "1In");
ylabel("Simulated \mu_{eff} [1/cm]");
xlabel("Calculated \mu_{eff} [1/cm]");

%%
close all;
clear all;
clc

meshWidth = [30, 40, 50];
meshDepth = [5, 10, 15, 20, 30, 40, 50, 60];

numOfWidth = length(meshWidth);
numOfDepth = length(meshDepth);

for k = 1:length(meshWidth)
    for j = 1:length(meshDepth)
        meshName = sprintf("1LayerSlabMesh-%dx%d-DxW-mm",  meshDepth(j), meshWidth(k));
        resDirName = sprintf("E:/AcoustoOptics/Simulations/Results/%dx%d-DxW-mm", meshDepth(j), meshWidth(k));
        fullResFilename = sprintf("%s/%s - FullResults.mat", resDirName, meshName);

        resTot(k,j) = load(fullResFilename);
%         resTot(k,j,i) = resTmp.resTot;

        muEffVec   = resTot(k,j).muEffVec;
        muEffPivot = resTot(k,j).muEffPivot;
        depthVec{k,j}   = resTot(k,j).depthVec;
        depthVecHR{k,j} = resTot(k,j).depthVecHR;

        slopeIdx(k,j) = floor(length(depthVecHR{k,j})/2);
        dl(k,j) = abs(depthVecHR{k,j}(2) - depthVecHR{k,j}(1));
        numOfmuEff = length(muEffVec);

        for i=1:numOfmuEff
            phiUni{k,j}(i,:)    = resTot(k,j).resArr{i}.phiHRMidUni;    %1D
            phiUniLog{k,j}(i,:) = resTot(k,j).resArr{i}.phiHRMidUniLog; %1D
            phi2In{k,j}(i,:)    = resTot(k,j).resArr{i}.phiHRMid2In;    %1D
            phi2InLog{k,j}(i,:) = resTot(k,j).resArr{i}.phiHRMid2InLog; %1D
            phi1In{k,j}(i,:)    = resTot(k,j).resArr{i}.phiHRMid1In;    %1D
            phi1InLog{k,j}(i,:) = resTot(k,j).resArr{i}.phiHRMid1InLog; %1D
        end

        phiUniNorm{k,j} = normMatf(phiUni{k,j},2);
        phiUniLog{k,j}  = log(phiUniNorm{k,j});
        phi2InNorm{k,j} = normMatf(phi2In{k,j},2);
        phi2InLog{k,j}  = log(phi2InNorm{k,j});
        phi1InNorm{k,j} = normMatf(phi1In{k,j},2);
        phi1InLog{k,j}  = log(phi1InNorm{k,j});

        gValUni(k,j,:) = abs(mean(diff(phiUniLog{k,j}(:,1:slopeIdx(k,j)),1,2),2))/dl(k,j);
        gVal2In(k,j,:) = abs(mean(diff(phi2InLog{k,j}(:,1:slopeIdx(k,j)),1,2),2))/dl(k,j);
        gVal1In(k,j,:) = abs(mean(diff(phi1InLog{k,j}(:,1:slopeIdx(k,j)),1,2),2))/dl(k,j);
    
    end
end

pivotIdx = 3;

%% Phi vs X for a given MuEFF with different phantoms width
figure();
subplot(1,3,1)
%(1) in a given depth for different width 
clear legStr
i = 3; j = 4;
plot(depthVec{k,j}, log(exp(-muEffVec(i)*depthVec{k,j}))); hold on
legStr(1) = "Reference";
for k = 1:numOfWidth
    plot(depthVec{k,j}, squeeze(phiUniLog{k,j}(i,:)) , '-+'); hold on
    legStr(k) = sprintf("%d cm", meshWidth(k));
end
legend(legStr, 'Location', 'northeast')
ylabel("Simulated \mu_{eff} [1/cm]");
xlabel("Calculated \mu_{eff} [1/cm]");
title(["Full Illumination(Log)"; sprintf("Phantom Depth: %d [cm]", meshDepth(j))])
subplot(1,3,2)
%(1) in a given depth for different width 
clear legStr
i = 3; j = 4;
plot(depthVec{k,j}, log(exp(-muEffVec(i)*depthVec{k,j}))); hold on
legStr(1) = "Reference";
for k = 1:numOfWidth
    plot(depthVec{k,j}, squeeze(phi2InLog{k,j}(i,:)) , '-+'); hold on
    legStr(k) = sprintf("%d cm", meshWidth(k));
end
legend(legStr, 'Location', 'northeast')
ylabel("Simulated \mu_{eff} [1/cm]");
xlabel("Calculated \mu_{eff} [1/cm]");
title(["2In Illumination(Log)"; sprintf("Phantom Depth: %d [cm]", meshDepth(j))])
subplot(1,3,3)
%(1) in a given depth for different width 
clear legStr
i = 3; j = 4;
plot(depthVec{k,j}, log(exp(-muEffVec(i)*depthVec{k,j}))); hold on
legStr(1) = "Reference";
for k = 1:numOfWidth
    plot(depthVec{k,j}, squeeze(phi1InLog{k,j}(i,:)) , '-+'); hold on
    legStr(k) = sprintf("%d cm", meshWidth(k));
end
legend(legStr, 'Location', 'northeast')
ylabel("Simulated \mu_{eff} [1/cm]");
xlabel("Calculated \mu_{eff} [1/cm]");
title(["1In Illumination(Log)"; sprintf("Phantom Depth: %d [cm]", meshDepth(j))])

%% Phi vs X for a given MuEFF with different phantoms depth
figure()
subplot(1,3,1)
clear legStr
i = 3; k = 2;
plot(depthVec{k,end}, log(exp(-muEffVec(i)*depthVec{k,end}))); hold on
legStr(1) = "Reference";
for j = 1:numOfDepth
    plot(depthVec{k,j}, squeeze(phiUniLog{k,j}(i,:)) , '-+'); hold on
    legStr(j+1) = sprintf("%d cm", meshDepth(j));
end
legend(legStr, 'Location', 'northeast')
ylabel("Simulated \mu_{eff} [1/cm]");
xlabel("Calculated \mu_{eff} [1/cm]");
title(["Full Illumination(Log)"; sprintf("Phantom Width: %d [cm]", meshWidth(k))])
subplot(1,3,2)
%(2) in given width for different depth
clear legStr
i = 3; k = 2;
plot(depthVec{k,end}, log(exp(-muEffVec(i)*depthVec{k,end}))); hold on
legStr(1) = "Reference";
for j = 1:numOfDepth
    plot(depthVec{k,j}, squeeze(phi2InLog{k,j}(i,:)) , '-+'); hold on
    legStr(j+1) = sprintf("%d cm", meshDepth(j));
end
legend(legStr, 'Location', 'northeast')
ylabel("Simulated \mu_{eff} [1/cm]");
xlabel("Calculated \mu_{eff} [1/cm]");
title(["2In Illumination(Log)"; sprintf("Phantom Width: %d [cm]", meshWidth(k))])
subplot(1,3,3)
%(2) in given width for different depth
clear legStr
i = 3; k = 2;
plot(depthVec{k,end}, log(exp(-muEffVec(i)*depthVec{k,end}))); hold on
legStr(1) = "Reference";
for j = 1:numOfDepth
    plot(depthVec{k,j}, squeeze(phi1InLog{k,j}(i,:)) , '-+'); hold on
    legStr(j+1) = sprintf("%d cm", meshDepth(j));
end
legend(legStr, 'Location', 'northeast')
ylabel("Simulated \mu_{eff} [1/cm]");
xlabel("Calculated \mu_{eff} [1/cm]");
title(["1In Illumination(Log)"; sprintf("Phantom Width: %d [cm]", meshWidth(k))])


%% Phantom Width Effect
figure();
subplot(1,3,1)
legStr(1) = "Reference";
plot(muEffVec./muEffPivot, muEffVec./muEffPivot, '-+'); hold on
i = 3; j = 2;
for k = 1:numOfWidth
    plot(muEffVec./muEffPivot, squeeze(gValUni(k,j,:)/gValUni(k,j,pivotIdx)) , '-+');
    legStr(k+1) = sprintf("%d cm", meshWidth(k));
end
legend(legStr, 'Location', 'northwest')
xlabel("\mu_{eff}^{calc} / \mu_{eff}^{calc- pivot}");
ylabel("\mu_{eff}^{sim} / \mu_{eff}^{sim- pivot}");
title(["Full Illumination(Log)"; sprintf("Phantom Depth: %d [cm]", meshDepth(j))])
subplot(1,3,2)
legStr(1) = "Reference";
plot(muEffVec./muEffPivot, muEffVec./muEffPivot, '-+'); hold on
i = 3; j = 4;
for k = 1:numOfWidth
    plot(muEffVec./muEffPivot, squeeze(gValUni(k,j,:)/gValUni(k,j,pivotIdx)) , '-+');
    legStr(k+1) = sprintf("%d cm", meshWidth(k));
end
legend(legStr, 'Location', 'northwest')
xlabel("\mu_{eff}^{calc} / \mu_{eff}^{calc- pivot}");
ylabel("\mu_{eff}^{sim} / \mu_{eff}^{sim- pivot}");
title(["Full Illumination(Log)"; sprintf("Phantom Depth: %d [cm]", meshDepth(j))])
subplot(1,3,3)
legStr(1) = "Reference";
plot(muEffVec./muEffPivot, muEffVec./muEffPivot, '-+'); hold on
i = 3; j = 6;
for k = 1:numOfWidth
    plot(muEffVec./muEffPivot, squeeze(gValUni(k,j,:)/gValUni(k,j,pivotIdx)) , '-+');
    legStr(k+1) = sprintf("%d cm", meshWidth(k));
end
legend(legStr, 'Location', 'northwest')
xlabel("\mu_{eff}^{calc} / \mu_{eff}^{calc- pivot}");
ylabel("\mu_{eff}^{sim} / \mu_{eff}^{sim- pivot}");
title(["Full Illumination(Log)"; sprintf("Phantom Depth: %d [cm]", meshDepth(j))])

%% Phantom Depth Effect
figure();
subplot(1,3,1)
legStr(1) = "Reference";
plot(muEffVec./muEffPivot, muEffVec./muEffPivot, '-+'); hold on
i = 3; k = 1;
for j = 1:numOfDepth
    plot(muEffVec./muEffPivot, squeeze(gValUni(k,j,:)/gValUni(k,j,pivotIdx)) , '-+');
    legStr(j+1) = sprintf("%d cm", meshDepth(j));
end
legend(legStr, 'Location', 'northwest')
xlabel("\mu_{eff}^{calc} / \mu_{eff}^{calc- pivot}");
ylabel("\mu_{eff}^{sim} / \mu_{eff}^{sim- pivot}");
title(["Full Illumination(Log)"; sprintf("Phantom Width: %d [cm]", meshWidth(k))])
subplot(1,3,2)
legStr(1) = "Reference";
plot(muEffVec./muEffPivot, muEffVec./muEffPivot, '-+'); hold on
i = 3; k = 2;
for j = 1:numOfDepth
    plot(muEffVec./muEffPivot, squeeze(gValUni(k,j,:)/gValUni(k,j,pivotIdx)) , '-+');
    legStr(j+1) = sprintf("%d cm", meshDepth(j));
end
legend(legStr, 'Location', 'northwest')
xlabel("\mu_{eff}^{calc} / \mu_{eff}^{calc- pivot}");
ylabel("\mu_{eff}^{sim} / \mu_{eff}^{sim- pivot}");
title(["Full Illumination(Log)"; sprintf("Phantom Width: %d [cm]", meshWidth(k))])
subplot(1,3,3)
legStr(1) = "Reference";
plot(muEffVec./muEffPivot, muEffVec./muEffPivot, '-+'); hold on
i = 3; k = 3;
for j = 1:numOfDepth-2
    plot(muEffVec./muEffPivot, squeeze(gValUni(k,j,:)/gValUni(k,j,pivotIdx)) , '-+');
    legStr(j+1) = sprintf("%d cm", meshDepth(j));
end
legend(legStr, 'Location', 'northwest')
xlabel("\mu_{eff}^{calc} / \mu_{eff}^{calc- pivot}");
ylabel("\mu_{eff}^{sim} / \mu_{eff}^{sim- pivot}");
title(["Full Illumination(Log)"; sprintf("Phantom Width: %d [cm]", meshWidth(k))])




%% How the length of the slope upon gradient is calculated affecting muEff recon
k = 2; i=3; j = 8;
slopeIdxVec = 1:120;
refVec = muEffVec(i)*ones(1,length(slopeIdxVec));
dVec = depthVec{k,end}(slopeIdxVec);
for j = 1:numOfDepth
    numOfSlope(j) = sum(slopeIdxVec < length(depthVec{k,j}));
    for m = 1:numOfSlope(j)
        gValUniSlopeIdx{j}(m) = abs(mean(diff(phiUniLog{k,j}(i,3:slopeIdxVec(m)),1,2),2))/dl(k,j);
        gVal2InSlopeIdx{j}(m) = abs(mean(diff(phi2InLog{k,j}(i,3:slopeIdxVec(m)),1,2),2))/dl(k,j);
        gVal1InSlopeIdx{j}(m) = abs(mean(diff(phi1InLog{k,j}(i,3:slopeIdxVec(m)),1,2),2))/dl(k,j);
    end
end

figure();
subplot(1,3,1)
plot(dVec(1:numOfSlope(end)),  refVec); hold on
legStr(1) = sprintf("Reference");
for j = 1:numOfDepth
    plot(dVec(1:numOfSlope(j)), gValUniSlopeIdx{j});
    legStr(j+1) = sprintf("%d cm", meshDepth(j));
end
legend(legStr);
xlabel( "Length of Slope [mm]")
ylabel(" Calculated \mu_{eff}")
title("Full Illumination")
ylim([0.15, 0.35])
subplot(1,3,2)
plot(dVec(1:numOfSlope(end)),  refVec); hold on
legStr(1) = sprintf("Reference");
for j = 1:numOfDepth
    plot(dVec(1:numOfSlope(j)), gVal2InSlopeIdx{j});
    legStr(j+1) = sprintf("%d cm", meshDepth(j));
end
legend(legStr);
xlabel( "Length of Slope [mm]")
ylabel(" Calculated \mu_{eff}")
title("2In Illumination")
ylim([0.15, 0.35])
subplot(1,3,3)
plot(dVec(1:numOfSlope(end)),  refVec); hold on
legStr(1) = sprintf("Reference");
for j = 1:numOfDepth
    plot(dVec(1:numOfSlope(j)), gVal1InSlopeIdx{j});
    legStr(j+1) = sprintf("%d cm", meshDepth(j));
end
legend(legStr);
ylim([0.15, 0.35])
xlabel( "Length of Slope [mm]")
ylabel(" Calculated \mu_{eff}")
title("1In Illumination")
