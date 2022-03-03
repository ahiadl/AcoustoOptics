close all
clear all
clc

addpath(genpath(pwd))

%% Mesh 
toastThreadCount(20);

%Create the mesh from the geometrical points defined in the .geo file
fprintf("Constructing Mesh\n");
if ~exist('./Simulations/Toast++Sims/1LayerSlabMesh2cm-Wide.msh', 'file')
    fprintf("Can't locate the mesh. Generating...\n"); 
    system ('gmsh -3 ./Simulations/Toast++Sims/1LayerBox.geo -o ./Simulations/Toast++Sims/1LayerSlabMesh2cm-Wide.msh');
else
    fprintf("Located the mesh. Loading without generating.\n");
end

% Load the mesh into toastMesh object.
% Notice: If the .msh file was already constructed the above code is
% redundant and the mesh should only be loaded with the code below.
fprintf("Loading Mesh\n");
tic;
mesh = toastMesh('./Simulations/Toast++Sims/1LayerSlabMesh2cm-Wide.msh','gmsh');
toc

dirName = "D:\Results\MuEff\Simulations\12-01-22-TopHatFirstTry\Simulations\20mmDepthWide";
if ~exist(dirName, 'dir')
    mkdir(dirName)
end
%% Define Sizes of Mesh Axes:
dl      = 0.5; % distance between 2 elements in the mesh in each one of the axes.

bbMesh = mesh.BoundingBox;

%Define axes lengths:
xLen    = bbMesh(2,1) - bbMesh(1,1);  % length of x axis of the mesh -[mm]
yLen    = bbMesh(2,2) - bbMesh(1,2);  % length of y axis of the mesh -[mm] 
zLen    = bbMesh(2,3) - bbMesh(1,3);  % length of z axis of the mesh -[mm]

% Direction with relation to the illumination:
% Lateral axis      - X [0, len]
% Transversal plane - YZ [-len/2, len/2]
% NOTICE: Illumination (planar) will be applied on the YZ and wil spread into the
% phantom along x axis

% Calculate middle lengths:
xMidLen = xLen/2; 
yMidLen = yLen/2;
zMidLen = zLen/2;

% Create axes vectors:
xVec = 0:dl:xLen;
yVec = -yMidLen:dl:yMidLen;
zVec = -zMidLen:dl:zMidLen;

% Calculate Indexes:
xSize = length(xVec);
ySize = length(yVec);
zSize = length(zVec);

% Index of the axis center:
xMid = floor(xSize/2); 
yMid = floor(ySize/2);
zMid = floor(zSize/2);

% Create mesh grid that will help cast the toastMesh onto.
[X, Y, Z] = meshgrid (yVec, xVec, zVec);

% Create rMat to describe distances of each point in the mesh from the
% illumination point on the boundary (relevant for point source)
rMat = sqrt(X.^2 + Y.^2 + Z.^2);

grd = [xSize, ySize, zSize];

%Create high resolution grid for uniform illumination simulation
dlHR = dl/10;
xVecHR = 0:dlHR:xLen;
yVecHR = -yLen/2:dlHR:yLen/2;
zVecHR = -zLen/2:dlHR:zLen/2;
% zVecHR = -zLen/2:dl:zLen/2;

[XHR, YHR, ZHR] = meshgrid (yVecHR, xVecHR, zVecHR);

xSizeHR = length(xVecHR);
ySizeHR = length(yVecHR);
zSizeHR = length(zVecHR);  

yMidHRIdx = floor(ySizeHR/2)+1;
zMidHRIdx = floor(zSizeHR/2)+1;

mask2In = sqrt(ZHR.^2 + XHR.^2) <  25.4;
mask1In = sqrt(ZHR.^2 + XHR.^2) <  12.7;

slopeIdx = 200;

%% Create Toast Basis
% Casting the coordinates defined in the previose section onto the mesh 
% generated with gMesh 

fprintf("Creating Basis\n");
tic
basis = toastBasis (mesh,grd);
toc

ne             = mesh.ElementCount;
nv             = mesh.NodeCount;
regidx         = mesh.Region;
regno          = unique(regidx);%sorted from small to largest.
firstLayerIdx  = find(regidx == regno(1));

%%
% ref  = 1.37;
% mus  = 0.7;  %[1/mm]
% mua  = 0.03; %[1/mm]

% muaVec   = [mua/4, mua/3, mua/2, mua/1.5, mua, mua*1.5, mua*2, mua*3, mua*4];
% musVec   = [mus/4, mus/3, mus/2, mus/1.5, mus, mus*1.5, mus*2, mus*3, mus*4];
% muEffVec = sqrt(3*muaVec.*musVec);

ref  = 1.37;
musP  = 0.4; %Prime(!) 1/mm
muaPivot  = 0.0184; % 1/mm  

muaVec   = [muaPivot/4, muaPivot/2, muaPivot, muaPivot*2, muaPivot*4];
musPVec   = musP*ones(1,length(muaVec));
muEffVec = sqrt(3*muaVec.*(muaVec+musPVec));

srcVec = [0,0,0];
mesh.SetQM(srcVec, []);
qvec = mesh.Qvec('Neumann','Point',0);
%%
for i =1 :length(muEffVec)
    fprintf("Starts a New Iteration\n");
    curMua  = muaVec (i);
    curMusP  = musPVec(i);
    
    startIdx  = 16;
    endIdx    = 46;

    muEffCalc     = sqrt(3*(curMua+curMusP)*curMua);

    ref                 = ones(ne,1);
    ref(firstLayerIdx)  = ref;

    musP                 = ones(ne,1);
    musP(firstLayerIdx)  = curMusP;

    mua                 = ones(ne,1);
    mua(firstLayerIdx)  = curMua;

    musPMat = zeros(xSize, ySize, zSize);
    musPMat( :, :, :) = curMusP;

    muaMat = zeros(xSize, ySize, zSize);
    muaMat( :, :, :) = curMua;

    % rhs (single source)
    fprintf("Creating SysMat\n");
    tic
    K = dotSysmat(mesh, mua ,musP, ref, 'EL');
    toc

    % Solve
    fprintf("Solving\n");
    tic
    phi = K\qvec;
    phiMapped  = basis.Map('M->B',phi);
    toc
    
    fprintf("Analyzing\n");
    % Extract Phi: Linear and Log
    phiGrd     = reshape(phiMapped, xSize, ySize, zSize);
    phiGrdLog  = log(abs(phiGrd));

    fprintf("Interpolating\n");
    tic
    phiHR       = interp3(X, Y, Z, phiGrd, XHR, YHR, ZHR);
    toc

    fprintf("Extracting\n");
    phiHRLog    = log(phiHR);
    phiHRMid    = phiHR(:, :, zMidHRIdx);
    phiHRMidLog = log(phiHRMid);
    
    phiHRMidUni    =(sum(sum(phiHR,3),2));
    phiHRMidUniNorm = normMatf(phiHRMidUniMat);
    phiHRMidUniLog = log(phiHRMidUniNorm);
    
    phiHRMid2In    = sum(sum(phiHR.*mask2In,3),2);
    phiHRMid2InNorm = normMatf(phiHRMid2In);
    phiHRMid2InLog = log(phiHRMid2InNorm);
    
    phiHRMid1In     = sum(sum(phiHR.*mask1In,3),2);
    phiHRMid1InNorm = normMatf(phiHRMid1In);
    phiHRMid1InLog  = log(phiHRMid1InNorm);

    [gVecUniLog] = gradient(phiHRMidUniLog);
    gValUni = abs(mean(gVecUniLog(1:slopeIdx)))/dlHR;
    [gVec2InLog] = gradient(phiHRMid2InLog);
    gVal2In = abs(mean(gVec2InLog(1:slopeIdx)))/dlHR;
    [gVec1InLog] = gradient(phiHRMid1InLog);
    gVal1In = abs(mean(gVec1InLog(1:slopeIdx)))/dlHR;
    
    res.phiGrd          = phiGrd;         %3D
    res.phiHR           = phiHR;          %3D
    res.phiHRLog        = phiHRLog;       %3D
    res.phiHRMid        = phiHRMid;       %2D
    res.phiHRMidLog     = phiHRMidLog;    %2D
    res.phiHRMidUni     = phiHRMidUni;    %1D
    res.phiHRMidUniNorm = phiHRMidUniNorm;%1D
    res.phiHRMidUniLog  = phiHRMidUniLog; %1D
    res.phiHRMid2In     = phiHRMid2In;    %1D
    res.phiHRMid2InNorm = phiHRMid2InNorm;%1D
    res.phiHRMid2InLog  = phiHRMid2InLog; %1D
    res.phiHRMid1In     = phiHRMid1In;    %1D
    res.phiHRMid1InNorm = phiHRMid1InNorm; %1D
    res.phiHRMid1InLog  = phiHRMid1InLog; %1D
    
    res.gVecUniLog        = gVecUniLog;
    res.gValUni           = gValUni;
    res.gVec2InLog        = gVec2InLog;
    res.gVal2In           = gVal2In;
    res.gVec1InLog        = gVec1InLog;
    res.gVal1In           = gVal1In;

%     figure();
%     subplot(2,1,1)
%     imagesc(yVecHR, xVecHR, phiHRMidLog)
%     axis tight equal
%     subplot(2,1,2)
%     plot(xVecHR, phiHRMidUniLog); hold on
%     plot(xVecHR, phiHRMid2InLog); 
%     plot(xVecHR, phiHRMid1InLog); 
%     legend("Infinite", "2in", "1in")

    filename = sprintf("%s/IncreasingMuEffSingleLayerUniformIllumination - %d.mat", dirName,i);
    
    fprintf("Saving\n");
    save(filename, 'res',  '-v7.3')
    
    resTot.resArr{i}.phiHRMidUni    = phiHRMidUni;    %1D
    resTot.resArr{i}.phiHRMidUniLog = phiHRMidUniLog; %1D
    resTot.resArr{i}.phiHRMid2In    = phiHRMid2In;    %1D
    resTot.resArr{i}.phiHRMid2InLog = phiHRMid2InLog; %1D
    resTot.resArr{i}.phiHRMid1In    = phiHRMid1In;    %1D
    resTot.resArr{i}.phiHRMid1InLog = phiHRMid1InLog; %1D
    resTot.resArr{i}.gVecLog        = gVecLog;
    resTot.resArr{i}.gVal           = gValUni;
    
    res = [];
    phiHR = [];
    phiHRLog = [];
end

resTot.ref       = ref;
resTot.muaVec    = muaVec; 
resTot.musVec    = musPVec;
resTot.muEffVec  = muEffVec;

resTot.depthVec   = xVec;
resTot.depthVecHR = xVecHR;

% collectResSingleLayer;
fprintf("Saving total results\n");
filename = sprintf("%s/IncreasingMuEffSingleLayerUniformIllumination - FullResuls.mat", dirName);
save(filename, 'resTot', '-v7.3')
fprintf("Done simulations\n");

%%

for i=1:length(muEffVec)
    phiHRMidUniMat(i,:)    = resTot.resArr{i}.phiHRMidUni;%1D
    phiHRMidUniLogMat(i,:) = resTot.resArr{i}.phiHRMidUniLog; %1D
    phiHRMid2InMat(i,:)    = resTot.resArr{i}.phiHRMid2In;    %1D
    phiHRMid2InLogMat(i,:) = resTot.resArr{i}.phiHRMid2InLog; %1D
    phiHRMid1InMat(i,:)    = resTot.resArr{i}.phiHRMid1In;    %1D
    phiHRMid1InLogMat(i,:) = resTot.resArr{i}.phiHRMid1InLog; %1D
end

%%

phiHRMidUniMatNorm = normMatf(phiHRMidUniMat);
phiHRMidUniMatNormLog = log(phiHRMidUniMatNorm);


phiHRMid2InMatNorm = normMatf(phiHRMid2InMat);
phiHRMid2InMatNormLog = log(phiHRMid2InMatNorm);


phiHRMid1InMatNorm = normMatf(phiHRMid1InMat);
phiHRMid1InMatNormLog = log(phiHRMid1InMatNorm);

for i=1:length(muEffVec)
    [gVecUni] = gradient(phiHRMidUniMatNormLog(i,:));
    gValUni(i) = abs(mean(gVecUni(1:200)))/dlHR;
    [gVec2In] = gradient(phiHRMid2InMatNormLog(i,:));
    gVal2In(i) = abs(mean(gVec2In(1:200)))/dlHR;
    [gVec1In] = gradient(phiHRMid1InMatNormLog(i,:));
    gVal1In(i) = abs(mean(gVec1In(1:200)))/dlHR;
end


% Fluence in each illumination size for different phantoms
figure();
subplot(2,3,1)
plot(phiHRMidUniMatNorm(1,:)); hold on
plot(phiHRMidUniMatNorm(2,:)); hold on
plot(phiHRMidUniMatNorm(3,:)); hold on
title("Full Illumination Over Simulated Area")
legend("Light", "Pivot", "Dark")
subplot(2,3,2)
plot(phiHRMid2InMatNorm(1,:)); hold on
plot(phiHRMid2InMatNorm(2,:)); hold on
plot(phiHRMid2InMatNorm(3,:)); hold on
title("2in Illumination Over Simulated Area")
legend("Light", "Pivot", "Dark")
subplot(2,3,3)
plot(phiHRMid1InMatNorm(1,:)); hold on
plot(phiHRMid1InMatNorm(2,:)); hold on
plot(phiHRMid1InMatNorm(3,:)); hold on
title("1in Illumination Over Simulated Area")
legend("Light", "Pivot", "Dark")
subplot(2,3,4)
plot(phiHRMidUniMatNormLog(1,:)); hold on
plot(phiHRMidUniMatNormLog(2,:)); hold on
plot(phiHRMidUniMatNormLog(3,:)); hold on
title("Full Illumination Over Simulated Area (Log)")
legend("Light", "Pivot", "Dark")
subplot(2,3,5)
plot(phiHRMid2InMatNormLog(1,:)); hold on
plot(phiHRMid2InMatNormLog(2,:)); hold on
plot(phiHRMid2InMatNormLog(3,:)); hold on
title("2in Illumination Over Simulated Area (Log)")
legend("Light", "Pivot", "Dark")
subplot(2,3,6)
plot(phiHRMid1InMatNormLog(1,:)); hold on
plot(phiHRMid1InMatNormLog(2,:)); hold on
plot(phiHRMid1InMatNormLog(3,:)); hold on
title("1in Illumination Over Simulated Area (Log)")
legend("Light", "Mid", "Dark")

% Fluence in Each phantom under different illumination size
figure();
subplot(2,3,1)
plot(phiHRMidUniMatNorm(1,:)); hold on
plot(phiHRMid2InMatNorm(1,:)); hold on
plot(phiHRMid1InMatNorm(1,:)); hold on
title("Light")
legend("Uni", "2In", "1In")
subplot(2,3,2)
plot(phiHRMidUniMatNorm(2,:)); hold on
plot(phiHRMid2InMatNorm(2,:)); hold on
plot(phiHRMid1InMatNorm(2,:)); hold on
title("Pivot")
legend("Uni", "2In", "1In")
subplot(2,3,3)
plot(phiHRMidUniMatNorm(3,:)); hold on
plot(phiHRMid2InMatNorm(3,:)); hold on
plot(phiHRMid1InMatNorm(3,:)); hold on
title("Dark")
legend("Uni", "2In", "1In")
subplot(2,3,4)
plot(phiHRMidUniMatNormLog(1,:)); hold on
plot(phiHRMid2InMatNormLog(1,:)); hold on
plot(phiHRMid1InMatNormLog(1,:)); hold on
title("Light (Log)")
legend("Uni", "2In", "1In")
subplot(2,3,5)
plot(phiHRMidUniMatNormLog(2,:)); hold on
plot(phiHRMid2InMatNormLog(2,:)); hold on
plot(phiHRMid1InMatNormLog(2,:)); hold on
title("Pivot (Log)")
legend("Uni", "2In", "1In")
subplot(2,3,6)
plot(phiHRMidUniMatNormLog(3,:)); hold on
plot(phiHRMid2InMatNormLog(3,:)); hold on
plot(phiHRMid1InMatNormLog(3,:)); hold on
title("Dark (Log)")
legend("Uni", "2In", "1In")

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
plot(muEffVec./muEffVec(2), muEffVec./muEffVec(2), '-+'); hold on
plot(muEffVec./muEffVec(2),gValUni./gValUni(2), '-+'); hold on
plot(muEffVec./muEffVec(2),gVal2In./gVal2In(2), '-+');
plot(muEffVec./muEffVec(2),gVal1In./gVal1In(2), '-+');
title(["Calculated \mu_{eff} vs. Simulated \mu_{eff}"; "Normalized by Pivot Value"]);
legend("Reference", "Uni", "2In", "1In");
ylabel("Simulated \mu_{eff} [1/cm]");
xlabel("Calculated \mu_{eff} [1/cm]");

function normMat = normMatf(mat)
    minVec = min(mat, [], 2);
    maxVec = max(mat, [], 2);
    span = maxVec - minVec;
    normMat = (mat - minVec) ./ span;
end


% To generate a simple density function for the adaptive process, we now 
% run a forward problem on this mesh, and use the resulting field to 
% determine the target node density. You can apply other criteria for the 
% adaptive remeshing, but the mechanism remains the same.