close all
clear all
clc

addpath(genpath(pwd));

%%
fprintf("%s Starting Simulator\n", datestr(now,'HH:MM:SS'));

intDepthVec = [5, 10, 15];
numOfMid = length(intDepthVec);
slabWidth = 30; % [mm] half Width
slabDepth = 40; % [mm] full Depth

for k = 1:numOfMid
    %% Generate and load the Mesh
    toastThreadCount(20);

    meshName = sprintf("2Layers-%dx%dx%d-DxWxI-mm", slabDepth, slabWidth, intDepthVec(k));
    meshPath = "./Meshes";
    meshFilename = sprintf("%s.msh", meshName);
    meshFullName = sprintf("%s/%s", meshPath, meshFilename);
    resDirName = sprintf("E:/AcoustoOptics/Simulations/Results/2Layers-%dx%dx%d-DxWxI-mm-lightInputLayer", slabDepth, slabWidth, intDepthVec(k));

    %Create the mesh from the geometrical points defined in the .geo file
    fprintf("%s Constructing Mesh\n", datestr(now,'HH:MM:SS')); 
    if ~exist(meshFullName, 'file')
        fprintf("%s Can't locate the mesh. Generating...\n", datestr(now,'HH:MM:SS')); 
        str(1) = sprintf("intDepth = %d;", intDepthVec(k));
        str(2) = sprintf("depth = %d;", slabDepth);
        str(3) = sprintf("width = %d;", slabWidth);
        replaceLineInFile("./2LayerBox.geo", 1:3,  str);
        system (sprintf("gmsh -3 ./2LayerBox.geo -o %s", meshFullName));
    else
        fprintf("%s Located the mesh. Loading without generating.\n", datestr(now,'HH:MM:SS'));
    end

    % Load the mesh into toastMesh object.
    % Notice: If the .msh file was already constructed the above code is
    % redundant and the mesh should only be loaded with the code below.
    fprintf("%s Loading Mesh\n", datestr(now,'HH:MM:SS'));
    mesh = toastMesh(meshFullName,'gmsh');

    if ~exist(resDirName, 'dir')
        mkdir(resDirName)
    end

    %% Define The Mesh Grid
    dl      = 0.5;
    intDepth = intDepthVec(k);
    bbMesh  = mesh.BoundingBox;
    
    %Define axes lengths:
    xLen    = bbMesh(2,1) - bbMesh(1,1);  % length of x axis of the mesh -[mm]
    yLen    = bbMesh(2,2) - bbMesh(1,2);  % length of y axis of the mesh -[mm] 
    zLen    = bbMesh(2,3) - bbMesh(1,3);  % length of z axis of the mesh -[mm]
    
    % Create axes vectors:
    xVec = bbMesh(1,1):dl:bbMesh(2,1);
    yVec = bbMesh(1,2):dl:bbMesh(2,2);
    zVec = bbMesh(1,3):dl:bbMesh(2,3);
    
    % Calculate Indexes:
    xSize = length(xVec);
    ySize = length(yVec);
    zSize = length(zVec);

%     % Index of the axis center:
%     xMid = floor(xSize/2); 
%     yMid = floor(ySize/2);
%     zMid = floor(zSize/2);

    % Create mesh grid that will help cast the toastMesh onto.
    [X, Y, Z] = meshgrid (yVec, xVec, zVec);
    
    % Create rMat to describe distances of each point in the mesh from the
    % illumination point on the boundary (relevant for point source)
    rMat = sqrt(X.^2 + Y.^2 + Z.^2);

    grd = [xSize, ySize, zSize];
    
    % Create high resolution grid for uniform illumination simulation
    dlHR = dl/10;
    xVecHR = xVec;
    yVecHR = bbMesh(1,2):dlHR:bbMesh(2,2);
    zVecHR = bbMesh(1,3):dlHR:bbMesh(2,3);

    [XHR, YHR, ZHR] = meshgrid (yVecHR, xVecHR, zVecHR);

    xSizeHR = length(xVecHR);
    ySizeHR = length(yVecHR);
    zSizeHR = length(zVecHR);  
    
    yMidHRIdx = floor(ySizeHR/2)+1;
    zMidHRIdx = floor(zSizeHR/2)+1;

    mask2In = sqrt(ZHR.^2 + XHR.^2) <=  25.4;
    mask1In = sqrt(ZHR.^2 + XHR.^2) <=  12.7;

    % Layers Length and Indeces
    firstLen  = intDepth;
    secondLen = xLen - intDepth;
    
    intIdx = round((intDepth/xLen) *xSizeHR);
    
    firstSize  = intIdx;
    secondSize = xSizeHR - firstSize;
    
    slopeFirstIdx  = floor(firstSize/2);
    slopeSecondIdx = floor(secondSize/2);

    %% Create Toast Basis
    % Casting the coordinates defined in the previose section onto the mesh 
    % generated with gMesh 
    fprintf("%s Creating Basis\n", datestr(now,'HH:MM:SS'));
    basis = toastBasis (mesh,grd);
    
    ne         = mesh.ElementCount;
    nv         = mesh.NodeCount;
    regidx     = mesh.Region;
    regno      = unique(regidx);%sorted from small to largest.
    
    firstLayerIdx  = find(regidx == regno(1));
    secondLayerIdx = find(regidx == regno(2));
    %% Optical Properties:
    
    %Refractive Index:
    refFirst  = 1.37;
    refSecond = 1.37;
    
    %MuS Prime:   [1/mm]
    musFirst  = 0.4;
    musSecond = 0.4;
    
    %MuA (Pivot): [1/mm]
    muaPivotFirst  = 0.0184/4;
    muaPivotSecond = 0.0184;

    % Simulated Values:
%     muaFirstVec   = muaPivotFirst * [1/4, 1/2, 1, 2, 4];
    muaFirstVec   = muaPivotFirst;
    musFirstVec   = musFirst*ones(1,length(muaFirstVec));
    muEffFirstVec = sqrt(3*muaFirstVec.*(muaFirstVec + musFirstVec));

    muaSecondVec   = muaPivotSecond * [1/4, 1/2, 1, 2, 4];
%     muaSecondVec   = muaPivotSecond;
    musSecondVec   = musSecond*ones(1,length(muaSecondVec));
    muEffSecondVec = sqrt(3*muaSecondVec.*(muaSecondVec + musSecondVec));
    
    pivotIdxFirst    = find(muaFirstVec  == muaPivotFirst);
    pivotIdxSecond   = find(muaSecondVec == muaPivotSecond);
    muEffPivotFirst  = muEffFirstVec(pivotIdxFirst);
    muEffPivotSecond = muEffSecondVec(pivotIdxSecond);
    
    %% Source Properties
    srcVec = [0,0,0];
    mesh.SetQM(srcVec, []);
    qvec = mesh.Qvec('Neumann','Point',0);
    
    
    for j = 1:length(muEffFirstVec)
        %% Loop Over First Layer Parameters
        fprintf("%s Starts a New Iteration\n", datestr(now,'HH:MM:SS'));
        curMuaFirst   = muaFirstVec (j);
        curMusFirst   = musFirstVec(j);
        curMuEffFirst = muEffFirstVec(j);
        
        
        for i = 1:length(muEffSecondVec)
            %% Loop Over Second Layer Parameters
            curMuaSecond = muaSecondVec(i);
            curMusSecond = musSecondVec(i);
            muEffSecond  = sqrt(3*curMuaSecond*curMusSecond);

            ref                 = ones(ne,1);
            ref(firstLayerIdx)  = refFirst;
            ref(secondLayerIdx) = refSecond;

            mus                 = ones(ne,1);
            mus(firstLayerIdx)  = curMusFirst;
            mus(secondLayerIdx) = curMusSecond;

            mua                 = ones(ne,1);
            mua(firstLayerIdx)  = curMuaFirst;
            mua(secondLayerIdx) = curMuaSecond;

            % Creating The System Matrix
            fprintf("%s Creating SysMat\n", datestr(now,'HH:MM:SS'));
            K = dotSysmat(mesh, mua ,mus, ref, 'EL');

            % Solve
            fprintf("%s Solving\n" , datestr(now,'HH:MM:SS'));
            phi = K\qvec;
            phiMapped  = basis.Map('M->B',phi);

            fprintf("%s Analyzing\n", datestr(now,'HH:MM:SS'));
            % Extract Phi: Linear and Log
            phiGrd     = reshape(phiMapped, xSize, ySize, zSize);

            fprintf("%s Interpolating\n", datestr(now,'HH:MM:SS'));
            phiHR       = interp3(X, Y, Z, phiGrd, XHR, YHR, ZHR);
            
            fprintf("%s Extracting\n", datestr(now,'HH:MM:SS'));
            phiHRMid    = phiHR(:, :, zMidHRIdx);
            phiHRMidLog = log(phiHRMid);

            phiHRMidUni     =(sum(sum(phiHR,3),2));
            phiHRMidUniNorm = normMatf(phiHRMidUni,1);
            phiHRMidUniLog  = log(phiHRMidUniNorm);

            phiHRMid2In     = sum(sum(phiHR.*mask2In,3),2);
            phiHRMid2InNorm = normMatf(phiHRMid2In,1);
            phiHRMid2InLog  = log(phiHRMid2InNorm);

            phiHRMid1In     = sum(sum(phiHR.*mask1In,3),2);
            phiHRMid1InNorm = normMatf(phiHRMid1In,1);
            phiHRMid1InLog  = log(phiHRMid1InNorm);
            
            [gVecUniLog] = gradient(phiHRMidUniLog);
            gValUni(1) = abs(mean(gVecUniLog(1:slopeFirstIdx)))/dlHR;
            gValUni(2) = abs(mean(gVecUniLog(intIdx:slopeSecondIdx)))/dlHR;
            [gVec2InLog] = gradient(phiHRMid2InLog);
            gVal2In(1) = abs(mean(gVec2InLog(1:slopeFirstIdx)))/dlHR;
            gVal2In(2) = abs(mean(gVec2InLog(intIdx:slopeSecondIdx)))/dlHR;
            [gVec1InLog] = gradient(phiHRMid1InLog);
            gVal1In(1) = abs(mean(gVec1InLog(1:slopeFirstIdx)))/dlHR;
            gVal1In(2) = abs(mean(gVec1InLog(intIdx:slopeSecondIdx)))/dlHR;

            res.phi3D      = phiGrd;         %3D
            res.phi2D      = phiHRMid;       %2D
            res.phi2DLog   = phiHRMidLog;    %2D
            res.phiUni     = phiHRMidUni;    %1D
            res.phiUniNorm = phiHRMidUniNorm;%1D
            res.phiUniLog  = phiHRMidUniLog; %1D
            res.phi2In     = phiHRMid2In;    %1D
            res.phi2InNorm = phiHRMid2InNorm;%1D
            res.phi2InLog  = phiHRMid2InLog; %1D
            res.phi1In     = phiHRMid1In;    %1D
            res.phi1InNorm = phiHRMid1InNorm;%1D
            res.phi1InLog  = phiHRMid1InLog; %1D

            res.gVecUniLog        = gVecUniLog;
            res.gValUni           = gValUni;
            res.gVec2InLog        = gVec2InLog;
            res.gVal2In           = gVal2In;
            res.gVec1InLog        = gVec1InLog;
            res.gVal1In           = gVal1In;
            
            filename = sprintf("%s/%s-1L-%d-2L-%d.mat", resDirName,meshName, j, i);

            fprintf("%s Saving\n", datestr(now,'HH:MM:SS'));
            save(filename, '-Struct', 'res', '-v7.3')
            
            resTot.resArr{i}.phiHRMidUni    = phiHRMidUni;    %1D
            resTot.resArr{i}.phiHRMidUniLog = phiHRMidUniLog; %1D
            resTot.resArr{i}.phiHRMid2In    = phiHRMid2In;    %1D
            resTot.resArr{i}.phiHRMid2InLog = phiHRMid2InLog; %1D
            resTot.resArr{i}.phiHRMid1In    = phiHRMid1In;    %1D
            resTot.resArr{i}.phiHRMid1InLog = phiHRMid1InLog; %1D

            res = [];
            phiHR = [];
            phiHRLog = [];
            
            pause(1);
        end
        
        resTot.ref              = ref;
        resTot.muaFirstVec      = muaFirstVec; 
        resTot.musFirstVec      = musFirstVec;
        resTot.muEffFirstVec    = muEffFirstVec;
        resTot.muEffPivotFirst  = muEffPivotFirst;
        resTot.muaSecondVec     = muaSecondVec; 
        resTot.musSecondVec     = musSecondVec;
        resTot.muEffSecondVec   = muEffSecondVec;
        resTot.muEffPivotSecond = muEffPivotSecond;
        
        resTot.depthVec   = xVec;
        resTot.depthVecHR = xVecHR;

        % collectResSingleLayer;
        fprintf("%s Saving total results\n", datestr(now,'HH:MM:SS'));
        filename = sprintf("%s/%s-1L-%d-FullResults.mat", resDirName, meshName, j );
        save(filename, '-Struct', 'resTot', '-v7.3')
        fprintf("%s Done all combinations\n", datestr(now,'HH:MM:SS'));

        resTot = [];

        pause(1);
    end
    
    clear mesh
    clear basis
    clear qvec
    
    fprintf("%s Done simulations\n", datestr(now,'HH:MM:SS'));
    
    pause(1);
end