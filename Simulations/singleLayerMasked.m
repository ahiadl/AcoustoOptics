close all;
clear all;
clc;

addpath(genpath(pwd))

%% Mesh

meshWidth = [30, 40, 50];
meshDepth = [5, 10, 15, 20, 30, 40, 50, 60];

for k = 1:length(meshWidth)
    for j = 1:length(meshDepth)

        toastThreadCount(20);

        meshName = sprintf("1LayerSlabMesh-%dx%d-DxW-mm",  meshDepth(j), meshWidth(k));
        meshPath = "./Meshes";
        meshFilename = sprintf("%s.msh", meshName);
        meshFullName = sprintf("%s/%s", meshPath, meshFilename);
        resDirName = sprintf("E:/AcoustoOptics/Simulations/Results/%dx%d-DxW-mm", meshDepth(j), meshWidth(k));

        %Create the mesh from the geometrical points defined in the .geo file
        fprintf("Constructing Mesh\n"); 
        if ~exist(meshFullName, 'file')
            fprintf("Can't locate the mesh. Generating...\n"); 
            str(1) = sprintf("depth = %d;", meshDepth(j));
            str(2) = sprintf("width = %d;", meshWidth(k));
            replaceLineInFile("./1LayerBox.geo", [1,2],  str);
            system (sprintf("gmsh -3 ./1LayerBox.geo -o %s", meshFullName));
        else
            fprintf("Located the mesh. Loading without generating.\n");
        end

        % Load the mesh into toastMesh object.
        % Notice: If the .msh file was already constructed the above code is
        % redundant and the mesh should only be loaded with the code below.
        fprintf("Loading Mesh\n");
        tic;
        mesh = toastMesh(meshFullName,'gmsh');
        toc
        fprintf("Loading Mesh\n");
        
        if ~exist(resDirName, 'dir')
            mkdir(resDirName)
        end
        %% Define The Mesh Grid
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

        % Create axes vectors:
        xVec = bbMesh(1,1):dl:bbMesh(2,1);
        yVec = bbMesh(1,2):dl:bbMesh(2,2);
        zVec = bbMesh(1,3):dl:bbMesh(2,3);

        % Calculate Indexes:
        xSize = length(xVec);
        ySize = length(yVec);
        zSize = length(zVec);

%         % Index of the axis center:
%         xMid = floor(xSize/2); 
%         yMid = floor(ySize/2);
%         zMid = floor(zSize/2);

        % Create mesh grid that will help cast the toastMesh onto.
        [X, Y, Z] = meshgrid (yVec, xVec, zVec);

        % Create rMat to describe distances of each point in the mesh from the
        % illumination point on the boundary (relevant for point source)
        rMat = sqrt(X.^2 + Y.^2 + Z.^2);

        grd = [xSize, ySize, zSize];

        %Create high resolution grid for uniform illumination simulation
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

        mask2In = sqrt(ZHR.^2 + XHR.^2) <  25.4;
        mask1In = sqrt(ZHR.^2 + XHR.^2) <  12.7;

        slopeIdx = floor(xSizeHR/2);

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

        %% Optical Properties:
        
        % Refractive Index:
        ref  = 1.37;
        
        % MuS Prime: [1/mm]
        mus  = 0.4; %Prime(!) 1/mm
        
        % MuA (Pivot): [1/mm]
        muaPivot  = 0.0184; % 1/mm  
        
        % Simulated Values:
        muaVec   = [muaPivot/4, muaPivot/2, muaPivot, muaPivot*2, muaPivot*4];
        musPVec   = mus*ones(1,length(muaVec));
        muEffVec = sqrt(3*muaVec.*(muaVec+musPVec));
        
        pivotIdx = find(muaVec == muaPivot);
        muEffPivot = muEffVec(pivotIdx);
        
        %% Source Properties
        srcVec = [0,0,0];
        mesh.SetQM(srcVec, []);
        qvec = mesh.Qvec('Neumann','Point',0);
        %%
        for i =1 :length(muEffVec)
            fprintf("Starts a New Iteration: %d/%d\n", i, length(muEffVec));
            curMua  = muaVec (i);
            curMus = musPVec(i);

            muEffCalc     = sqrt(3*(curMua+curMus)*curMua);

            ref                 = ones(ne,1);
            ref(firstLayerIdx)  = ref;

            mus                 = ones(ne,1);
            mus(firstLayerIdx)  = curMus;

            mua                 = ones(ne,1);
            mua(firstLayerIdx)  = curMua;

            % Creating The System Matrix
            fprintf("Creating SysMat\n");
            tic
            K = dotSysmat(mesh, mua ,mus, ref, 'EL');
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

            fprintf("Interpolating\n");
            tic
            phiHR       = interp3(X, Y, Z, phiGrd, XHR, YHR, ZHR);
            toc

            fprintf("Extracting\n");
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
            gValUni = abs(mean(gVecUniLog(1:slopeIdx)))/dlHR;
            [gVec2InLog] = gradient(phiHRMid2InLog);
            gVal2In = abs(mean(gVec2InLog(1:slopeIdx)))/dlHR;
            [gVec1InLog] = gradient(phiHRMid1InLog);
            gVal1In = abs(mean(gVec1InLog(1:slopeIdx)))/dlHR;

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

            filename = sprintf("%s/%s - %d.mat", resDirName,meshName, i);

            fprintf("Saving\n");
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

        resTot.ref        = ref;
        resTot.muaVec     = muaVec; 
        resTot.musVec     = musPVec;
        resTot.muEffVec   = muEffVec;
        resTot.muEffPivot = muEffPivot;
        resTot.depthVec   = xVec;
        resTot.depthVecHR = xVecHR;

        % collectResSingleLayer;
        fprintf("Saving total results\n");
        filename = sprintf("%s/%s - FullResults.mat", resDirName, meshName);
        save(filename, '-Struct', 'resTot', '-v7.3')
        fprintf("Done simulations\n");

        resTot = [];

        pause(1);
    end
    
    clear mesh
    clear basis
    clear qvec
    
    pause(1);
end