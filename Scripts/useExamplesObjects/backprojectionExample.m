close all
clear all
clc;

%% Load Data
% addpath('C:\Users\sahiadl.EED\OneDrive - Technion\Graduate\Measurements\OA-Finger\ToOmer')
load('C:\Users\sahiadl.EED\OneDrive - Technion\Graduate\Measurements\OA-Finger\fingerNew.mat');
% load('finger.mat');

% addpath('E:\TH\Nov-2021 Knot Targets')
% load('30-Nov-2021 21-29-54-Tamar_knot2_z_255_agar_50micron5_5_bin232');

for i=1:330
    sigMat(:,:,i) = res{i}.sigMat;
end


%% Init BP
clc
bp = BackProjection();


% Config
uVarsBP = bp.createUserVars();

uVarsBP.c  = 1440;
uVarsBP.fs = 40e6;
uVarsBP.geo = 'Circular';
uVarsBP.imageWidth = 30e-3; %[m]
uVarsBP.mode = 3;
uVarsBP.tVec = (1/40e6)*(0:1:2029);

idx = 512;

% uVarsBP.xSensorPos = csVars.scanVecBin; % horizontal position of the transducer
% uVarsBP.ySensorPos = csVars.disc1Vec;
% uVarsBP.zSensorPos = csVars.tVec(idx)*uVarsBP.c;

bp.setVars(uVarsBP);
vars = bp.getVars();

h = figure();
hAx = axes();
hIm = imagesc('XData', vars.recon.reconAxis*1e3,...
              'YData', vars.recon.reconAxis*1e3, ...
              'CData', zeros(length(vars.recon.reconAxis)));
hTit = title("title");

axis equal tight

j=1;
N = size(sigMat, 3);
dn = 2;  
recon = zeros(vars.recon.n, vars.recon.n, N);

sigMat2 = permute(sigMat, [4,5,6,2,1,3]);
% sigMat2 = sigMat;
for i= 1:N
%     uVarsBP.zSensorPos = uVarsBP.tVec(i)*uVarsBP.c;
%     bp.setVars(uVarsBP);
    recon(:,:,i) = bp.calcBP(sigMat2(:,:,:,:,:,i));
%     recon(:,:,i) = bp.calcBP(sigMat(:,:,i));
    set(hIm, 'CData', recon(:,:,i));
    set(hTit, 'String', sprintf("Image: %d", i));
    drawnow();
    pause(0.2);
%     j=j+1;
end

% reconTot = max(recon, [], 3);

% hf = figure;
% axRD = subplot(1,2,1);
% hRD = imagesc(axRD, 'XData', 1:140, 'YData', uVarsBP.tVec*1e6, 'CData', squeeze(resMat(85,:,1,1,:))');
% hTitRD = title(axRD, 'Sinogram');
% axis(axRD, 'tight')
% axBP = subplot(1,2,2);
% hBP = imagesc(axBP, 'XData', vars.recon.reconAxis ,'YData', vars.recon.reconAxis,'CData', reconTot);
% hTitBP = title(axBP, 'Recon');
% axis(axBP, 'equal', 'tight');

% A = squeeze(resCs(70,:,1,1,:))';
% A(:,66) = 0.5*(A (:, 67)+ A (:, 65));
% A(:,36) = 0.5*(A (:, 37)+ A (:, 35));
% figure();
% imagesc(A)

% 
% for i =1:size(finger.sigMat, 4)
%     sigMat = single(squeeze(finger.sigMat(:,:,1,i)));
%     set(hRD, 'CData', sigMat);
% 
%     set(hBP, 'CData', recon);
% %     pause(1)
%     drawnow()
% end