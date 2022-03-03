close all
clear all
clc;

%% Load Data
addpath('C:\Users\sahiadl.EED\Desktop\TH')
load('14-Oct-2021 23-09-59-Tamar_knot_ink');

%% Init BP
bp = BackProjection();


%% Config
uVarsBp = bp.createUserVars();

uVarsBP.c  = 1440;
uVarsBP.fs = csVars.fs;
uVarsBP.geo = 'Planar';
uVarsBP.imageWidth = 10e-3; %[m]
uVarsBP.mode = 3;
uVarsBP.tVec = csVars.tVec;

idx = 512;

uVarsBP.xSensorPos = csVars.scanVecBin; % horizontal position of the transducer
uVarsBP.ySensorPos = csVars.disc1Vec;
uVarsBP.zSensorPos = csVars.tVec(idx)*uVarsBP.c;

bp.setVars(uVarsBP);
vars = bp.getVars();


h = figure();
hAx = axes();
hIm = imagesc('XData', vars.recon.reconAxis*1e3,...
              'YData', vars.recon.reconAxis*1e3, ...
              'CData', zeros(length(vars.recon.reconAxis)));
axis equal tight

j=1;
N = 170;
dn = 2;  
recon = zeros(vars.recon.n, vars.recon.n, floor(N/dn));
for i= 1:dn:N
    uVarsBP.zSensorPos = csVars.tVec(i)*uVarsBP.c;
    bp.setVars(uVarsBP);
    recon(:,:,j) = bp.calcBP(resCs);
    set(hIm, 'CData', recon(:,:,j));
    drawnow();
    j=j+1;
end

reconTot = max(recon, [], 3);

hf = figure;
axRD = subplot(1,2,1);
hRD = imagesc(axRD, 'XData', 1:140, 'YData', csVars.tVec*1e6, 'CData', squeeze(resCs(85,:,1,1,:))');
hTitRD = title(axRD, 'Sinogram');
axis(axRD, 'tight')
axBP = subplot(1,2,2);
hBP = imagesc(axBP, 'XData', vars.recon.reconAxis ,'YData', vars.recon.reconAxis,'CData', reconTot);
hTitBP = title(axBP, 'Recon');
axis(axBP, 'equal', 'tight');

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