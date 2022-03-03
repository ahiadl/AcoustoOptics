close all
clear all
clc

res  = load('D:\Results\06-Dec-2021 18-35-28-Scan2D\2DScan-Results.mat');
vars = load('D:\Results\06-Dec-2021 18-35-28-Scan2D\2DScan-Vars.mat');

debug = false;
%%
depth = 75;

phiPDRaw  = squeeze(res.phi(:,:,:,:,1));
phiPMTRaw = squeeze(res.phi(:,:,:,:,2));

usAx    = vars.grid.depthVec*1e3;
scanAx  = vars.grid.scanVec;

if debug
    figure(); 
    ax1 = subplot(2,1,1); 
    imagesc(usAx, scanAx, db(phiPDRaw));
    colorbar
    axis tight equal 
    ax2 = subplot(2,1,2); 
    imagesc(usAx, scanAx, db(phiPMTRaw))
    colorbar
    axis tight equal 
    linkaxes([ax1, ax2])
end
%% Pre Processing

% Normalizing
phiPDNorm = ( phiPDRaw - min(phiPDRaw(:)) ) / ( max(phiPDRaw(:)) - min(phiPDRaw(:)) );
phiPMTNorm = ( phiPMTRaw - min(phiPMTRaw(:)) ) / ( max(phiPMTRaw(:)) - min(phiPMTRaw(:)) );

if debug
    figure(); 
    ax1 = subplot(2,1,1); 
    imagesc(usAx, scanAx, db(phiPDNorm));
    colorbar
    axis tight equal 
    ax2 = subplot(2,1,2); 
    imagesc(usAx, scanAx, db(phiPMTNorm))
    colorbar
    axis tight equal 
    linkaxes([ax1, ax2])

    figure(); 
    ax1 = subplot(2,1,1); 
    imagesc( db(phiPDNorm));
    colorbar
    axis tight equal 
    ax2 = subplot(2,1,2); 
    imagesc( db(phiPMTNorm))
    colorbar
    axis tight equal 
    linkaxes([ax1, ax2])
end

%AC Coupling

phiPDAC  = phiPDRaw  - mean(phiPDRaw(:, 59:116),2);
phiPMTAC = phiPMTRaw - mean(phiPMTRaw(:, 59:116),2);

phiPDACNorm  = ( phiPDAC  - min(phiPDAC(:)) )  / ( max(phiPDAC(:))  - min(phiPDAC(:)) );
phiPMTACNorm = ( phiPMTAC - min(phiPMTAC(:)) ) / ( max(phiPMTAC(:)) - min(phiPMTAC(:)) );

% phiPDAC = phiPDRaw;
% phiPMTAC = phiPMTRaw;

climLog(1) = -50;
climLog(2) = max (max(phiPDACNorm(:)), max(phiPMTACNorm(:)));

if debug
    figure(); 
    ax1 = subplot(2,1,1); 
    imagesc(usAx, scanAx, db(phiPDACNorm), climLog);
    colorbar
    axis tight equal 
    ax2 = subplot(2,1,2); 
    imagesc(usAx, scanAx, db(phiPMTACNorm), climLog)
    colorbar
    axis tight equal 
    linkaxes([ax1, ax2])
    colormap hot
end

% Interpulation

[X,Y] = meshgrid(usAx, scanAx);

xq = usAx(1) : 0.05 : usAx(end);
yq = scanAx(1) : 0.05 : scanAx(end);

[Xq,Yq] = meshgrid(xq, yq);

phiPDInterp = interp2(X, Y, phiPDACNorm, Xq, Yq, 'cubic'); 
phiPMTInterp = interp2(X, Y, phiPMTACNorm, Xq, Yq, 'cubic'); 

phiPDInterpNorm  = ( phiPDInterp  - min(phiPDInterp(:)) )  / ( max(phiPDInterp(:))  - min(phiPDInterp(:)) );
phiPMTInterpNorm = ( phiPMTInterp - min(phiPMTInterp(:)) ) / ( max(phiPMTInterp(:)) - min(phiPMTInterp(:)) );

climLog(1) = -50;
climLog(2) = max (max(phiPDInterpNorm(:)), max(phiPMTInterpNorm(:)));

if debug
    figure(); 
    ax1 = subplot(2,1,1); 
    imagesc(usAx, scanAx, db(phiPDInterpNorm), climLog);
    colorbar
    axis tight equal 
    ax2 = subplot(2,1,2); 
    imagesc(usAx, scanAx, db(phiPMTInterpNorm), climLog)
    colorbar
    axis tight equal 
    linkaxes([ax1, ax2])
    colormap hot
end

% Rotating
angle = 5;
phiPDRot  = imrotate(phiPDInterp, angle);
phiPMTRot = imrotate(phiPMTInterp, angle);

if debug
    figure(); 
    ax1 = subplot(2,1,1); 
    imagesc(db(phiPDRot));
    colorbar
    axis tight equal 
    ax2 = subplot(2,1,2); 
    imagesc(db(phiPMTRot))
    colorbar
    axis tight equal 
    linkaxes([ax1, ax2])
    colormap hot
end

% Cropping
xMin = 80;
xMax = 3986;
yMin = 350;
yMax = 898;

phiPDCrop1 = phiPDRot(yMin:yMax,xMin:xMax) ;
phiPMTCrop1 = phiPMTRot(yMin:yMax,xMin:xMax) ;

if debug
    figure(); 
    ax1 = subplot(2,1,1); 
    imagesc(db(phiPDCrop1));
    colorbar
    axis tight equal 
    ax2 = subplot(2,1,2); 
    imagesc(db(phiPMTCrop1))
    colorbar
    axis tight equal 
    linkaxes([ax1, ax2])
end

xMin = 175;
xMax = 825;
yMin = 103;
yMax = size(phiPDCrop1,1);

phiPDCrop2  = phiPDCrop1(yMin:yMax,xMin:xMax) ;
phiPMTCrop2 = phiPMTCrop1(yMin:yMax,xMin:xMax) ;

phiPDCrop2Norm  = ( phiPDCrop2  - min(phiPDCrop2(:)) )  / ( max(phiPDCrop2(:))  - min(phiPDCrop2(:)) );
phiPMTCrop2Norm = ( phiPMTCrop2 - min(phiPMTCrop2(:)) ) / ( max(phiPMTCrop2(:)) - min(phiPMTCrop2(:)) );

usAxHR = xq(1:size(phiPDCrop2Norm,2));
usAxHR = (usAxHR - min(usAxHR))*cos(deg2rad(angle));
scanAxHR = yq(1:size(phiPDCrop2Norm,1));
scanAxHR = scanAxHR - min(scanAxHR)*cos(deg2rad(angle));

%% Plot
climLog(1) = -50;
climLog(2) = max (max(phiPDCrop2Norm(:)), max(phiPMTCrop2Norm(:)));

x1 = 1;

hFig1 = figure(); 
ax1 = axes();
imagesc(ax1, usAxHR, scanAxHR, db(phiPDCrop2Norm), climLog);
hCB = colorbar;
axis tight equal 
x = [x1, x1+5];
y = [20, 20];
line(x,y,'LineWidth',3,'Color',[1,1,1]);
text(x1,y(1)-1,'5 mm','FontWeight','bold','FontSize', 18,'Color',[1,1,1]);
text(x1,y(1)+1.5,'Homodyne AOI','FontWeight','bold','FontSize', 18,'Color',[1,1,1]);
text(x1,1.75,'(b)','FontWeight','bold','FontSize', 22,'Color',[1,1,1]);
colormap hot
set(ax1, 'XTick', [], 'YTick', [])
set(ax1, 'FontSize', 15)
ylabel(hCB,'Normalized Fluence [dB]')

hFig2 = figure(); 
ax2 = axes();
imagesc(ax2, usAxHR, scanAxHR, db(phiPMTCrop2Norm), climLog);
% hCB = colorbar;
axis tight equal 
x = [x1, x1+5];
y = [20, 20];
line(x,y,'LineWidth',3,'Color',[1,1,1]);
text(x1,y(1)-1,'5 mm','FontWeight','bold','FontSize', 18,'Color',[1,1,1]);
text(x1,y(1)+1.5,'Conventional AOI','FontWeight','bold','FontSize', 18,'Color',[1,1,1]);
text(x1,1.75,'(a)' ,'FontWeight','bold','FontSize', 22,'Color',[1,1,1]);
colormap hot
set(ax2, 'XTick', [], 'YTick', [])
% set(ax2, 'FontSize', 15)
% ylabel(hCB,'dB')

filename = '..\..\Publications\Homodyne AO\Figures\2DScanPDNew.pdf';
set(hFig1, 'Units', 'centimeters');
figPosCM = get(hFig1, 'Position');
hFig.PaperUnits = 'centimeters';
hFig.PaperType = '<custom>';
hFig.PaperSize = [figPosCM(3), figPosCM(4)];
hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
print(hFig1, filename, '-dpdf', '-r0')

filename = '..\..\Publications\Homodyne AO\Figures\2DScanPMTNew.pdf';
set(hFig2, 'Units', 'centimeters');
figPosCM = get(hFig2, 'Position');
hFig.PaperUnits = 'centimeters';
hFig.PaperType = '<custom>';
hFig.PaperSize = [figPosCM(3), figPosCM(4)];
hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
print(hFig2, filename, '-dpdf', '-r0')


hFig1 = figure(); 
ax1 = axes();
imagesc(ax1, usAxHR, scanAxHR, db(phiPDCrop2Norm), climLog);
hCB = colorbar;
axis tight equal 
x = [x1, x1+5];
y = [20, 20];
line(x,y,'LineWidth',3,'Color',[1,1,1]);
text(x1,y(1)-1,'5 mm','FontWeight','bold','FontSize', 18,'Color',[1,1,1]);
% text(x1,y(1)+1.5,'Homodyne AOI','FontWeight','bold','FontSize', 18,'Color',[1,1,1]);
% text(x1,1.75,'(b)','FontWeight','bold','FontSize', 22,'Color',[1,1,1]);
colormap hot
set(ax1, 'XTick', [], 'YTick', [])
set(ax1, 'FontSize', 15)
ylabel(hCB,'Normalized Fluence [dB]')

hFig2 = figure(); 
ax2 = axes();
imagesc(ax2, usAxHR, scanAxHR, db(phiPMTCrop2Norm), climLog);
% hCB = colorbar;
axis tight equal 
x = [x1, x1+5];
y = [20, 20];
line(x,y,'LineWidth',3,'Color',[1,1,1]);
text(x1,y(1)-1,'5 mm','FontWeight','bold','FontSize', 18,'Color',[1,1,1]);
% text(x1,y(1)+1.5,'Conventional AOI','FontWeight','bold','FontSize', 18,'Color',[1,1,1]);
% text(x1,1.75,'(a)' ,'FontWeight','bold','FontSize', 22,'Color',[1,1,1]);
colormap hot
set(ax2, 'XTick', [], 'YTick', [])
% set(ax2, 'FontSize', 15)
% ylabel(hCB,'dB')

filename = '..\..\Publications\Homodyne AO\Figures\2DScanPDNew2.pdf';
set(hFig1, 'Units', 'centimeters');
figPosCM = get(hFig1, 'Position');
hFig.PaperUnits = 'centimeters';
hFig.PaperType = '<custom>';
hFig.PaperSize = [figPosCM(3), figPosCM(4)];
hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
print(hFig1, filename, '-dpdf', '-r0')

filename = '..\..\Publications\Homodyne AO\Figures\2DScanPMTNew2.pdf';
set(hFig2, 'Units', 'centimeters');
figPosCM = get(hFig2, 'Position');
hFig.PaperUnits = 'centimeters';
hFig.PaperType = '<custom>';
hFig.PaperSize = [figPosCM(3), figPosCM(4)];
hFig.PaperPosition = [0, 0, figPosCM(3), figPosCM(4)];
print(hFig2, filename, '-dpdf', '-r0')
