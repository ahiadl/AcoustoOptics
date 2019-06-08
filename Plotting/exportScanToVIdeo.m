close all;
clear all;
clc;

scan = load('./Results/Results-12-Jul-2018 18-25-19.mat');
load('./Results/Results-12-Jul-2018 12-52-34-params.mat');
videoFWriter = vision.VideoFileWriter('scan.avi','FrameRate',2,'FileFormat','AVI');
videoFWriter.Quality = 100;
fig = figure();
fig.Position = [100, 200, 780, 640]
ax = gca();
% scan.I = scan.I./(max(max(max(scan.I))));
axisName = ["X axis"; "Y axis"; "Z axis"];
limits = [min(min(min(scan.I))), max(max(max(scan.I)))];
dims = [params.Stages.Xlen, params.Stages.Ylen, params.Sample.Zlen];
dims(1) = 15;
dims(2) = 31;
params.Sample.Zbar = params.Sample.Zbar*1000;
for i = 1:3
    value = axisName(i);
    cla(ax);
    if value == 'X axis'
        data = flip(permute(scan.I(:,1,:), [3,1,2]),1);
        data = data./max(max(data));
        h = imagesc(ax,'XData', params.Stages.Ypos,'YData', params.Sample.Zbar,...
             'CData', data,  limits);
        xlabel(ax, 'Y[mm]');
        ylabel(ax, 'Z[mm]');
        title('Scan Along X axis')
        xlim([params.Stages.Ypos(1)-1, params.Stages.Ypos(end)+1])
        ylim([params.Sample.Zbar(1)-1, params.Sample.Zbar(end)+1])
    elseif value == 'Y axis'
        data = flip(permute(scan.I(1,:,:), [3,2,1]),1);
        data = data./max(max(data));
        h = imagesc(ax, 'XData', params.Stages.Xpos,'YData', params.Sample.Zbar,...
             'CData',data, limits);
        xlabel(ax, 'X[mm]');
        ylabel(ax, 'Z[mm]');
        title('Scan Along Y axis')
        xlim([params.Stages.Xpos(end)-1, params.Stages.Xpos(1)+1])
        ylim([params.Sample.Zbar(1)-1, params.Sample.Zbar(end)+1])
    else
        data =  permute(scan.I(:,:,1), [2,1,3]);
        data = data./max(max(data));
        h = imagesc(ax, 'XData', params.Stages.Ypos, 'YData', params.Stages.Xpos,...
             'CData', data,  limits);
        xlabel(ax, 'Y[mm]');
        ylabel(ax, 'X[mm]');
        title('Scan Along Z axis')
        xlim([params.Stages.Ypos(1)-1,   params.Stages.Ypos(end)+1])
        ylim([params.Stages.Xpos(end)-1, params.Stages.Xpos(1)+1])
    end
    axis (ax,'fill');
    ax.FontSize = 24;
    c = colorbar(ax);
    c.FontSize = 28;
    
    for j = 1:dims(i)
        if     value == 'X axis'
            set(h, 'CData', flip(permute(scan.I(:,j,:), [3,1,2]),1));
            title({'Scan Along X axis'; ['X=' num2str(params.Stages.Xpos(j)),'[mm]']})
        elseif value == 'Y axis'
            set(h, 'CData', flip(permute(scan.I(j,:,:), [3,2,1]),1));
            title({'Scan Along Y axis'; ['Y=' num2str(params.Stages.Ypos(j)),'[mm]']})
        else
            set(h, 'CData', permute(scan.I(:,:,j), [2,1,3]));
            title({'Scan Along Z axis'; ['Z=' num2str(params.Sample.Zbar(j)),'[mm]']})
        end
 
        c.Limits = limits;
        F = getframe(fig);
        step(videoFWriter,F.cdata);
    end
end

release(videoFWriter);                 
close(fig);