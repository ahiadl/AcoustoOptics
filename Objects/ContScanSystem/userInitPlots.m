function plots = userInitPlots(vars)

    plots.hFig = figure();
    plots.ax1 = subplot(1,2,1);
    plots.ax2 = subplot(1,2,2);
    
    tmp = zeros(1,vars.numOfSamples);
    plots.hP1 = plot(plots.ax1, vars.tVec*1e6, tmp);
    xlabel("t[\mu s]");
    ylabel ("Signal [mV]");
%     xlim(plots.ax1, [30,90])
    
    plots.hP2 = plot(plots.ax2, vars.userAux.fVec/1e6, tmp);
                    
    xlabel("frequency [MHz]");
    ylabel ("PS [AU]");
    xlim(plots.ax2, [-5,5])

end