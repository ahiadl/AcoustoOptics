function userPlotResults(plots, userData, vars)
%     set(plots.hP1, 'CData', userData);
%     set(plots.hP1, 'YData', userData.sig);
%     set(plots.hP2, 'CData', userData.p2p');

if sum(size(userData.sig)>1)==1
    set(plots.hP1, 'YData', userData.sig);
    set(plots.hP2, 'YData', userData.sigFFT);
    
end