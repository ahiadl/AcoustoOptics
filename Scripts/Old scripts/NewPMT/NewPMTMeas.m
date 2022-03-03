maxVal = [];
avgBkg = [];
stdBkg = [];
SNR = [];

peakIdx = 80;
len = 127;
minIdx = 22;

peakIdxShort = peakIdx - minIdx;

figure();
h=stem(1,1);

%%
aoVars = ao.getVars();
clear res
resMeas = [];
maxVal  = [];
avgBkg  = [];
stdBkg  = [];
SNR     = [];

if ~isgraphics(h)
    figure();
    h=stem(1,1);
end

% dcLevel = [0.5:0.25:1, 1.5:0.5:2.5];
% dcLevel = [-0.05:0.05:0.2, 0.25:0.25:1, 1.5:0.5:2.5];
quant = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05];
dcLevel = 0;
for j=1:length(dcLevel)
%     d = dialogBox(dcLevel(j));
    uiwait(d)
    
    maxVal = [];
    avgBkg = [];
    stdBkg = [];
    SNR = [];

    for i=1:30
%         fprintf("Meas: %d\n", );
        res(j,i)      = ao.runAcoustoOptics();
        phi = res(j,i).rawPhi(minIdx:end);
        bkg = phi (1:peakIdxShort);
        maxVal(j,i)   = max(phi);
        avgBkg(j,i)   = mean(bkg);
        stdBkg(j,i)   = std(bkg);
        SNR(j,i)      = (maxVal(j,i)-avgBkg(j,i))/stdBkg(j,i);
        normNoise(j,i) = res(1,i).normNoise;
        set(h, 'XData', 1:length(phi), 'YData', phi);
        drawnow;
        fprintf("Meas: %.2f : %d, SNR: %.2f\n", dcLevel(j), i, SNR(j,i));
    end

    resMeas(j).res = res(j,:);
    resMeas(j).maxVal  = maxVal(j,:);
    resMeas(j).avgBkg  = avgBkg(j,:);
    resMeas(j).stdBkg = stdBkg(j,:);
    resMeas(j).SNR = SNR(j,:);
    resMeas(j).avgSNR = mean(SNR(j,:));
    resMeas(j).stdSNR = std(SNR(j,:));
    resMeas(j).normNoise = mean(normNoise(j,:));

    SNRAvg(j) = resMeas(j).avgSNR;
    SNRstd(j) = resMeas(j).stdSNR;
    fprintf("Results: Average SNR: %.2f, STD: %.2f\n", SNRAvg(j), SNRstd(j));

    timeStamp = strrep(datestr(datetime('now')),':','-');
    filename = sprintf('../Measurements/PDFinal-Exp/%s-PMT-2s-q.%.3f.mat', timeStamp, qunat(j));
    save(filename, 'resMeas', 'aoVars', '-v7.3');
end 

timeStamp = strrep(datestr(datetime('now')),':','-');
filename = sprintf('../Measurements/PDFinal-Exp/%s-PMT-2s-q.%.3f.mat', timeStamp, q);
save(filename, 'resMeas', 'aoVars', '-v7.3');

%%
figure()
errorbar(dcLevel, SNRAvg,SNRstd)

% save('../Measurements/07-07-21-PD-FreeSpace-2s.mat', 'resPD')