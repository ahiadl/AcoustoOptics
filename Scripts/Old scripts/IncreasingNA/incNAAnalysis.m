dirName = 'D:\ResultsToKeep\Increasing NA';
thick = [200, 300, 400, 800, 1000];

%%
for j = 1:length(thick)
    
    measName = sprintf("%dum-woHad-2s", thick(j));

    filename = sprintf("%s/%s/AO-Vars.mat", dirName, measName);
    vars = load(filename);
    for i =1:30
        filename = sprintf("%s/%s/LiveAOResults/AO-Results-%d.mat", dirName, measName, i);
        tmpRes = load(filename);
        tmpRes.rawData = [];
        res(j,i) = tmpRes;
    end

    for i = 1:30
       idx(j,i)    = res(j,i).SNR.tailIdx;
       SNR(j,i)    = res(j,i).SNR;
       SNRval(j,i) = res(j,i).SNR.val;
    end
end

avgSNR = mean(SNRval, 2);
stdSNR = std(SNRval, 0, 2);

figure()
errorbar(thick, avgSNR, stdSNR)
xlim([100, 1200])
xlabel ("Core Thickness [\mum]");
ylabel ("Avg. SNR");
title("Single Pulse AOI")

%%
for j = 1:length(thick)
    
    measName = sprintf("%dum-wHad-2s", thick(j));

    filename = sprintf("%s/%s/AO-Vars.mat", dirName, measName);
    varsHad = load(filename);
    for i =1:30
        filename = sprintf("%s/%s/LiveAOResults/AO-Results-%d.mat", dirName, measName, i);
        tmpRes = load(filename);
        tmpRes.rawData = [];
        resHad(j,i) = tmpRes;
    end

    for i = 1:30
       idxHad(j,i)    = resHad(j,i).SNR.tailIdx;
       SNRHad(j,i)    = resHad(j,i).SNR;
       SNRHadval(j,i) = resHad(j,i).SNR.val;
    end
end

avgSNRHad = mean(SNRHadval, 2);
stdSNRHad = std(SNRHadval, 0, 2);

figure()
errorbar(thick, avgSNRHad, stdSNRHad)
xlim([100, 1200])
xlabel ("Core Thickness [\mum]");
ylabel ("Avg. SNR");
title("CT AOI");

%%
NA = [0.22 0.39 0.5];

for j = 1:length(NA)
    
    measName = sprintf("400um-woHad-2s-%.2f", NA(j));

    filename = sprintf("%s/%s/AO-Vars.mat", dirName, measName);
    vars = load(filename);
    for i =1:30
        filename = sprintf("%s/%s/LiveAOResults/AO-Results-%d.mat", dirName, measName, i);
        tmpRes = load(filename);
        tmpRes.rawData = [];
        resNA(j,i) = tmpRes;
    end

    for i = 1:30
       idxNA(j,i)    = resNA(j,i).SNR.tailIdx;
       SNRNA(j,i)    = resNA(j,i).SNR;
       SNRvalNA(j,i) = resNA(j,i).SNR.val;
    end
end

avgSNRNA = mean(SNRvalNA, 2);
stdSNRNA = std(SNRvalNA, 0, 2);

figure()
errorbar(NA, avgSNRNA, stdSNRNA)
xlim([0, 0.8])
xlabel ("Core Thickness [\mum]");
ylabel ("Avg. SNR");
title("Single Pulse AOI")

%%
NA = [0.22 0.39 0.5];

for j = 1:length(NA)
    
    measName = sprintf("400um-wHad-%.2f", NA(j));

    filename = sprintf("%s/%s/AO-Vars.mat", dirName, measName);
    varsHad = load(filename);
    for i =1:30
        filename = sprintf("%s/%s/LiveAOResults/AO-Results-%d.mat", dirName, measName, i);
        tmpRes = load(filename);
        resNAHad(j,i) = tmpRes;
    end

    for i = 1:30
       idxNAHad(j,i)    = resNAHad(j,i).SNR.tailIdx;
       SNRNAHad(j,i)    = resNAHad(j,i).SNR;
       SNRvalNAHad(j,i) = resNAHad(j,i).SNR.val;
    end
end

avgSNRNAHad = mean(SNRvalNAHad, 2);
stdSNRNAHad = std(SNRvalNAHad, 0, 2);

figure()
errorbar(NA, avgSNRNAHad, stdSNRNAHad)
xlim([0, 0.8])
xlabel ("FIber's NA [\mum]");
ylabel ("Avg. SNR");
title("CT AOI")

for i=1:30
    NA22AvgPhi(i,:) = resNAHad(1,i).phi;
    NA39AvgPhi(i,:) = resNAHad(2,i).phi;
    NA50AvgPhi(i,:) = resNAHad(3,i).phi;
end

avgPhi22 = mean(NA22AvgPhi,1);
avgPhi39 = mean(NA39AvgPhi,1);
avgPhi50 = mean(NA50AvgPhi,1);

figure();
subplot(1,3,1)
stem(avgPhi22)
title("NA=0.22");
subplot(1,3,2)
stem(avgPhi39)
title("NA=0.39");
subplot(1,3,3)
stem(avgPhi50)
title("NA=0.50");

for i=1:3
    for j = 1:30
        stdVec(i,j) = std(resNAHad(i,j).phi);
    end
end

stdVecAvg = mean(stdVec, 2);
stdVecStd = std(stdVec, 0, 2);

figure()
errorbar(NA, stdVecAvg, stdVecStd);
xlim([0, 0.8])
xlabel("NA");
ylabel("Raw Data STD");