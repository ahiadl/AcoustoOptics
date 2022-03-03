function [res] = SNRanalysis(sig, time)
dts = time.dt;
fs = 1/dts;
len = size(sig,2);

sigFFT = fftshift(fft(sig,len,2),2);

fBar = (fs/len)*((-len/2):1:(len/2)-1);
DCIdx = len/2+1;

avgFFT = mean((2/len)*abs(sigFFT),1);
avgFFTNoDC = avgFFT;
avgFFTNoDC(DCIdx) = 0;
avgFFTNoDCNorm = (avgFFTNoDC - min(avgFFTNoDC))/abs(max(avgFFTNoDC) - min(avgFFTNoDC));
avgFFTDB =  db(avgFFTNoDCNorm);

usIdx   = (abs(fBar - 1.25e6) < 1);
USVal   = avgFFTNoDC(usIdx);
USValDB = avgFFTDB(usIdx);
stdUS = std((2/len)*abs(sigFFT( :, usIdx )));

hFig = figure();
subplot(1,2,1)
plot(fBar/1e6, avgFFTNoDC);hold on
plot(fBar(usIdx)/1e6,USVal, 'g+')
xlabel("f[MHz]")
ylabel("Spectral Response (Linear)")
subplot(1,2,2)
plot(fBar/1e6, avgFFTDB);hold on
plot(fBar(usIdx)/1e6, USValDB, 'g+')
xlabel("f[MHz]")
ylabel("Spectral Response (dB)")

res.sig            = sig; 
res.sigFFT         = sigFFT;
res.avgFFT         = avgFFT;
res.avgFFTNoDC     = avgFFTNoDC;
res.avgFFTNoDCNorm = avgFFTNoDCNorm;
res.avgFFTDB       = avgFFTDB;
res.fBar           = fBar;
res.fs             = fs;
res.dts            = dts;
res.len            = len;
res.DCIdx          = DCIdx;
res.usIdx          = usIdx;
res.USVal          = USVal;
res.USValDB        = USValDB;
res.stdUS          = stdUS;
res.fig            = hFig;

end

