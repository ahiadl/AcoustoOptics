function [I, V, fftRes] = signalProcessing(U, Sample,Plotting)
% U - (z, ptsPerPos, ch)
fftRes = 2*fft(U,[],2)./Sample.samplesPerPos;
% sig = I(Plotting.Zcoor,:,Plotting.channelToPlot);
V = abs(fftRes(:,Sample.USidx,:));
I = mean(V,3);

% figure()
% plot(Sample.frequencyBarShifted, fftshift(abs(I(:,:,1)),2));

end

