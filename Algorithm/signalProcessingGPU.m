function [ZpowerSpectrum] = signalProcessing(U, sample)
%SIGNALPROCESSING Summary of this function goes here
%   Detailed explanation goes here
I = fft(U,2);
ZpowerSpectrum = mean(abs(I(:,sample.USidx)).^2,3);
end

