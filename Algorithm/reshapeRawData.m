function [U] = reshapeRawData(rawData, Sample)
% rawData should be given as:
% cols - complete period
% rows - samples
U = zeros(Sample.Npos, Sample.actualPeriodsSampled*Sample.Npulse, Sample.channels);

%Chopping the samples representing the 1st pulse before it got into the
%phantom and the corresponding residue

for i=1:Sample.Npos
    U(i,:,:) = reshape(rawData((i-1)*Sample.Npulse+1:i*Sample.Npulse,:,:),1,Sample.Npulse*Sample.actualPeriodsSampled ,Sample.channels);
end

end

