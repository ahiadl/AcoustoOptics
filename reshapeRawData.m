function [Zres, U] = reshapeRawData(rawData, dim)
% rawData should be given as:
% rows - complete period
% cols - samples

for i=1:Npos
    U(i,:) = reshape(rawData(:,(i-1)*dim.Npulse+1:i*dim.Npulse),1,dim.Npos*N);
end

end

