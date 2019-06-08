function [U] = reshapeRawData(rawData, dim, ch)
% rawData should be given as:
% rows - complete period
% cols - samples
U = zeros(Npos, dim.Npos*dim.Npulse, ch);
for i=1:dim.Npos
    U(i,:,:) = reshape(rawData(:,(i-1)*dim.Npulse+1:i*dim.Npulse,:),1,dim.Npos*N,ch);
end

end

