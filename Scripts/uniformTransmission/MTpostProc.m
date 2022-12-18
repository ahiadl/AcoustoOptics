function ppData = postProc(userAux, data, vars)
%     slice    = data.resMat(:, userAux.idx, :, :, :);
%     sliceAC  = slice - mean(slice, 5);
%     sliceP2P = peak2peak(sliceAC, 5);
%     slicePa  = squeeze(sliceP2P / userAux.sensitivity);
%     ppData.slicePa = slicePa;
        
%     ppData = squeeze(peak2peak(data.resMat(:,:,:,:,195:end),5));
%     idx = userAux.idx;
%     
%     % AC Coupling
%     ppData.ac = permute(data.resMat(:,:,vars.idxs(2),:,:)...
%              - mean(data.resMat(:,:,vars.idxs(2),:,:), 5), [1,2,4,5,3]);
%     
%     % Extract Channel
%     ppData.chSignals = permute(ppData.ac(:,:,idx,:), [1,2,4,3]);
% %     ppData.chSignals = squeeze(ppData.ac(idx,:));
%     % Choose Pos
%     ppData.p2p = peak2peak(ppData.chSignals(:, :, 250:end), 3);
% %     [~, I] = max(ppData.p2p(:,vars.idxs(1)));
%     I=1;
%     ppData.sig = squeeze(ppData.chSignals(I,vars.idxs(1),:));
%       ppData.sig    = squeeze(data.ppData);
%       ppData.sigFFT = abs(ifftshift(fft(squeeze(data.ppData)))).^2;
%       ppData.dcMat =


end