    %%
    sClk = 100e6;
    N = 65536;
    sigTime = N/sClk;
    
    cycles = 20;
    sincFactor = 0.001;
    fsin = 1.25e6;
    % Create a sin
    tVecSin  = (0:1:((sClk/fsin)-1))*(1/sClk);
    Nsin = length(tVecSin);
    
    tVecSinPad = (0:1:N-1)./sClk;
    fVecSig = (sClk/N) *  ( (-N/2) : 1 : (N/2)-1 ); 
    sig = sin(2*pi*fsin*tVecSin);
    
    delay = 64.4e-6;
    
    delaySamples = floor(delay * sClk);
    padLen = N-length(sig)-delaySamples;
    sig = [zeros(1,delaySamples), sig, zeros(1,padLen)];
    
    figure();
    subplot(3,2,1);
    plot(tVecSinPad, sig)
    subplot(3,2,2);
    plot(fVecSig/1e6, abs(fftshift(fft(sig))).^2);
    xlim([-2.5, 2.5])
    %%
    % Create a sinc
    fsinc    = sincFactor*fsin;
%     tVecSinc = (0:1:(cycles*(sClk/fsinc)-1))*(1/sClk);
%     tVecSinc = tVecSinc-tVecSinc(floor(length(tVecSinc)/2));
    tVecSinc = (0:1:(Nsin-1))*(1/sClk);
    Nsinc    = length(tVecSinc);
    fVecSinc = (sClk/Nsinc) *  ( (-Nsinc/2) : 1 : (Nsinc/2)-1 ); 
    sinc     = sin(2*pi*fsinc*(tVecSinc))./(2*pi*fsinc*(tVecSinc));
    
    idx = find(isnan(sinc));
    sinc(idx)= 1;
    
    subplot(3,2,3);
    plot(tVecSinc, sinc)
    subplot(3,2,4)
    plot(fVecSinc/1e6, abs(fftshift(fft(sinc))).^2);
    xlim([-5, 5])
    %Create deformed signal


    
    
    
%     convSig = conv(sig, sinc, 'same');
    convSig = sig.*sinc;

    subplot(3,2,5);
    plot(tVecSinPad, convSig)
    subplot(3,2,6)
    plot(fVecSin/1e6, abs(fftshift(fft(convSig))).^2);
    xlim([-5, 5])
    

%     sig = [zeros(1,500), ones(1,1), zeros(1,500)];
%     sig = [zeros(1,500), 0.25, 1,  0.25, zeros(1,500)];
%     figure();
%     plot(sig)
%     padLen = 65600-length(sig);

h = squeeze(resCs);
H = fftshift(fft(h));










