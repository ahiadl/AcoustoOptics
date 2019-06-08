    figure(); 
    plot(linspace(0, tTrainTrue, samplesInTrain)*1e6, pulseData)
    set(gca, 'FontSize', 20)
    title('Input US pulse')
    xlabel('t[\mus]')
    ylabel('Amplitude')
    xlim([0, 5])
    
    figure(); 
    plot(linspace(0, tTrainTrue, samplesInTrain)*1e6, rand(1,samplesInTrain))
    set(gca, 'FontSize', 20)
    title('Sampled Signal at Z_16')
    xlabel('t[\mus]')
    ylabel('Amplitude')
    xlim([48, 51.2])
    
    
    fs = 5e6;
    T = 3.2e-6;
    rep = 100;
    signalLength = T*rep;
    numOfSamples = ceil(signalLength*fs);
 
    fbar = linspace(-fs/2,fs/2,numOfSamples);
    signal = (((1/1)*sin(2*pi*1.25e6*linspace(0,signalLength,numOfSamples)))+rand(1,numOfSamples));
    
    figure(); 
    plot(fbar*1e-6 , fftshift( (1/numOfSamples)*fft( signal ) ))
    set(gca, 'FontSize', 24)
    title('F\{s_{zn}(t)\}', 'FontSize', 20)
    xlabel('f[MHz]', 'FontSize', 20)
    ylabel('Amplitude','FontSize', 20)
    xlim([-2, 2])