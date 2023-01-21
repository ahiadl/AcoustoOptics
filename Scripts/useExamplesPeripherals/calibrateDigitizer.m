close all
clear all
clc
instrreset

fg = fGen();
fg.connect();
digi = Digitizer();
digi.connect();

%% Config the Function Generator Signal:
fgClk = 100e6;

fGenVars = fg.uVarsCreate();

fGenVars.ch{1}.daqFilter     = 50;
fGenVars.ch{1}.amp           = 1;
fGenVars.ch{1}.bias          = 0;
fGenVars.ch{1}.triggerOwner  = true;
fGenVars.ch{1}.triggerWidth  = 512;
fGenVars.ch{1}.useExtSclkSrc = false;

fGenVars.ch{2}.daqFilter     = 50;
fGenVars.ch{2}.amp           = 1;
fGenVars.ch{2}.bias          = 0;
fGenVars.ch{2}.triggerOwner  = false;
fGenVars.ch{2}.triggerWidth  = 16;
fGenVars.ch{2}.useExtSclkSrc = false;

fGenVars.ch{1}.Sclk = fgClk;
fGenVars.ch{2}.Sclk = fgClk;

dataCh1 = fg.getChDataStruct();
dataCh2 = fg.getChDataStruct();


dataCh1.data = createSig(1);
dataCh1.dataLen = length(dataCh1.data);

%% Set Digitizer generic vars
digiVars = Digitizer.uVarsCreate();

timeToSample = 0.0002; %[s]
% bytesPerSample = 2;

digiVars.mode      = 'TS'; 
digiVars.useGPU    = false;

digiVars.triggerDelay = 0;
digiVars.extClk = true;
digiVars.extTrig = true;
digiVars.draw    = false;

digiVars.bufferSizeBytes = 8*(2^20);


digiVars.exportCropped = false;

fs = [5e6, 20e6; 20e6, 5e6];
ch = [1, 2, 4, 8, 16];

marker = ["-o", "-x", "-v", "-*", "-."];

%%
close all;

for i=1:size(fs,1)
    for j = 1:size(fs,2)
        curFS = fs(i,j);
        digiVars.samplesPerMeas  = timeToSample*curFS;
        figure();
        ax = axes(); hold(ax, 'on');
        title(sprintf("FS = %d[MHz]:  Phase: %d", curFS*1e-6, i));
        xlim(ax, [1,50]);
        for k = 1:5
            curCh = ch(k);
            dataCh2.data = createClkSig(curFS);
            
            dataCh2.dataLen = length(dataCh2.data);
            
            fg.reset();
            fg.setProperties(fGenVars.ch{1}, fGenVars.ch{2});
            fg.setData(dataCh1, dataCh2);
            fg.configChannel(1); % this will enable the channel
            fg.configChannel(2); % this will enable the channel
            
            digiVars.fs        = curFS;
            digiVars.channels  = curCh;
            
            digi.setVars(digiVars);
            digi.configure();
            
            curRes = digi.acquire();

            res{i,j,k,:} = curRes;
% 
%             set(hP, 'YData', curRes(1,:))
            plot(ax, 1:digiVars.samplesPerMeas, curRes(1,:), marker(k));
            legStr{k} = sprintf("Ch =%d", k);
        end
        legend(legStr);
        
    end
end


%%
close all;
fs = 5e6;

dataCh2.data = createClkSig(fs);
dataCh2.dataLen = length(dataCh2.data);

% digiVars.fs        = fs;
% digiVars.channels  = 1;
% digiVars.samplesPerMeas  = timeToSample*fs;
% 
% digi.setVars(digiVars);
% digi.configure();

cyc = 1:10;

dataCh1.data = createSig(cyc(1));
dataCh1.dataLen = length(dataCh1.data);  
fg.reset();
fg.setProperties(fGenVars.ch{1}, fGenVars.ch{2});
fg.setData(dataCh1, dataCh2);
fg.configChannel(1); % this will enable the channel
fg.configChannel(2); % this will enable the channel

figure();
ax2 = axes();
hp = plot(zeros(1,length(createSig(1))));
xlim(ax2, [0, 1000])

figure();
ax = axes(); hold(ax, 'on');
title(sprintf("FS = %d[MHz]", fs*1e-6));
xlim(ax, [1,500]);
for k = 1:length(cyc)
%     dataCh1.data = createSig(cyc(k));
%     dataCh1.dataLen = length(dataCh1.data);            
    
%     set(hp, 'YData', dataCh1.data)
% 
%     fg.reset();
%     fg.setProperties(fGenVars.ch{1}, fGenVars.ch{2});
%     fg.setData(dataCh1, dataCh2);
%     fg.configChannel(1); % this will enable the channel
%     fg.configChannel(2); % this will enable the channel

    curRes = digi.acquire();

    resCyc(k,:) = curRes;

    plot(ax, 1:digiVars.samplesPerMeas, curRes);
    legStr{k} = sprintf("Cyc =%d", cyc(k));
end
legend(legStr);

%%
function clkData = createClkSig(fs)
    % The data loaded to the AFG must be a multiplication of 16 in manners of data length.
    % creating sClk with 16 cycles meets this requirement no matter
    % the fFgClk used.
    fgClk = 100e6;
    sClkCycles  = 16;
    fsNaive     = fs;
    sClkDutyCyc = 50/100;
    
    fgClkSamplesPerSClkCyc           = fgClk/fs;
    
    % If there is no round number of fgClk samples in sClk cycle, 
    % slow the sClk to the period where fgClk fits in.
    if mod(fgClkSamplesPerSClkCyc,1)~=0
        fgClkSamplesPerSClkCyc    =  ceil(fgClkSamplesPerSClkCyc);
        fs                   = fgClk / fgClkSamplesPerSClkCyc;
        fprintf("Notice: The sampling frequence you have chosen cannot be genrated by the AFG. \n the closest sampling frequency is: %d\n", fs);
    end
    
    sClkT                = 1/fs;
    fgClkSamplesPerFsSig = fgClkSamplesPerSClkCyc * sClkCycles;
    
    %Create sClk Data
    cycleData                   = ones(fgClkSamplesPerSClkCyc,1);
    dutyCycleIdx                = floor(fgClkSamplesPerSClkCyc*(1-sClkDutyCyc))+1;
    cycleData(dutyCycleIdx:end) = 0;
    clkData                     = repmat(cycleData, sClkCycles, 1);
end

function sig = createSig(cyc)
    fgClk = 100e6;
    N = 20080;
    % AO Signal
    fsin = 1.25e6;
    % fSqnc = 5e3;
    delay = 0;
    
    sigTime = N/fgClk;
    NSin = fgClk/(fsin/cyc);
    
    % Create a sin
    tVecSin  = (1:1:NSin)*(1/fgClk);
    Nsin = length(tVecSin);
    
    tVecSinPad = (0:1:N-1)./fgClk;
    fVecSig = (fgClk/N) *  ( (-N/2) : 1 : (N/2)-1 ); 
    sig = sin(2*pi*fsin*tVecSin);
    
    delaySamples = floor(delay * fgClk);
    padLen = N-length(sig)-delaySamples;
    sig = [zeros(1,delaySamples), sig, zeros(1,padLen)];
    
    % PreTrigger Test Signal
    % sig = linspace(-1,1,N);
end
