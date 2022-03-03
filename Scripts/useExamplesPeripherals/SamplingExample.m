close all
clear all
clc
instrreset
%% Function Generator

fg = fGen();
fg.connect();
digi = Digitizer();
digi.connect();
%%
fg.reset();
fGenVars = fg.uVarsCreate();

sClk = 100e6;

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

fGenVars.ch{1}.Sclk = sClk;
fGenVars.ch{2}.Sclk = sClk;

dataCh1 = fg.getChDataStruct();
dataCh2 = fg.getChDataStruct();

% NOTE: number of samples in the data for each channel should be a
% multiplication of 16.
N = 2^16;
dataCh1.data = [-1*ones(1,N - 60000), linspace(-1,1,60000)];
% dataCh1.data = linspace(-1,1,65600);
dataCh1.dataLen = length(dataCh1.data);

fg.setProperties(fGenVars.ch{1}, fGenVars.ch{2});
fg.setData(dataCh1, dataCh2);

fg.configChannel(1); % this will enable the channel
% figure();
% plot(dataCh1.data)
%%
timeToSample = N /sClk;
fClkVec = flip([5e6, 20e6]);
channelsVec = [1,2,4,8,16];
figure();
ind = 0;
numAvg = 100;
for i = 1:length(fClkVec)
    fClk = fClkVec(i);

    samplesInClk = floor(sClk/fClk);
    fClkReal     = sClk/samplesInClk;
    if mod(samplesInClk,2)
        samplesInDCHigh  = floor(0.5*samplesInClk)+1;
        samplesInDCLow   = samplesInClk - samplesInDCHigh;
    else
        samplesInDCHigh = 0.5*samplesInClk;
        samplesInDCLow = samplesInDCHigh;
    end
    fprintf("Requested sampling rate %.2f [MHz]\n", fClk/1e6);
    fprintf("Configured sampling rate %.2f [MHz]\n", fClkReal/1e6);

    dataCh2.data = repmat([ones(1,samplesInDCHigh), zeros(1,samplesInDCLow)], 1, 16); % 5Mhz clock
    dataCh2.dataLen = length(dataCh2.data);

%     figure();
%     subplot(1,2,1)
%     plot(dataCh1.data)
%     title("Signal")
%     subplot(1,2,2)
%     plot(dataCh2.data)
%     title("Clock")

    fg.setData(dataCh1, dataCh2);
    fg.configChannel(2); % this will enable the channel


    for j=1:length(channelsVec)
        ind = ind +1;
        channels = channelsVec(j);

        digiVars = Digitizer.uVarsCreate();

        fs = fClkReal;

        digiVars.mode      = 'TS'; 
        digiVars.fs        = fs;
        digiVars.useGPU    = false;
        digiVars.channels  = channels; 

        digiVars.triggerDelay = 0;
        digiVars.extClk  = false;
        digiVars.extTrig = true;

        digiVars.bufferSizeBytes = 1*(2^20);
        digiVars.samplesPerMeas  = ceil(timeToSample*fs);

        digiVars.exportCropped = false;
        digiVars.draw = true;

        digi.setVars(digiVars);
        digi.configure();
        for k=1:numAvg
            if ~mod(k,10)
                fprintf("%d/%d\n", k, numAvg)
            end
            tmpRes(:,:,k) = digi.acquire();
        end
        res{ind} = mean(tmpRes,3);
        clear('tmpRes');
        ax(i,j) = subplot(2,5,ind);
        plot( ax(i,j), res{ind}')
        drawnow();
        title(sprintf("%d Channels at %.2f[MHz]", digiVars.channels, fClkReal/1e6));
        linkaxes(ax(i,:));
    end
end

%% Monitor
% digi = Digitizer();
% digi.connect();
% 
% digiVars = Digitizer.uVarsMonitorCreate();
%  
% digiVars.fs        = 5e6;
% digiVars.channels  = 1; 
% 
% digiVars.triggerDelay = 0;
% digiVars.extClk  = true;
% digiVars.extTrig = true;
% 
% digiVars.timeToSample = 50e-6;
% digiVars.avgNum       = 1;
% 
% 
% digi.monitor(digiVars);

