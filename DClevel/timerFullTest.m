close all;
clear all;
clc;
addpath(genpath(pwd));
[app.boardHandle, app.bytesPerSample, ~] = initCardToDCLevel(); 
% app.pBufferA = libpointer('uint16Ptr',zeros(5008,1));
% app.pBufferB = libpointer('uint16Ptr',zeros(5008,1));
% app.pBufferC = libpointer('uint16Ptr',zeros(5008,1));
% app.pBufferD = libpointer('uint16Ptr',zeros(5008,1));

voltsRange = 2; %[V]
shiftFact  = 4; %[2 bits]
bitsRange  = 2^14 / 2;
for i = 1:100000
    tic
    rawData = acquireDCLevelCS(app);
%     samples = [app.pBufferA.Value, app.pBufferB.Value, app.pBufferC.Value, app.pBufferD.Value];
    data = convertUnsignedSampleToVolts(rawData, voltsRange);

    % Calc Rms
    RMS = sqrt(mean(data.^2 , 1));
    
    Value = sum(RMS) ;
    toc
    disp(RMS)
    pause(0.5)

end

clear ('app.pBufferA')
clear ('app.pBufferB')
clear ('app.pBufferC')
clear ('app.pBufferD')
