function [] = initFlowNew(userParams, initiator, app)
% addpath(genpath([pwd,'\..\']]));
if (initiator ~= initType.initByGUI)
    app = [];
%     addpath(genpath(pwd));
    %-------------------------------
    % Flow Parameters
    %-------------------------------
    % Positioning Paramteres:
    userParams.Stages.port = 'COM3';
    userParams.Stages.startPosX = 35; %[mm]
    userParams.Stages.startPosY = 95; %[mm]
    userParams.Stages.endPosX   = 95; 
    userParams.Stages.endPosY   = 28;
    userParams.Stages.Xstride   = 0.5; %[mm]
    userParams.Stages.Ystride   = 0.5; %[mm]

    % Ultrasound parameters
    userParams.US.freqSin       = 1.25e6;  %[Hz]
    userParams.US.freqTrain     = 19e3;    %[Hz]
    userParams.US.freqTrueTrain = 0;       %[Hz], init value only, true value will be recieived from AFGinit.
    userParams.US.speed         = 1550;    %[m/s]
    userParams.US.sinPeriods    = 4;

    % Function Generator Parameters
    userParams.AFG.freqSin       = userParams.US.freqSin;
    userParams.AFG.sinPeriods    = userParams.US.sinPeriods;
    userParams.AFG.freqTrain     = userParams.US.freqTrain;
    userParams.AFG.Sclk          = 100e6; %[Hz]
    userParams.AFG.triggertWidth = '16';  %[Samples]
    userParams.AFG.amp           = '1';   %[V]
    userParams.AFG.bias          = '0';   %[V]
    userParams.AFG.daqFilter     = '25';  %[Mhz], discrete values [25,50]
    userParams.AFG.extClkFreq    = 5e6;   %[S/s] 
    userParams.AFG.extClkAmp     = '2';   %[V]
    userParams.AFG.extClkOff     = '-1';   %[V]
    userParams.AFG.extClkDcyc    = 50;    %[%]

    % Laser Parameters
    userParams.Laser.waveLen = 785e-9; %[m]

    % Sampling Parameters
    userParams.Sample.channels              =  4;
    userParams.Sample.freqSample            =  userParams.AFG.extClkFreq; %[S/s]
    userParams.Sample.periodsToSample       =  8192;  %[#]
    userParams.Sample.trueSampleFreq        =  0;     %[#], init value
    userParams.Sample.recordsPerBuffer      =  4096;  %[#]
    userParams.Sample.buffersPerAcquisition =  2;     %[#]
    userParams.Sample.voltsRange            =  2;     %[volts]
    
    %DEBUG
    userParams.AFG.debug    = 1;
    userParams.Stages.debug = 0;
    
    %General Properties
    userParams.General.path = [ pwd, '/Results/'];
    userParams.General.filename = strrep(['Results-',datestr(datetime('now'))], ':', '-');
    userParams.General.useGPU = false;
    userParams.General.plotFFT = true;
    userParams.General.plotSpatialResult = true;
    userParams.General.fftPlotHandle = figure();
    userParams.General.spatialPowerSpectrumHandle = figure();
    userParams.General.Zcoor = 8;
    userParams.General.channelToPlot = 1;
    userParams.General.initiator = initiator;
end
% Verify Number of Channels
if (userParams.Sample.channels <1 || userParams.Sample.channels >16)
    fprintf('Error: Invalid channel mask %08X\n', channelMask);
    return
end

[Results.I, Results.params, status] = fullMeasMode(userParams, app);

if (status)
    save([userParams.General.path,userParams.General.filename],'-struct', 'Results');
%     save([userParams.General.path,userParams.General.filename, '-params'], 'params');
end

end

