close all;
clear all;
clc;

%-------------------------------
% Flow Parameters
%-------------------------------

% Positioning Paramteres:
userParams.Stages.port = 'COM3';
userParams.Stages.startPosX = 65; %[mm]
userParams.Stages.startPosY = 95; %[mm]
userParams.Stages.endPosX   = 123; 
userParams.Stages.endPosX   = 28;

% Ultrasound parameters
userParams.US.ferqSin       = 1.25e6;  %[Hz]
userParams.US.freqTrain     = 19e3;    %[Hz]
userParams.US.freqTrueTrain = 0; %[Hz], init value only, true value will be recieived from AFGinit.
userParams.US.speed         = 1550;   %[m/s]
userParams.US.sinPeriods    = 4;

% Function Generator Parameters
userParams.AFG.freqSin       = userParams.US.ferqSin;
userParams.AFG.numOfPeriods  = 4;
userParams.AFG.freqTrain     = userParams.US.freqTrain;
userParams.AFG.Sclk          = 100e6; %[Hz]
userParams.AFG.triggertWidth = '500'; %[Samples]
userParams.AFG.amp           = '1'; %[V]
userParams.AFG.bias          = '0'; %[V]
userParams.AFG.daqFilter     = '25'; %[Mhz], discrete values [25,50]

% Laser Parameters
userParams.Laser.waveLen = 785e-9; %[m]

% Sampling Parameters
userParams.Sample.freqSample = 100e6; %[S/s]
userParams.Sample.periodsToSample =  8192; %[#]

%DEBUG
userParams.AFG.debug    = 0;
userParams.Stages.debug = 0;


%-------------------------------
% Init the Instruments
%-------------------------------

status = 0;
% releaseing all connections to innstruments from outside matlab
instrreset 

%------ Arbitrary Function generator -------
fprintf("Initiating the Arbitrary Function Generator");
[params.US, status] = AFGinit(userParams.AFG);
if status
    printf("Arbitrary Function Generator is no valid") 
    return
end
fprintf("----------Done INIT AFG------------");

%-------------- Zaber Stages ----------------
fprintf("Starting The Stages Control");
[params.Stages, status] = StagesInit(userParams.Stages);
if status
    printf("Stages init is not valid") 
    return
end
fprintf("----------Done INIT Stages------------");

%------------- I/O controller ---------------

%------------- Sampling Card ----------------


%-------------------------------
% Scan, Sample and Analyze
%-------------------------------

fprintf("Calculating Problem Dimensions according to User Params");
[params.dimensions, params.timing] = calcDimensions(params.US, params.Sample.freqSample);
fprintf("----------Done Calculate Dimensions------------");

fprintf("Starting Scanning and Sampling");
for Xidx = 1:lengthX
    for Yidx = 1:lengthY
        % move to point
        % start puls train and sampling - I/O and Sampling cards
        % stop pulse and sampling - I/O and Sampling cards
        % yield results - Sampling Cards
        U = reshapeRawData(rawData, params.dimensions);          %reshape data
        I(Yidx,Xidx,:) = signalProcessing(U, params.dimensions); %fft for each Z pos on this (X,Y) pos
    end
end








