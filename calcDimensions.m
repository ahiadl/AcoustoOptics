function [dim, timing] = calcDimensions(US, freqSample)

Tsine     = 1/US.freqSin;              %[s], duration of one period of the sin
Tpulse    = US.sinPeriods*US.Tsin;     %[s], duration of the pulse
Ttrain    = 1/US.freqTrain;            %[s], duration of one period of the train
Tsample   = 1/freqSample;              %[s], duration of a single sample

Npulse    = Tpulse/Tsample;          %[#], Amount of samples in one pulse duration
Ntrainraw = Ttrain/Tsample;          %[#], Amount of samples in one train cycle duration, raw because can be not a multiplier of Npulse
Npos      = floor(Ntrainraw/Npulse); %[#], Amount of positions we can yield in Z axis
Ntrain    = Npos*Npulse;             % this makes sure that Ntrain divided by p
Zres      = Npulse*Tsample*US.speed; %[m], possible resolution in Z axis 

dim.Npulse    = Npulse;
dim.Ntrainraw = Ntrainraw;
dim.Ntrain    = Ntrain;
dim.Npos      = pos;
dim.Zres      = Zres;

timing.Tsine   = Tsine;     
timing.Tpulse  = Tpulse;    
timing.Ttrain  = Ttrain;
timing.Tsample = Tsample;

end