function [freqTrainTrue, status] = AFGinit(AFG)
%AFGINIT - initiate the communication with the Arbitrary Function
%Generator and set it's static configurations.

%% Extract Struct variables to workspace
save('temp.mat','-struct','AFG');
load('temp.mat');
delete('temp.mat');
freqTrainTrue = 0;
status = 0;
%% Instrument Connection
% Find a GPIB object.
    obj1 = instrfind('Type', 'gpib', 'BoardIndex', 7, 'PrimaryAddress', 4, 'Tag', '');

% Create the GPIB object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = gpib('AGILENT', 7, 4);
else
    fclose(obj1);
    obj1 = obj1(1);
end
blocksize = 10000;
set(obj1, 'OutputBufferSize', blocksize);

% Connect to instrument object, obj1.
fopen(obj1);
%% Configure Instrument
% Clearing all the old errors from error's buffer.
fprintf(obj1, '*CLS');

% Reseting the device.
% fprintf(obj1, '*RST');

% Whoami
idn = query(obj1, '*IDN?');
fprintf('Instrument IDN: %s', idn);

% fprintf(obj1, ':OUTP ON');
% fprintf(obj1, ':OUTP:FILT: 25M');

%% Construct the Pulse
% ------------Create 500KHz sine wave for testing------------------------------------------------
% Number of waveform samples - shift the traing frequency so will have a
% number of samples divided by 16.
numOfSamples = max(192, ceil((Sclk/freqTrain)/16)*16);
freqTrainTrue = Sclk/numOfSamples;

% Setting the output frequency
dt = 1/Sclk;
t = 0:dt:(sinPeriods/freqPulse)- dt;

data = zeros(numOfSamples, 1);
data(1:length(t)) = sin(2*pi*freqPulse*t);

if (AFG.debug)
    figure()
    stem(data)
    title('Signal Downloaded to AFG');
end

% Convert data to 14 bit, as described in manual section 4-66
data = uint16(((2^13-1)*data)+2^13); % 14 bit data
% Convert data to byte array [low byte 0, high byte 0, low byte 1, ...]
bytes = zeros(2*numOfSamples,1, 'uint8');
% Low byte
bytes(1:2:2*length(data)) = uint8(bitand(data, uint16(hex2dec('00FF'))));
% High byte
bytes(2:2:2*length(data)) = uint8(bitshift(bitand(data, uint16(hex2dec('FF00'))), -8));


%% Send the Pulse to the AFG

% ------ Configure the Output ----------
% Set a filter of 25Mhz at the output
fprintf(obj1, [':OUTP:FILT ', daqFilter, 'M']);
% Set a sync signal of a 16 samples bit
fprintf(obj1, ':OUTP:SYNC:SOUR BIT');
fprintf(obj1, [':OUTP:SYNC:WIDT ', triggertWidth]);
fprintf(obj1, ':OUTP:SYNC ON');

% ----- Define the Samples CLK ---------
% Set arbitrary waveform SCLK to 2Gs/s
fprintf(obj1, [':FREQ:RAST ', num2str(Sclk, '%10.4e')]);
% Select internal source for the SCLK
fprintf(obj1, ':FREQ:RAST:SOUR INT');

% ----- Set Signal Properties ---------
% Set signal source to as given by user
fprintf(obj1, ':FUNC:MODE USER');
% Setting the amplitude and bias
fprintf(obj1, [':VOLT:LEV:AMPL ', amp] );
fprintf(obj1, [':VOLT:LEV:OFFS ', bias]);
% Set continous mode on
fprintf(obj1, ':INIT:CONT 1');
% Select channel 1 as the active channel
fprintf(obj1, ':INST:SEL CH1');

% ------ Data Source ---------
% Delete all segments from memory
fprintf(obj1, ':TRAC:DEL:ALL');
% Define trace
fprintf(obj1, [':TRAC:DEF 1,' num2str(numOfSamples)]);
cLength = strcat(num2str(length(num2str(2*numOfSamples))),num2str(2*numOfSamples));
samples = sprintf('#%s',cLength);
% Select trace
fprintf(obj1, ':TRAC:SEL 1');
% Write data blockwise
fprintf(obj1, '%s', [':TRAC:DATA' samples]);
n = numel(bytes);
num_blocks = ceil(n/blocksize);

for k=1:num_blocks
pause(0.2)
fprintf('downloading chunk %d\n', k);
fwrite(obj1, bytes((k-1)*blocksize+1:min(n,k*blocksize)), 'uint8');
end

% Connect the Output
fprintf(obj1, ':OUTP ON');
% Checking for Errors
query(obj1, ':SYST:ERR?')
% Wait for completion
query(obj1, '*OPC?');
% Enable
fprintf(obj1, ':ENAB');

end
