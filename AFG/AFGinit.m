function [status] = AFGinit(AFG)
%AFGINIT - initiate the communication with the Arbitrary Function
%Generator and set it's static configurations.

%% Extract Struct variables to workspace
save('temp.mat','-struct','AFG');
load('temp.mat');
delete('temp.mat');
status = true;
%% Instrument Connection
% Find a GPIB object.
    obj1 = instrfind('Type', 'gpib', 'BoardIndex', 8, 'PrimaryAddress', 4, 'Tag', '');

% Create the GPIB object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = gpib('AGILENT', 8, 4);
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
fprintf(obj1, '*RST');

% Whoami
idn = query(obj1, '*IDN?');
fprintf('Instrument IDN: %s', idn);

%% Construct the Signals
% Convert data to 14 bit, as described in manual section 4-66
pulseData = uint16(((2^13-1)* pulseData)+2^13); % 14 bit data
% Convert data to byte array [low byte 0, high byte 0, low byte 1, ...]
bytes = zeros(2*samplesInTrain,1, 'uint8');
% Low byte
bytes(1:2:2*length(pulseData)) = uint8(bitand(pulseData, uint16(hex2dec('00FF'))));
% High byte
bytes(2:2:2*length(pulseData)) = uint8(bitshift(bitand(pulseData, uint16(hex2dec('FF00'))), -8));
err = zeros(1,4);
%% Send the Pulse to the AFG

while ~strcmp(err(1:4), '-100')
    % ------ Configure the Output ----------
    fprintf(obj1, ':INSTrument 1');
    fprintf(obj1, ':OUTP OFF');
    % Set a filter of 25Mhz at the output
    fprintf(obj1, [':OUTP:FILT ', daqFilter, 'M']);
    % Set a sync signal of a 16 samples bit
    fprintf(obj1, ':OUTP:SYNC:SOUR BIT');
    fprintf(obj1, [':OUTP:SYNC:WIDT ', triggertWidth]);
    disp(query(obj1, ':OUTP:SYNC:WIDT?'))
    fprintf(obj1, ':OUTP:SYNC ON');

    % ----- Define the Samples CLK ---------
    % Set arbitrary waveform SCLK to 2Gs/s
    fprintf(obj1, [':FREQ:RAST ', num2str(Sclk, '%10.4e')]);
    disp(query(obj1, ':FREQ:RAST?'))
    
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
    fprintf(obj1, [':TRAC:DEF 1,' num2str(samplesInTrain)]);
    cLength = strcat(num2str(length(num2str(2*samplesInTrain))),num2str(2*samplesInTrain));
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
    err = query(obj1, ':SYST:ERR?');
    % Wait for completion
    query(obj1, '*OPC?');
% Enable
end 
fprintf(obj1, ':ENAB');
%%
% ------------Create The Extrernak Clock Signal------------------------------------------------
% Number of waveform samples - shift the traing frequency so will have a
% number of samples divided by 16.
% dotsPerCycle = ceil(Sclk/extClkFreq);
% minCycInSig = (floor(192/dotsPerCycle)+1)*16; %this is to make sure signal is longer than 192 samples and can be divided by 16.
% extClkSamples = dotsPerCycle*minCycInSig;
% Convert data to 14 bit, as described in manual section 4-66
clkData = uint16(((2^13-1)*clkData)+2^13); % 14 bit data
% Convert data to byte array [low byte 0, high byte 0, low byte 1, ...]
bytes = zeros(2*extClkSamples,1, 'uint8');
% Low byte
bytes(1:2:2*length(clkData)) = uint8(bitand(clkData, uint16(hex2dec('00FF'))));
% High byte
bytes(2:2:2*length(clkData)) = uint8(bitshift(bitand(clkData, uint16(hex2dec('FF00'))), -8));

%%
err = zeros(1:4);
while ~strcmp(err(1:4), '-100')
    % % ------ Configure the Output ----------
    fprintf(obj1, ':INSTrument 2');
    fprintf(obj1, ':OUTP OFF');
    % Set a filter of 25Mhz at the output
    fprintf(obj1, ':OUTP:FILT 50M');

    % ----- Define the Samples CLK ---------
    % Set arbitrary waveform SCLK to 2Gs/s
    fprintf(obj1, [':FREQ:RAST ', num2str(Sclk, '%10.4e')]);
    % Select internal source for the SCLK
    fprintf(obj1, ':FREQ:RAST:SOUR INT');

    % ----- Set Signal Properties ---------
    % Set signal source to as given by user
    fprintf(obj1, ':FUNC:MODE USER');
    % Setting the amplitude and bias
    fprintf(obj1, [':VOLT:LEV:AMPL ', extClkAmp] );
    fprintf(obj1, [':VOLT:LEV:OFFS ', '0']);
    % Set continous mode on
    fprintf(obj1, ':INIT:CONT 1');
    % Select channel 1 as the active channel
    % fprintf(obj1, ':INST:SEL CH2');

    % ------ Data Source ---------
    % Delete all segments from memory
    fprintf(obj1, ':TRAC:DEL:ALL');
    % Define trace
    fprintf(obj1, [':TRAC:DEF 1,' num2str(extClkSamples)]);
    cLength = strcat(num2str(length(num2str(2*extClkSamples))),num2str(2*extClkSamples));
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
    err =query(obj1, ':SYST:ERR?');
    % Wait for completion
    query(obj1, '*OPC?');
    % Enable
    fprintf(obj1, ':ENAB');
end
end
