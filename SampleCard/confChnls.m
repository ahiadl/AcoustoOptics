function [channelMask, status] = confChnls(boardHandle, numOfChannels)

AlazarDefs;
status = true;
inputRange = INPUT_RANGE_PM_1_V;

retCode = AlazarInputControlEx( boardHandle, CHANNEL_A, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_B, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_C, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_D, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_E, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_F, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_G, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_H, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_I, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_J, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_K, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_L, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_M, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_N, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_O, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end
retCode = AlazarInputControlEx( boardHandle, CHANNEL_P, DC_COUPLING, inputRange, IMPEDANCE_50_OHM);
if (retCode ~= ApiSuccess); fprintf('Error: AlazarInputControlEx failed -- %s\n', errorToText(retCode)); status = false; return; end

% TODO: verify this mask -Verified!
channelMask = 0;
for i = 1:numOfChannels
    channelMask = bitor(bitshift(channelMask,1),1);
end

end

