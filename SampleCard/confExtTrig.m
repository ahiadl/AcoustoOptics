function [status] = confExtTrig(boardHandle, samplesPerSec)

AlazarDefs;
status = true;
triggerLevel_volts = 3; % trigger level
triggerRange_volts = 5; % input range
triggerLevel_code =(128 + 127 * triggerLevel_volts / triggerRange_volts);

% Disable internal Triggers and set External Trigger
retCode = ...
    AlazarSetTriggerOperation( ...
        boardHandle,            ... % HANDLE -- board handle
        TRIG_ENGINE_OP_J,       ... % U32 -- trigger operation
        TRIG_ENGINE_J,          ... % U32 -- trigger engine id
        TRIG_EXTERNAL,          ... % U32 -- trigger source id
        TRIGGER_SLOPE_POSITIVE, ... % U32 -- trigger slope id
        triggerLevel_code,      ... % U32 -- trigger level from 0 (-range) to 255 (+range)
        TRIG_ENGINE_K,          ... % U32 -- trigger engine id
        TRIG_DISABLE,           ... % U32 -- trigger source id for engine K
        TRIGGER_SLOPE_POSITIVE, ... % U32 -- trigger slope id
        128                     ... % U32 -- trigger level from 0 (-range) to 255 (+range)
        );
if retCode ~= ApiSuccess
    fprintf('Error: AlazarSetTriggerOperation failed -- %s\n', errorToText(retCode));
    status = false;
    return
end

% Set External trigger Configurations
retCode = ...
    AlazarSetExternalTrigger( ...
        boardHandle,        ... % HANDLE -- board handle
        DC_COUPLING,        ... % U32 -- external trigger coupling id
        ETR_TTL             ... % U32 -- external trigger range id
        );
if retCode ~= ApiSuccess
    fprintf('Error: AlazarSetExternalTrigger failed -- %s\n', errorToText(retCode));
    status = false;
    return
end
    
% Trigger Delay set to 0
triggerDelay_sec = 0;
triggerDelay_samples = uint32(floor(triggerDelay_sec * samplesPerSec + 0.5));
retCode = AlazarSetTriggerDelay(boardHandle, triggerDelay_samples);
if retCode ~= ApiSuccess
    fprintf('Error: AlazarSetTriggerDelay failed -- %s\n', errorToText(retCode));
    status = false;
    return;
end

% Set Trigger Timeout 
triggerTimeout_sec = 0;
triggerTimeout_clocks = uint32(floor(triggerTimeout_sec / 10.e-6 + 0.5));
retCode = ...
    AlazarSetTriggerTimeOut(    ...
        boardHandle,            ... % HANDLE -- board handle
        triggerTimeout_clocks   ... % U32 -- timeout_sec / 10.e-6 (0 == wait forever)
        );
if retCode ~= ApiSuccess
    fprintf('Error: AlazarSetTriggerTimeOut failed -- %s\n', errorToText(retCode));
    status = false;
    return
end    
    
end

