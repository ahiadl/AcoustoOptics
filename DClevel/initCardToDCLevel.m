function [boardHandle, bytesPerSample, status] = initCardToDCLevel()
%--------------------------------
%--- Init Variables and board ---
%--------------------------------

AlazarDefs;

if ~alazarLoadLibrary(); fprintf('Error: ATSApi library not loaded\n'); status = false; return; end

systemId = int32(1);
boardId = int32(1);

boardHandle = AlazarGetBoardBySystemID(systemId, boardId);
setdatatype(boardHandle, 'voidPtr', 1, 1);

%----------------------------
%--- Configure Clk & Trig ---
%----------------------------

retCode = ...
    AlazarSetCaptureClock(  ...
        boardHandle,        ... % HANDLE -- board handle
        INTERNAL_CLOCK,     ... % U32 -- clock source id
        SAMPLE_RATE_100MSPS, ... % U32 -- sample rate id
        CLOCK_EDGE_RISING,  ... % U32 -- clock edge id
        0                   ... % U32 -- clock decimation
        );
if retCode ~= ApiSuccess
    fprintf('Error: AlazarSetCaptureClock failed -- %s\n', errorToText(retCode));
    return
end

[channelMask, status] = confChnls(boardHandle, 4);
if(~status); return; end

retCode = ...
    AlazarSetTriggerOperation( ...
        boardHandle,        ... % HANDLE -- board handle
        TRIG_ENGINE_OP_J,   ... % U32 -- trigger operation
        TRIG_ENGINE_J,      ... % U32 -- trigger engine id
        TRIG_CHAN_A,        ... % U32 -- trigger source id
        TRIGGER_SLOPE_POSITIVE, ... % U32 -- trigger slope id
        150,                ... % U32 -- trigger level from 0 (-range) to 255 (+range)
        TRIG_ENGINE_K,      ... % U32 -- trigger engine id
        TRIG_DISABLE,       ... % U32 -- trigger source id for engine K
        TRIGGER_SLOPE_POSITIVE, ... % U32 -- trigger slope id
        128                 ... % U32 -- trigger level from 0 (-range) to 255 (+range)
        );
if retCode ~= ApiSuccess
    fprintf('Error: AlazarSetTriggerOperation failed -- %s\n', errorToText(retCode));
    return
end

retCode = ...
    AlazarSetExternalTrigger( ...
        boardHandle,        ... % HANDLE -- board handle
        DC_COUPLING,        ... % U32 -- external trigger coupling id
        ETR_TTL              ... % U32 -- external trigger range id
        );
if retCode ~= ApiSuccess
    fprintf('Error: AlazarSetExternalTrigger failed -- %s\n', errorToText(retCode));
    return
end

triggerDelay_sec = 0;
triggerDelay_samples = 0;
retCode = AlazarSetTriggerDelay(boardHandle, triggerDelay_samples);
if retCode ~= ApiSuccess
    fprintf('Error: AlazarSetTriggerDelay failed -- %s\n', errorToText(retCode));
    return;
end

retCode = ...
    AlazarSetTriggerTimeOut(    ...
        boardHandle,            ... % HANDLE -- board handle
        1   ... % U32 -- timeout_sec / 10.e-6 (0 == wait forever)
        );
if retCode ~= ApiSuccess
    fprintf('Error: AlazarSetTriggerTimeOut failed -- %s\n', errorToText(retCode));
    return
end

AlazarSetParameter( ...
boardHandle, ...% HANDLE -- board handle
0, ...% U8 -- channel Id (not used)
SET_DATA_FORMAT, ...% U32 -- parameter to set
DATA_FORMAT_UNSIGNED ...% long -- value (0 = unsigned, 1 = signed)
);

[retCode, boardHandle, ~, bitsPerSample] = AlazarGetChannelInfo(boardHandle, 0, 0);
if retCode ~= ApiSuccess; fprintf('Error: AlazarGetChannelInfo failed -- %s\n', errorToText(retCode)); status = false; return; end

bytesPerSample = floor((double(bitsPerSample) + 7) / double(8));

end

