function [status] = confExtClk(boardHandle, extClkRate)
%CONFEXTCLK Summary of this function goes here
%   Detailed explanation goes here
status = true;
AlazarDefs

% TODO: add logic to convert value to macro

% SAMPLE_RATE_USER_DEF
% SAMPLE_RATE_5MSPS
% EXTERNAL_CLOCK
% INTERNAL_CLOCK
retCode = ...
    AlazarSetCaptureClock(  ...
        boardHandle,        ... % HANDLE -- board handle
        EXTERNAL_CLOCK,  ... % U32 -- clock source id
        SAMPLE_RATE_USER_DEF,  ... % U32 -- sample rate id
        CLOCK_EDGE_RISING,  ... % U32 -- clock edge id
        0                   ... % U32 -- clock decimation
        );
if retCode ~= ApiSuccess
    fprintf('Error: AlazarSetCaptureClock failed -- %s\n', errorToText(retCode));
    status = false;
    return
end

end

