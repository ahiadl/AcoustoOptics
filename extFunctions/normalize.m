function [ret, spanX, offset] = normalize(X)
%NORMALIZE Summary of this function goes here
%   Detailed explanation goes here
    spanX = max(X(:)) - min(X(:));
    offset = min(X(:));
    ret = ( X - offset )/ spanX;
end

