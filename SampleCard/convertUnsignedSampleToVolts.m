function [castedRawData] = convertUnsignedSampleToVolts(rawData, voltsRange )
   byteSize = 8;
   lowByteMask = uint16(31);
   shiftFact = 2;
   bitsRange = 2^16/2; 
   
%    rawData = bitshift(rawData,-2);
%    castingTime = tic;
   castedRawData =  cast(rawData, 'double');
%    castingTime = toc(castingTime)
   castedRawData = (castedRawData- bitsRange)*voltsRange / (bitsRange * shiftFact);
 end