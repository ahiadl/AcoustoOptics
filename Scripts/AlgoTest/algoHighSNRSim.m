clear all
close all
clc

numOfCh      = 1;
numOfPos     = 2;
numOfQuant   = 2;
trainPerQuant   = 2;
samplesPerTrain = 4;

for i=1:numOfCh
    for j = 1:numOfPos
        for k = 1:numOfQuant
            for m = 1:trainPerQuant
                for n = 1:samplesPerTrain
                    sig(i,(j-1)*numOfQuant*trainPerQuant*samplesPerTrain+(k-1)*trainPerQuant*samplesPerTrain+(m-1)*samplesPerTrain+n) = sprintf("C%dP%dQ%dT%dS%d", i,j,k,m,n);
                end
            end
        end
    end
end

totalSamples =numOfCh*numOfPos;



A = sparse(eye(samplesPerTrain););
B = sparse(eye(samplesPerTrain));
C = kron(A,B)

%% 

ao = acoustoOptics();

