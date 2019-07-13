clear A 

ch = 5;
numOfTrains = 2;
samplesPerPulse = 3;
numOfPos = 4;
samplesPerTrain = numOfPos*samplesPerPulse;
samplesPerSignal = samplesPerTrain * numOfTrains;
samplesPerPos = numOfTrains * samplesPerPulse;
chanelCode = ["A", "B", "C", "D", "E"];
for i = 1:ch
    cc = chanelCode(i);
    curCode = [];
    for j = 1:numOfTrains
        for k = 1:numOfPos
            for p = 1:samplesPerPulse
                curCode = [curCode, sprintf("%s%s%s%s", cc, num2str(j), num2str(k), num2str(p))];
            end
        end
    end
    A(i,:) = curCode;
end

        
% A = [ "A11", "A12", "A13", "A14", "A21", "A22", "A23", "A24", "A31", "A32", "A33", "A34", "411", "A42", "A43", "A44", "A51", "A52", "A53", "A54";...
%       "B11", "B12", "B13", "B14", "B21", "B22", "B23", "B24", "B31", "B32", "B33", "B34", "B11", "B42", "B43", "B44", "B51", "B52", "B53", "B54"];
  
B = A';

C = reshape(A, ch, samplesPerTrain , numOfTrains);

D = permute(C, [1, 3, 2]);

E = reshape(D, ch, samplesPerPos, numOfPos);