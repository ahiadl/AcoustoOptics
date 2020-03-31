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

%% Test on real data


%% New Algorithm
newData = reshape(this.data, this.digitizer.channels,...
                           this.samples.numOfQuant,...
                           this.samples.samplesPerQuant);


%Bring quant-dim to be first (averaging convinient)
newData = permute(newData, [2,1,3]);

newData = reshape(newData, this.samples.numOfQuant,...
                           this.digitizer.channels,...
                           this.samples.samplesPerPulse,...
                           this.samples.trainsPerQuant * this.samples.numOfPos);

newData = reshape(newData, this.samples.numOfQuant,...
                           this.digitizer.channels,... 
                           this.samples.samplesPerPulse,...
                           this.samples.numOfPos,this.samples.trainsPerQuant);

newData = permute(newData, [1,2,3,5,4]);

newData = reshape(newData, this.samples.numOfQuant,...
                           this.digitizer.channels,...
                           this.samples.samplesPerPos,...
                           this.samples.numOfPos);
 


%%Old Algorithm


oldData = reshape(this.data, this.digitizer.channels,...
                                  this.samples.numOfQuant,...
                                  this.samples.samplesPerQuant);

%Bring quant-dim to be first (averaging convinient)
oldData = permute(oldData, [2,1,3]);

%separate trains
oldData = reshape(oldData, this.samples.numOfQuant,...
                               this.digitizer.channels,...
                               this.samples.samplesPerTrain,...
                               this.samples.trainsPerQuant);

for i=1:this.samples.numOfPos
    reshapedOld(:,:,:,i) = reshape(oldData(:, :, ((i-1)*this.samples.samplesPerPulse+1) : (i*this.samples.samplesPerPulse), :),...
                 this.samples.numOfQuant,...
                 this.digitizer.channels,...
                 this.samples.samplesPerPos,...
                 1);
end

oldData = reshapedOld;