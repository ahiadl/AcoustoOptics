close all;
clear all;
clc;

%% Load Results:

for i=1:8
    filename = sprintf("D:/DCLevel/DCLvl-%d/AO-Results.mat", i-1);
    res(i) = load(filename);
    rawPhi(i,:) = res(i).rawPhi;
    analysis(i) = res(i).analysis.rawPhi;
end

intVec = 1.1*(0.75.^(0:7));


%% 
rawPhiCut = rawPhi(:, 25:99);

avgBkg = mean(rawPhiCut(:, 1:40), 2);
stdBkg = std(rawPhiCut(:, 1:40), 0 , 2);

peak = max(rawPhiCut, [], 2);

signal = peak-avgBkg;


figure(); 
plot(rawPhiCut')
xlabel("Idx")
ylabel("Fluence[AU]")
for i=1:8; str(i) = sprintf("Power = %.2f", intVec(i)); end
legend (str, 'Location', 'northwest')

figure();
plot(intVec, signal, '-+')
xlabel("Power [W]")
ylabel("Signal [AU]")



%%
%--------------------
% Diverged Reflection
%--------------------

%% Load Results:

for i=1:6
    filename = sprintf("D:/DCLevel/divergedReflection/DCLvl-%d/AO-Results.mat", i-1);
    res(i) = load(filename);
    rawPhi(i,:) = res(i).rawPhi;
    analysis(i) = res(i).analysis.rawPhi;
end

intVec = 1.1*(0.75.^(0:5));


%% 
rawPhiCut = rawPhi(:, 25:99);

avgBkg = mean(rawPhiCut(:, 1:40), 2);
stdBkg = std(rawPhiCut(:, 1:40), 0 , 2);

peak = max(rawPhiCut, [], 2);

signal = peak-avgBkg;


figure(); 
plot(rawPhiCut')
xlabel("Idx")
ylabel("Fluence[AU]")
for i=1:6; str(i) = sprintf("Power = %.2f", intVec(i)); end
legend (str, 'Location', 'northwest')

figure();
plot(intVec, signal, '-+')
xlabel("Power [W]")
ylabel("Signal [AU]")

%%
%--------------------
% Ultrasound Beam
%--------------------

%% Load Results:

for i=1:8
    filename = sprintf("D:/USLevel/USLvl-%d/AO-Results.mat", i-1);
    resUS(i) = load(filename);
    rawPhiUS(i,:) = resUS(i).rawPhi;
    analysisUS(i) = resUS(i).analysis.rawPhi;
end

intVecUS = 100*(0.75.^(0:7));


%% 
rawPhiUSCut = rawPhiUS(:, 25:99);

avgBkgUS = mean(rawPhiUSCut(:, 1:40), 2);
stdBkgUS = std(rawPhiUSCut(:, 1:40), 0 , 2);

peakUS = max(rawPhiUSCut, [], 2);

signalUS = peakUS-avgBkgUS;

figure(); 
plot(rawPhiUSCut')
xlabel("Idx")
ylabel("Fluence[AU]")
for i=1:8; str(i) = sprintf("Power = %.0f[%%]", intVecUS(i)); end
legend (str, 'Location', 'northwest')

figure();
plot(intVecUS, signalUS, '-+')
xlabel("Power [%]")
ylabel("Signal [AU]")