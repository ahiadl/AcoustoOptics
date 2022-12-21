close all
clear all
clc

pd = PhantomDesigner();

uVars.wl   = 785;
uVars.Vph  = 200;
uVars.musP = 0.01*2.^(1:5); %[mm^-1]
uVars.mua  = 0.0023 * (2.^[1:5]); %[mm^-1]
uVars.numLayers = 1;
uVars.numPh = 5;
uVars.inkType = 'RoyalTalens';
uVars.layerThick = 56; %[mm]

pd.createPhantomSet(uVars);

%% Graded Absorption:
clc;
uVars.musP = 0.4; %[mm^-1]
uVars.mua  = 2.3e-3 * 2.^(1:5); %[mm^-1]
pd.gradedAbsorption(uVars);

%% Graded Scattering:
clc;
uVars.musP = 0.01*2.^(1:5); %[mm^-1]
uVars.mua  = 2*4.6e-3;
pd.gradedScattering(uVars);