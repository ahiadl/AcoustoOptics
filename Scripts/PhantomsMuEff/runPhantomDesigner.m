close all
clear all
clc

pd = PhantomDesigner();

uVars.wl   = 785;
uVars.Vph  = 200;
uVars.musP = 0.01*2.^(1:5); %[mm^-1]
uVars.mua  = 2*[4.6e-3, 4.6e-3, 4.6e-3, 4.6e-3, 4.6e-3]  ; %[mm^-1]
uVars.inkType = 'RoyalTalens';
uVars.layerThick = 56; %[mm]

pd.createPhantomSet(uVars)