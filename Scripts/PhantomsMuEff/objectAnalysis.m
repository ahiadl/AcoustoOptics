close all;
clear all;
clc;

muObj = muEffAnalysis();


simFilename = "C:\Users\sahiadl.EED\OneDrive - Technion\Graduate\AcoustoOpticSystem\Measurements\MuEff\2Layers-40x30x10-DxWxI-mm-1L-1-FullResults.mat";

muObj.loadSingleLayerSim(simFilename, 1);

