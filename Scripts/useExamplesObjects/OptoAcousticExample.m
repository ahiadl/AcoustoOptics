close all
clear all
clc;
%%
stages = stages('PI');
stages.connect();
% stages.moveAbsAx('Z', 173);
%%
oa = optoAcoustic();

uVars = oa.uVarsCreate();

uVars.activeChannels  = 1:256; %Mask
uVars.numOfSamples    = 2030;
uVars.numOfAvg        = 10;
uVars.numOfFrames     = 1;

uVars.c = 1430;
uVars.fs = 40e6;
uVars.geometry = 'Circular';
uVars.imageWidth = 50e-3; %[m]
uVars.BPmode = 3;

uVars.scanStart  = 155;
uVars.scanStride = 1;
uVars.scanEnd    = 170;

uVars.firstAxis  = 'X';
uVars.secondAxis = 'Y';
uVars.scanAxis   = 'Z';

uVars.figs.firstAxType  = 'Normal';
uVars.figs.secondAxType = 'Normal';
uVars.figs.scanAxType   = 'Normal';
uVars.figs.timeAxType   = 'Normal';

uVars.figs.scanAxCoor   = 218;

% uVars.figs.displayAxAsIdx = false;
uVars.figs.useExtClims = false;
uVars.figs.reopenFigures = false;
uVars.figs.intExt = 'int';

uVars.figs.validStruct.sinogram = true;
uVars.figs.validStruct.recon    = true;
uVars.figs.validStruct.navSinogram  = true;
uVars.figs.validStruct.navRecon     = true;
           
uVars.fileSystem.saveResults  = false;
uVars.fileSystem.saveFigs     = false; 

uVars.fileSystem.dirPath     = "C:\Users\sahiadl.EED\OneDrive - Technion\Graduate\AcoustoOpticSystem\Measurements";
uVars.fileSystem.projName    = "fsTest";
uVars.fileSystem.resDirName  = "OAResults";

uVars.fileSystem.extProject        = false;
uVars.fileSystem.stackAllSubObjRes = false;
uVars.fileSystem.dontSaveVars      = false;

uVars.fileSystem.useExtVarsPath    = false;
uVars.fileSystem.extVarsPath       = [];

oa.setVars(uVars);
vars = oa.getVars();
oa.configPeripherals();
%%
%%
res = oa.runLiveOA();
% res = oa.scan();
%%
zVec = 207:1:226;

stages.moveAbsAx('Z', 240);
resTmp = oa.runOA();

for i=1:length(zVec)
    stages.moveAbsAx('Z', zVec(i));
    res(i) = oa.runOA();
end

figure()
for i=1:length(zVec)
    imagesc(res(i).recon);
    title(sprintf("Z=%d", zVec(i)));
    pause
end

figure()
imagesc(resRough.res.recon(:,:,20))
