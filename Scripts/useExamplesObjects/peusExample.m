close all;
clear all;
clc;

addpath(genpath(pwd))

instrreset;

%% init the object
Scope_visa_address = 'USB0::0x1AB1::0x04CE::DS1ZA172115085::0::INSTR';
scope = Rigol_DS1074z('ni',Scope_visa_address);

stages = stages('PI');
stages.connect();

peus = pulseEchoUS(scope, stages);

%% Prepare the user variables
uVars = pulseEchoUS.uVarsCreate();

uVars.distFromCenter = 9.2e-2;
uVars.imageWidth     = 30e-3;
uVars.c              = 1440;
uVars.usFreq         = 1.25e6;

% Scope
uVars.signalChannel = 1;
uVars.trigCh        = 2;
uVars.avg           = 128;
uVars.usRepRate     = 100;
uVars.pulseIntVal   = 12;

% Grid & Stages
uVars.scanStart  = 10;
uVars.scanStride = 1;
uVars.scanEnd    = 15;

uVars.scanAxis   = 'X';  
uVars.depthAxis  = 'Y';
uVars.thirdAxis  = 'Z';

uVars.thirdPos   = 75;

uVars.moveSecondStage     = true;
uVars.stagesOwnedByParent = false;

% Graphics
uVars.figs.intExt        = 'int';
uVars.figs.useExtClims   = false;
uVars.figs.reopenFigures = false;

uVars.figs.depthAxType   = 'Normal';
uVars.figs.scanAxType    = 'Normal';

uVars.figs.validStruct.rawData = true;
uVars.figs.validStruct.recon   = true;

% FileSystem
uVars.fileSystem.saveResults  = true;
uVars.fileSystem.saveFigs     = false; 

uVars.fileSystem.dirPath     = "C:\Users\ahiad\OneDrive - Technion\Graduate\AcoustoOpticSystem\Measurements";
uVars.fileSystem.projName    = "fsTest";
uVars.fileSystem.resDirName  = "PEUSResults";

uVars.fileSystem.extProject        = false;
uVars.fileSystem.stackAllSubObjRes = false;
uVars.fileSystem.dontSaveVars      = false;

uVars.fileSystem.useExtVarsPath    = false;
uVars.fileSystem.extVarsPath       = [];

%% Work with the object:
peus.setVars(uVars);
vars = peus.getVars();

peus.configureScan();

res = peus.scan();

