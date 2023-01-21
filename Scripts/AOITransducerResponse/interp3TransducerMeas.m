%% Load Data
dir = "C:/Users\sahiadl.EED/OneDrive - Technion/Graduate/AcoustoOpticSystem/Measurements/Transducer Pressure Field/Calibrated/Focused";
resDepth = load(sprintf("%s/%s", dir, "1.25MHz-Depth-Axis"));
resTrans = load(sprintf("%s/%s", dir, "1.25MHz-Transversal-Axis"));
resPulse = load(sprintf("%s/%s", dir, "1.25MHz-FocalPoint"));
resHad   = load(sprintf("%s/%s", dir, "1.25MHz-FocalPoint-HadamardSqnc"));

csVars = resTrans.csVars;
csVars.axDisc2Span   = resDepth.csVars.axScanSpan;
csVars.axDisc2Span   = resDepth.csVars.axScanSpan;
csVars.axDisc2Stride = resDepth.csVars.axScanStride;
csVars.axDisc2Ref    = resDepth.csVars.axDisc2Ref;
csVars.disc2Vec      = resDepth.csVars.scanVec;
csVars.disc2VecBin   = resDepth.csVars.scanVecBin;


close all

% resDepth = 
%%

