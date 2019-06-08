function [AFG, Stages, IOdev, digitizer, Sample] = peripheralsInitFlow(vars)
%------ Arbitrary Function generator -------
status = true;
fprintf("Initiating the Arbitrary Function Generator\n");

[AFG] = calcTrueAFGParams(vars.AFG);
if (~vars.Debug.dontActivateAFG)
    [status] = AFGinit(AFG);
end
US = vars.US;
US.freqTrueTrain = AFG.trainTrueFreq; 
vars.Sample.trueSampleFreq = AFG.extClkTrueFreq;
if (~status); printf("Arbitrary Function Generator is no valid\n"); return; end
fprintf("----------Done INIT AFG------------\n");

%-------------- Zaber Stages ----------------
fprintf("Starting The Stages Control\n");
[Stages, status] = StagesInit(vars.Stages);
if ~status; printf("Stages init is not valid\n"); return; end
fprintf("----------Done INIT Stages------------\n");

%------------- I/O controller ---------------
fprintf("Starting The I/O controller\n");
if ~vars.Debug.dontOpenIO; 
    IOdev = initIOControl(); closeIO(IOdev); 
else
    IOdev = [];
end
fprintf("Done INIT I/O controller\n");

%------------- Dimensions ---------------
fprintf("Calculating Problem Dimensions according to User Params\n");
[Sample, digitizer, status] = calcDimensions(US, vars.Sample, vars.SamplingCard);
if(~status); return; end
fprintf("----------Done Calculate Dimensions------------\n");

%------------- Sampling Card ----------------
fprintf("Starting The Sampling Card\n");
[digitizer, status] = digitizerInit(digitizer);
if(~status); return; end
fprintf("Done INIT The Sampling Card\n");
end

