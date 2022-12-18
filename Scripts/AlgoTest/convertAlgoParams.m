function [newParams] = convertAlgoParams(oldParams)
    newParams.fSin              = oldParams.fSin;              
    newParams.fSqnc             = oldParams.fTrain;
    newParams.cycPerPulse       = oldParams.cycInPulse;
    newParams.channels          = oldParams.channels;

    newParams.c                 = oldParams.c;  
    newParams.distFromPhantom   = oldParams.distFromPhantom;

    newParams.fs                = oldParams.fExtClk;

    newParams.fgClk             = oldParams.fSclk;
    newParams.sClkDcyc          = oldParams.extClkDcyc;

    newParams.timeToSample      = oldParams.timeToSample;
    newParams.frameTime         = oldParams.quantTime;
    newParams.envDC             = oldParams.envDC;
    newParams.envUS             = oldParams.envUS;
    newParams.useFrame         = oldParams.useQuant;

    newParams.useGPU              = oldParams.useGPU;
    newParams.useHadamard         = oldParams.useHadamard;
    newParams.contHadamard        = false;
    newParams.highResAO           = oldParams.highResAO;
    newParams.analyzeSingleCh     = false;
    newParams.contSpeckleAnalysis = false;

    newParams.cutArtfct           = false;
    newParams.artfctIdxVec        = [];

    newParams.export.meas           = false;
    newParams.export.signal         = false;
    newParams.export.deMul          = false;
    newParams.export.reshaped       = false;
    newParams.export.fft            = false;
    newParams.export.usCompCmplx    = false;
end