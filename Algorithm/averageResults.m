function [phiSets, phiFrame, phiSetsStd, phiFrameStd] = averageResults( phiQuant)

phiSets = mean(phiQuant, 4);
phiSetsStd = std(phiQuant, 0, 4);
phiFrame = mean(phiQuant, 3);
phiFrameStd = std(phiQuant, 0, 3);

end