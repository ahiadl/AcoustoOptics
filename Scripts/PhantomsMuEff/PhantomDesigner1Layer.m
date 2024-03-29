%% Phantom Designer
clear all
close all
clc
fprintf("\n");
% References:
% [1] Xu Chen et al. "Development and characterization of tissue-mimicking phantom
% materials for photoacoustic imaging", 2019
% [2] Kent F. Palmer et al. "Optical properties of water in the near
% infrared" 1974
% [3] René Michels ET AL. "Optical properties of fat emulsions", 2008

muaWater      = 0.023; % [cm^-1] from [2]

wl            = 785; % [nm]
Vph           = 200;
musPPh        = 4;
numOfPhantoms = 5;
muaMin        = muaWater;


%% Absorption Coefficients of Ink-Base 1:100 - Royal Talens India Ink
% Interpolate to find the dependency of the phantom extinction coefficient
% vs. its intrinsic concentration of the 1:100 ink base.
% (based on figure 4.5 in [1] )
wl1 = 690; a1 = 0.76;
wl2 = 800; a2 = 0.63;

m = (a2-a1)/(wl2-wl1);
n = a1 - m*wl1;
a = m*wl+n; % [1/cm] Intrinsic absotbtion coefficient of ink per %

muaInkBase100 = a*100; % absorption coefficient of the Ink Base 1:100 solution

fprintf("The absorption coeff. of 1:100 ink-base (Royal Talens) in 785nm is mu_a = %.3f [1/cm] [1]\n", muaInkBase100);

fprintf("\n"); 

%% Scattering of Lipid 20%:
% Anisotropy (g): (Apdx. A, Tab 3, [3])
y0 = 1.09;
a  = -6.812e-4;
g  = y0 + a * wl;

fprintf("The anisotropy coeff. of 20%% lipid in 785nm is g = %.3f [3]\n", g);

% Scattering Coefficient (\mu_s): (Apdx. A, Tab 4, [3])
a = 3.873e8;
b = -2.397e0;

musLipid20 = a * (wl^b) *10;
fprintf("The reduced scattering coeff. of 20%% lipid in 785nm is mu_s = %.3f [1/cm] [3]\n", musLipid20);

% Reduced Scattering Coefficient (\mu_s'): (Apdx. A, Tab 5, [3])
y0 = 8.261e1;
a  = -1.288e-1;
b  = 6.093e-5;
wl = 780;

musPLipid20 = y0+a*wl+b*(wl^2); %[1/mm], (lambda in nm)
musPLipid20 = musPLipid20*10; %[1/cm] 

fprintf("The reduced scattering coeff. of 20%% lipid in 785nm is mu_s^p = %.3f [1/cm] [3]\n", musPLipid20);
fprintf("Reminder: Mus\' = Mus(1-g). Calculated Mus\'=%.2f\n", musLipid20*(1-g));
fprintf("\n");

MusPh = musPPh / (1-g);

%% Scattering Base
Clipid = musPPh ./ musPLipid20;
Vlipid20 = Clipid .* Vph;
VWaterPh   = Vph - Vlipid20;

Vlipid20Total = Vlipid20 * (numOfPhantoms+1);
VWaterTotal   = VWaterPh   * (numOfPhantoms+1);

fprintf("To create the scattering base mix: \n%.2f [m] of 20%% lipid with \n%.2f[ml] pure water.\n", Vlipid20Total,VWaterTotal);  
fprintf("\n");

%% Create a 20% Lipid solution using a 38% milk cream
% Explanation:
% fatF = VfMC * CfMC;    % weight of lipid in final solution
% ViMC = fatF/CiMC;      % volume of MC to use to get desired weight of lipids
% ViMC = VfMC*CfMC/CiMC;

VfMC = Vlipid20Total*2;   % Final [ml]
CiMC = 0.38; % Initial fat %
CfMC = 0.2;  % Final fat %
ViMC = (CfMC/CiMC) * VfMC; % Milk cream volume that need to mix in base % [ml]
ViW  = VfMC - ViMC;   % volume of water to mix with MC; 

fprintf("To create a %.2f [ml] of 20%% lipid solution using a 38%% milk cream: \nMix a %.2f [ml] of milk cream with %.2f[ml] of water\n", VfMC, ViMC, ViW);
fprintf("\n");

%% Phantoms
muaPh   = muaMin* 2.^[1:(numOfPhantoms)];
muEffPh = sqrt(3*muaPh.*(muaPh + musPPh));
% muaPh = muEffPh.^2 ./(3.*musPPh); %[1/cm]

% Explanation:
% muaPh = Cinkbase * muaInkBase100 + (1-Cinkbase)*muaWater;
% muaPh = Cinkbase * muaInkBase100 + muaWater-Cinkbase*muaWater;
% Cinkbase * muaInkBase100 - Cinkbase*MuaWater = muaPh - muaWater
% Cinkbase *(muaInkBase100 - muaWater) = muaPh - muaWater
% Cinkbase  = (muaPh - muaWater) / (muaInkBase100 - muaWater);
Cinkbase =  (muaPh - muaWater) ./(muaInkBase100 - muaWater);

% Explanation:
% Cinkbase = Vinkbase/(Vinkbase+Vbase)
% Vinkbase = Cinkbase*Vinkbase + Cinkbase*Vbase
% Vinkbase(1-Cinkbase) =  Cinkbase*Vbase;
Vinkbase = Cinkbase*Vph ./ (1-Cinkbase);
VtotalPh = Vinkbase + Vph;

%% Display Results
for i = 1:numOfPhantoms
    phName{i}  = num2str(i);
end

MuA         = muaPh';
MuS         = MusPh*ones(numOfPhantoms, 1);
MuSP        = musPPh*ones(numOfPhantoms, 1);
MuEff       = muEffPh';
MFP         = 1./(MuA+MuSP);
VWater      = VWaterPh*ones(numOfPhantoms,1);
Vlipid      = Vlipid20*ones(numOfPhantoms,1);
VInkBase    = Vinkbase';
TotalVolume = VtotalPh';

T = table(phName',MuA,MuS,MuSP,MuEff,MFP);
disp(T)
fprintf("All Table Fields are in [1/cm]\n");
fprintf("MFP Field are in [cm]\n");
fprintf("\n");

T = table(phName',VWater,Vlipid,VInkBase,TotalVolume);
disp(T)
fprintf("All Table Fields are in [ml]\n");
fprintf("\n");