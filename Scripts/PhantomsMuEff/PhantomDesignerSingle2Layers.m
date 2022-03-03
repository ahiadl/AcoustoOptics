%% Phantom Designer
clear all
close all
clc
% fprintf("\n");
% References:
% [1] Xu Chen et al. "Development and characterization of tissue-mimicking phantom
% materials for photoacoustic imaging", 2019
% [2] Kent F. Palmer et al. "Optical properties of water in the near
% infrared" 1974
% [3] René Michels ET AL. "Optical properties of fat emulsions", 2008

%Mold Parameters:
Amold   = 60*60; %[mm^2]
Dmold   = 50;    %[mm]
Dlayer1 = 5;    %[mm];

muaWater      = 0.023;    %[cm^-1] from [2]
wl            = 785;      %[nm]
musPPh        = 4;        %[cm^-1]
muaMin = muaWater;

mua1 = 0.046;
mua2 = 0.736;

%% Calculate Each layer volume in ml
VPh     = Amold * Dmold   * 1e-3; 
Vlayer1 = Amold * Dlayer1 * 1e-3;
Vlayer2 = VPh - Vlayer1;

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
fprintf("The absorption coeff. of water in 785nm is mu_a = %.3f [1/cm] [1]\n", muaWater);

%% Scattering Coefficient of Lipid 20%
% The following parameters are from appendix A, tab 5 in [3]
y0 = 8.261e1;
a  = -1.288e-1;
b  = 6.093e-5;
wl = 780;

% The following formula is according to eq. 11 in [3]
musPLipid20 = y0+a*wl+b*wl^2; %[1/mm], (lambda in nm)
musPLipid20 = musPLipid20*10; %[1/cm] 

fprintf("The reduced scattering coeff. of 20%% lipid in 785nm is mu_s^p = %.3f [1/cm] [3]\n\n", musPLipid20);


%% Scattering Base
Clipid = musPPh ./ musPLipid20;
Vlipid20 = Clipid .* VPh;
VWaterPh   = VPh - Vlipid20;

Vlipid20Total = Vlipid20 * 2;
VWaterTotal   = VWaterPh   * 2;

%% Create a 20% Lipid solution using a 38% milk cream 
VfMC = Vlipid20Total*2;   % Final [ml]
CiMC = 0.38; % Initial fat %
CfMC = 0.2;  % Final fat %
ViMC = (CfMC/CiMC) * VfMC; % Milk cream volume that need to mix in base % [ml]
ViW  = VfMC - ViMC;   % volume of water to mix with MC; 

% Explanation:
% fatF = VfMC * CfMC;    % weight of lipid in final solution
% ViMC = fatF/CiMC;      % volume of MC to use to get desired weight of lipids
% ViMC = VfMC*CfMC/CiMC; 

agarWeight = (VWaterTotal+Vlipid20Total)/100*1.5;

fprintf("To create a %.2f [ml] of 20%% lipid solution using a 38%% milk cream mix: \n", VfMC);
fprintf("%.2f [ml] of milk cream \n%.2f [ml] of water\n\n", ViMC, ViW);

fprintf("To create the scattering base mix: \n");
fprintf("%.2f [m] of 20%% lipid \n%.2f[ml] pure water \n%.2f [g] of agar.\n\n", Vlipid20Total,VWaterTotal, agarWeight);  

%% Phantoms
muEffLayer1 = sqrt(3*mua1.*(mua1 + musPPh));
muEffLayer2 = sqrt(3*mua2.*(mua2 + musPPh));

% muaPh = muEffPh.^2 ./(3.*musPPh); %[1/cm]

% Explanation:
% muaPh = Cinkbase * muaInkBase100 + Cwater * muaWater;
% muaPh = Cinkbase * muaInkBase100 + (1-Cinkbase)*muaWater;
% muaPh = Cinkbase * muaInkBase100 + muaWater-Cinkbase*muaWater;
% Cinkbase * muaInkBase100 - Cinkbase*MuaWater = muaPh - muaWater
% Cinkbase *(muaInkBase100 - muaWater) = muaPh - muaWater
% Cinkbase  = (muaPh - muaWater) / (muaInkBase100 - muaWater);
Cinkbase1 =  (mua1 - muaWater) ./(muaInkBase100 - muaWater);
Cinkbase2 =  (mua2 - muaWater) ./(muaInkBase100 - muaWater);


% Explanation:
% Cinkbase = Vinkbase/Vph
% Cinkbase = Vinkbase/(Vinkbase+Vscaterbase)
% Vinkbase = Cinkbase*Vinkbase + Cinkbase * Vscaterbase;
% Vinkbase(1-Cinkbase) =  Cinkbase * Vscaterbase;
% Vinkbase = ( Cinkbase / (1-Cinkbase) ) * Vscaterbase ;
VinkLayer1   = (Cinkbase1 ./ (1-Cinkbase1)) * Vlayer1 ;
VtotalLayer1 = VinkLayer1 + Vlayer1;

VinkLayer2   = (Cinkbase2 ./ (1-Cinkbase2)) * Vlayer2 ;
VtotalLayer2 = VinkLayer2 + Vlayer2;

fprintf("To create layer 1:\n")
fprintf("Mix %.2f ml from the total base \nwith %.4f [ul] of ink base.\n\n", Vlayer1, VinkLayer1*1e3);


fprintf("To create layer 2:\n")
fprintf("Mix %.2f ml from the total base \nwith %.4f [ul] of ink base.\n\n", Vlayer2, VinkLayer2*1e3);

%% Display Results

% MuA         = muaLayer1;
% MuSP        = musPPh;
% MuEff       = muEffLayer1;
% fprintf("Layer 1 Properties: \n");
% T = table("Layer 1", MuA, MuSP, MuEff, ...
%            'VariableNames',{'Index','MuA [1/cm]', 'MuS Prime [1/cm]', 'Mu Eff [[1/cm]'});
% disp(T)
% 
% MuA         = muaLayer2';
% MuSP        = musPPh*ones(numOfPhantoms, 1);
% MuEff       = muEffLayer2';
% 
% fprintf("Layer 2 Properties: \n");
% T = table((1:numOfPhantoms)', MuA, MuSP, MuEff,...
%           'VariableNames',{'Index','MuA [1/cm]', 'MuS Prime [1/cm]', 'Mu Eff [[1/cm]'});
% disp(T)
% 
% fprintf("Layer 2 Volumes: \n");
% Vscatbase   = Vlayer2*ones(numOfPhantoms,1);
% VInkBase    = VinkLayer2'*1e3;
% TotalVolume = VtotalLayer2';
% T = table((1:numOfPhantoms)',Vscatbase,VInkBase,TotalVolume,...
%            'VariableNames',{'Index', 'Scat. base [ml]', 'Ink Base [ul]', 'Total[ml]'});
% disp(T)