classdef PhantomDesigner < handle
    % References:
    % [1] Xu Chen et al. "Development and characterization of tissue-mimicking phantom
    % materials for photoacoustic imaging", 2019
    % [2] Kent F. Palmer et al. "Optical properties of water in the near
    % infrared" 1974
    % [3] René Michels ET AL. "Optical properties of fat emulsions", 2008
    % [4] L.V. wang, "Biomedical Optics Principles", 2009 (Book)
    properties
        vars
    end
    
    methods
        function this = PhantomDesigner()
            this.vars.muaWater = 0.0023; % [mm^-1] from [2]
            this.vars.muaMin   = this.vars.muaWater;
            this.vars.moldSize = 60; % mm;
            fprintf("The absorption coeff. of water in 785nm is mu_a = %.3f [1/cm] [2]\n", this.vars.muaWater);
        end
        
        function setVars(this, uVars)
           this.vars.wl   = uVars.wl;
           this.vars.Vph  = uVars.Vph;
           this.vars.musP = uVars.musP;
           this.vars.mua  = uVars.mua;
           
           this.vars.muEff = sqrt(3*this.vars.mua.*(this.vars.mua + this.vars.musP)); % [4], p.99 (108)
           this.vars.MFP   = 1./(this.vars.mua + this.vars.musP);
           
           this.vars.inkType = uVars.inkType;
           
           this.vars.layerThick = uVars.layerThick; %[mm]
           this.vars.layerVol = this.vars.layerThick*(this.vars.moldSize^2) * 1e-3; % [ml]
           
           this.vars.numLayers  = size(this.vars.mua,1);
           this.vars.numPhantom = size(this.vars.mua, 2);
           
           this.calcInkProperties(this.vars.wl);
           this.calcLipid20Properties(this.vars.wl);
           
           this.vars.mus = this.vars.musP / (1-this.vars.g);
        end
        
        function calcInkProperties(this, wl)
            switch this.vars.inkType
                case 'RoyalTalens'
                    % Interpolate to find the dependency of the phantom extinction coefficient
                    % vs. its intrinsic concentration of the 1:100 ink base.
                    % (based on figures 4.4 in [1] )
                    
                    % Linear Interpolation:
                    wl1 = 690; 
                    a1  = 0.076;
                    wl2 = 800; 
                    a2  = 0.063;
                    
                    % For the given wavelength, and 1:100 base: What is the
                    % extinction coefficient when the *base* is further 
                    % dilluted with water: mu_a = a*Cbase;
                    m = (a2-a1)/(wl2-wl1);
                    n = a1 - m*wl1;
                    a = m*wl+n; 
                    
                    % Calculate the absorption coefficient of pure 1:100
                    % base:
                    this.vars.muaInkBase = a*100; %[mm^-1]
                    
                    fprintf("The absorption coeff. of 1:100 ink-base (Royal Talens) in %d nm is mu_a = %.3f [1/cm] [1]\n\n", wl, this.vars.muaInkBase);
            end
        end
        
        function calcLipid20Properties(this, wl)
            % Anisotropy (g): (Apdx. A, Tab 3, [3])
            y0 = 1.09;
            a  = -6.812e-4;
            g  = y0 + a * wl;
            
            this.vars.lipid.model.anisotropy.y0 = y0;
            this.vars.lipid.model.anisotropy.a  = a;
            this.vars.g = g;
            
            fprintf("The anisotropy coeff. of 20%% lipid in 785nm is g = %.3f [3]\n", g);
            
            % Scattering Coefficient (\mu_s): (Apdx. A, Tab 4, [3])
            a = 3.873e8;
            b = -2.397e0;
            
            musLipid20 = a * (wl^b);

            this.vars.lipid.model.mus.a  = a;
            this.vars.lipid.model.mus.b = b;
            this.vars.musLipid20 = musLipid20;
            
            fprintf("The reduced scattering coeff. of 20%% lipid in 785nm is mu_s = %.3f [1/mm] [3]\n", musLipid20);
            
            % Reduced Scattering Coefficient (\mu_s'): (Apdx. A, Tab 5, [3])
            y0 = 8.261e1;
            a  = -1.288e-1;
            b  = 6.093e-5;

            musPLipid20 = y0+a*wl+b*(wl^2); %[1/mm], (lambda in nm)
            
            this.vars.lipid.model.musP.a  = a;
            this.vars.lipid.model.musP.b  = b;
            this.vars.lipid.model.musP.y0 = y0;
            
            this.vars.musPLipid20 = musPLipid20;
            
            fprintf("The reduced scattering coeff. of 20%% lipid in 785nm is mu_s^p = %.3f [1/mm] [3]\n", musPLipid20);
            fprintf("Reminder: Mus\' = Mus(1-g). Calculated Mus\'=%.2f\n\n", musLipid20*(1-g));
        end

        function res = calcLipid20VolLayer(this, musP, vLayer)
            cLipid   = musP ./ this.vars.musPLipid20;
            vLipid20 = cLipid .* vLayer;
            vWater   = vLayer - vLipid20;

            res.vLipid = vLipid20;
            res.vWater = vWater;
        end
        
        function res = calcInkVolLayer(this, mua, vLayer)
            muaWater = this.vars.muaWater;
            muaInkBase = this.vars.muaInkBase;
            
            % Assumption: muA_Lipid20|_lambda>500 = mua_Water
            
            % Explanation:
            % mua = cInkbase * muaInkBase + (1-Cinkbase)*muaWater;
            % mua = cInkbase * muaInkBase + muaWater-Cinkbase*muaWater;
            % cInkbase * muaInkBase - cInkbase*MuaWater = mua - muaWater
            % cInkbase *(muaInkBase - muaWater) = mua - muaWater
            % cInkbase  = (mua - muaWater) / (muaInkBase - muaWater);
            cInkBase =  (mua - muaWater) ./(muaInkBase - muaWater);

            % Explanation:
            % cInkbase = vInkbase/(vInkbase + vLayer)
            % Vinkbase = cInkbase*vInkbase + cInkbase*vLayer
            % Vinkbase(1-cInkbase) =  cInkbase * vBase;
            vInkBase = cInkBase * vLayer ./ (1-cInkBase);
            vTotal = vInkBase + vLayer;
            
            res.cInkBase = cInkBase;
            res.vInkBase = vInkBase;
            res.vTotal   = vTotal;
            
        end
        
        function calcLipid20Base(this, vLipid20)
            % Explanation:
            % fatF = vLipid20 * CLipF;    % weight of lipid in final solution
            % vMC = fatF/CLipI;           % volume of MC to use to get desired weight of lipids
            % vMC = vLipid20*CLipF/CLipI; % volume of MC to use calculated with relation to desired volume of lipid20 solution.
            
            vLipid20 = 2*vLipid20;
            
            CLipI = 0.38; % Initial [%] of lipids in MC
            CLipF = 0.2;  % Final [%] of lipid in solution
            vMC = (CLipF/CLipI) * vLipid20; % Milk cream volume that need to mix in base % [ml]
            vWater  = vLipid20 - vMC;   % volume of water to mix with MC; 
            
            fprintf("To create a %.2f [ml] of 20%% lipid solution using a 38%% milk cream: \nMix a %.2f [ml] of milk cream with %.2f[ml] of water\n\n", vLipid20, vMC, vWater);
        end
        
        function createPhantomSet(this, uVars)
            this.setVars(uVars)

            for j = 1:this.vars.numLayers
                curLayerVol = this.vars.layerVol(j);
                for i = 1:this.vars.numPhantom
                    curMusP = this.vars.musP(j,i);
                    curMua = this.vars.mua(j,i);
                    
                    resScat = this.calcLipid20VolLayer(curMusP, curLayerVol);
                    resAbs  = this.calcInkVolLayer(curMua, curLayerVol);
                    
                    vLipid20(j,i) = resScat.vLipid;
                    vWater(j,i)   = resScat.vWater;
                    vInkBase(j,i) = resAbs.vInkBase;
                    
                    totalVol(j,i) = vWater(j,i) + vInkBase(j,i) + vLipid20(j,i);

                end
            end
            vLipidTotal = sum(vLipid20(:));
            vWaterTotal = sum(vWater(:));
            this.calcLipid20Base(vLipidTotal);
            
            fprintf("Total volume of water needed: %.2f [ml] mixed with %.2f [g] of agar \n\n", vWaterTotal*1.2, vWaterTotal*1.2/100 *1.5);
            
            varNamesVols = {'Phantom idx' ; 'V-Water'; 'V-Lipid 20%%'; 'V-InkBase'; 'Total Volume'};
            varNamesOP =  {'Phantom idx' ; 'MuA'; 'MuS'; 'MuSP'; 'MuEff'; 'MFP'};
            phIdx = 1:this.vars.numPhantom;
            for j=1:this.vars.numLayers
                fprintf("Layer %d:\n", j)
                T = table(phIdx', this.vars.mua(j,:)', this.vars.musP(j,:)', this.vars.mus(j,:)', this.vars.muEff(j,:)', this.vars.MFP(j,:)', 'VariableNames', varNamesOP)
                T = table(phIdx', vWater(j,:)', vLipid20(j,:)', vInkBase(j,:)', totalVol(j,:)', 'VariableNames', varNamesVols)
                fprintf("\n");
            end
        end
    end
end

