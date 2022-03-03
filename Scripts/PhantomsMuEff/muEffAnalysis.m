classdef muEffAnalysis
    %MUEFFANALYSIS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vars
    end
    
    methods
        function this = muEffAnalysis()

            
        end
        
        function loadSingleLayerSet(this, simulationFile)
            sim = load(simulationFile);
            this.vars.layers   = layers;
            this.vars.depthVec = sim.depthVecHR;
            numOfPhantoms = length(sim.resArr);
            
            if layers == 1
                for i=1:numOfPhantoms
                   phi(:,i) = 
                end
                
                
            elseif layers ==2
                
            else
                fprintf("Does not support number of layers larger than 2/n.");
                return;
            end
            
            for i = 1:numOfPhantoms
                phi(i,:) = sim.resArr{i}.phi;
            end
             
        end
       
       function analayseSingleLayerSet()

       function loadMeasurement(this)

 
       end
       
       
    end
end

