classdef ThermalInternalEnergyRock < StateFunction
   
    properties
        thermal = true; 
    end
    
    methods
        function gp = ThermalInternalEnergyRock(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'pressure, temperature'},'state');
        end
        
        function uR = evaluateOnDomain(prop,model, state)
            r     = model.rock;
            p     = model.getProp(state,'pressure');
            T     = model.getProp(state,'temperature');
            CpR   = prop.evaluateFunctionOnDomainWithArguments(r.CpR, p, T);
            rhoR  = prop.evaluateFunctionOnDomainWithArguments(r.rhoR, p, T);
            uR    = T.*CpR.*rhoR; 
        end 
    end
    
end

