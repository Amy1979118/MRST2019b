classdef HeatTransFluid < StateFunction
    properties
        thermal = true;
    end
    
    methods
        function pp = HeatTransFluid(model)
            pp@StateFunction(model);
            if isfield(model.fluid, 'transMultH') || isfield(model.fluid, 'pvMult')
                pp = pp.dependsOn({'pressure, temperature'}, 'state');
            end
        end
        
        function THF = evaluateOnDomain(prop, model, state)
            THF = model.operators.T_heat_F;
            p   = model.getProp(state, 'pressure');
            t   = model.getProp(state, 'temperature');
            if isfield(model.fluid, 'transMultH')  
                THF = model.fluid.transMultH(p,t).*THF;
            end
            
            if isfield(model.fluid, 'pvMult')  
                THF = model.fluid.pvMult(p,t).*THF;
            end
        end
    end
end

