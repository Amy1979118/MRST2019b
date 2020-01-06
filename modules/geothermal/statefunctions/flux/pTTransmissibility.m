classdef pTTransmissibility < StateFunction
    properties
        thermal = true;
    end
    
    methods
        function pp = pTTransmissibility(model)
            pp@StateFunction(model);
            if isfield(model.fluid, 'transMult')
                pp = pp.dependsOn({'pressure, temperature'}, 'state');
            end
        end
        
        function T = evaluateOnDomain(prop, model, state)
            T = model.operators.T;
            if isfield(model.fluid, 'transMult')
                p = model.getProp(state, 'pressure');
                t = model.getProp(state, 'temperature');
                T = model.fluid.transMult(p,t).*T;
            end
        end
    end
end