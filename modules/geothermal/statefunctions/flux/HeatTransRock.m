classdef HeatTransRock < StateFunction
    properties
        thermal = true;
    end
    
    methods
        function pp = HeatTransRock(model)
            pp@StateFunction(model);
            if isfield(model.rock, 'transMultH') || isfield(model.fluid, 'pvMult')
                pp = pp.dependsOn({'pressure, temperature'}, 'state');
            end
        end
        
        function THR = evaluateOnDomain(prop, model, state)
            THR = model.operators.T_heat_R;
            p   = model.getProp(state, 'pressure');
            t   = model.getProp(state, 'temperature');
            if isfield(model.rock, 'transMultH')  
                THR = model.rock.transMultH(p,t).*THR;
            end
            
            if isfield(model.fluid, 'pvMult')  
                THR = model.fluid.pvMult(p,t).*THR;
            end
        end
    end
end

