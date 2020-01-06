classdef HeatFluxFluid < StateFunction
   
    properties
        thermal = true; 
    end
    
    methods
        function gp = HeatFluxFluid(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'temperature'},'state');
            gp = gp.dependsOn({'HeatTransFluid'});
        end
        
        function HFlux_f = evaluateOnDomain(prop,model, state)
            op = model.operators;
            T  = model.getProp(state,'temperature');
            T_heat_F = prop.getEvaluatedDependencies(state, 'HeatTransFluid');
            HFlux_f  = -T_heat_F.*op.Grad(T);
        end 
    end
    
end

