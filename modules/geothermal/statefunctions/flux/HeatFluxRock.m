classdef HeatFluxRock < StateFunction
  
    properties
        thermal = true; 
    end
    
    methods
        function gp = HeatFluxRock(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'temperature'},'state');
            gp = gp.dependsOn({'HeatTransRock'});
        end
        
        function HFlux_r = evaluateOnDomain(prop,model, state)
            op = model.operators;
            T  = model.getProp(state,'temperature');
            T_heat_R = prop.getEvaluatedDependencies(state, 'HeatTransRock');
            HFlux_r  = -T_heat_R.*op.Grad(T);
        end 
    end   
end

