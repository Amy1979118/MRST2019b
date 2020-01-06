
classdef ThermalEnthalpy < StateFunction
        
    properties
        thermal = true;
        NaCl    = false;
    end
    
    methods
        function gp = ThermalEnthalpy(model, varargin)
            gp@StateFunction(model, varargin{:});
            
            gp = gp.dependsOn({'PhasePressures'});
            gp = gp.dependsOn({'temperature'},'state');
            if isprop(model,'NaCl')
                gp = gp.dependsOn({'NaCl'},'state');
            end
        end
        
        function h = evaluateOnDomain(prop,model, state)
            p_phase = prop.getEvaluatedDependencies(state, 'PhasePressures');
            T       = model.getProp(state,'temperature');
            f       = model.fluid;
            
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            h   = cell(1, nph);
            
            if model.water
                wix = phInd == 1;
                pw  = p_phase{wix};
                if isprop(model,'NaCl')
                    c = model.getProp(state,'NaCl');
                    h{wix} = prop.evaluateFunctionOnDomainWithArguments(f.hW, pw, T, c);
                else
                    h{wix} = prop.evaluateFunctionOnDomainWithArguments(f.hW, pw, T);
                end
            end
            if model.oil
                oix = phInd == 2;
                po  = p_phase{oix};
                h{oix} = prop.evaluateFunctionOnDomainWithArguments(f.hO, po, T);
            end
            if model.gas
                gix = phInd == 3;
                pg  = p_phase{gix};
                h{gix} = prop.evaluateFunctionOnDomainWithArguments(f.hG, pg, T);
            end
        end
    end
end

