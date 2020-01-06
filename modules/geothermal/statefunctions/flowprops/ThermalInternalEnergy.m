classdef ThermalInternalEnergy < StateFunction
       
    properties
        thermal = true;
        NaCl    = false;
    end
    
    methods
        function gp = ThermalInternalEnergy(model, varargin)
            gp@StateFunction(model, varargin{:});
            
            gp = gp.dependsOn({'PhasePressures'});
            gp = gp.dependsOn({'temperature'},'state');
            if isprop(model,'NaCl')
                gp = gp.dependsOn({'NaCl'},'state');
            end
        end
        
        function u = evaluateOnDomain(prop,model, state)
            p_phase = prop.getEvaluatedDependencies(state, 'PhasePressures');
            T       = model.getProp(state,'temperature');
            f       = model.fluid;
            
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            u   = cell(1, nph);
            
            if model.water
                wix = phInd == 1;
                pw  = p_phase{wix};
                if isprop(model,'NaCl')
                    c = model.getProp(state,'NaCl');
                    u{wix} = prop.evaluateFunctionOnDomainWithArguments(f.uW, pw, T, c);
                else
                    u{wix} = prop.evaluateFunctionOnDomainWithArguments(f.uW, pw, T);
                end
            end
            if model.oil
                oix = phInd == 2;
                po  = p_phase{oix};
                u{oix} = prop.evaluateFunctionOnDomainWithArguments(f.uO, po, T);
            end
            if model.gas
                gix = phInd == 3;
                pg  = p_phase{gix};
                u{gix} = prop.evaluateFunctionOnDomainWithArguments(f.uG, pg, T);
            end
        end       
    end
end

