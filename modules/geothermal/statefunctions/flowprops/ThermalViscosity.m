classdef ThermalViscosity < StateFunction
        
    properties
        thermal = true;
        NaCl    = false;
    end
    
    methods
        function gp = ThermalViscosity(model, varargin)
            gp@StateFunction(model, varargin{:});
            
            gp = gp.dependsOn({'PhasePressures'});
            gp = gp.dependsOn({'temperature'},'state');
            if isprop(model,'NaCl')
                gp = gp.dependsOn({'NaCl'},'state');
            end
        end
        
        function mu = evaluateOnDomain(prop,model, state)
            p_phase = prop.getEvaluatedDependencies(state, 'PhasePressures');
            T       = model.getProp(state,'temperature');
            f       = model.fluid;
            
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            mu  = cell(1, nph);
            
            if model.water
                wix = phInd == 1;
                pw  = p_phase{wix};
                if isprop(model,'NaCl')
                    c = model.getProp(state,'NaCl');
                    mu{wix} = prop.evaluateFunctionOnDomainWithArguments(f.muW, pw, T, c);
                else
                    mu{wix} = prop.evaluateFunctionOnDomainWithArguments(f.muW, pw, T);
                end
            end
            if model.oil
                oix = phInd == 2;
                po  = p_phase{oix};
                mu{oix} = prop.evaluateFunctionOnDomainWithArguments(f.muO, po, T);
            end
            if model.gas
                gix = phInd == 3;
                pg  = p_phase{gix};
                mu{gix} = prop.evaluateFunctionOnDomainWithArguments(f.muG, pg, T);
            end
        end
    end   
end
