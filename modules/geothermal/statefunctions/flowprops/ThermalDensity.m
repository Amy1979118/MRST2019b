classdef ThermalDensity < StateFunction

    properties
        thermal = true;
        NaCl    = false;
    end

    methods
        function gp = ThermalDensity(model, varargin)
            gp@StateFunction(model, varargin{:});

            gp = gp.dependsOn({'PhasePressures'});
            gp = gp.dependsOn({'temperature'},'state');
            if isprop(model,'NaCl')
                gp = gp.dependsOn({'NaCl'},'state');
            end
        end

        function rho = evaluateOnDomain(prop, model, state)

            p_phase = prop.getEvaluatedDependencies(state, 'PhasePressures');
            T       = model.getProp(state,'temperature');
            f       = model.fluid;
            
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            rho = cell(1, nph);

            if model.water
                wix = phInd == 1;
                pw  = p_phase{wix};
                if isprop(model,'NaCl')
                    c = model.getProp(state,'NaCl');
                    rho{wix} = prop.evaluateFunctionOnDomainWithArguments(f.rhoW, pw, T, c);
                else
                    rho{wix} = prop.evaluateFunctionOnDomainWithArguments(f.rhoW, pw, T);
                end
            end
            if model.oil
                oix = phInd == 2;
                po  = p_phase{oix};
                rho{oix} = prop.evaluateFunctionOnDomainWithArguments(f.rhoO, po, T);
            end
            if model.gas
                gix = phInd == 3;
                pg  = p_phase{gix};
                rho{gix} = prop.evaluateFunctionOnDomainWithArguments(f.rhoG, pg, T);
            end
        end
    end
end

