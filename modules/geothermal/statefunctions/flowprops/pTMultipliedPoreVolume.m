classdef pTMultipliedPoreVolume < StateFunction
    
    properties
    thermal = true;
    end
    
    methods
        function gp = pTMultipliedPoreVolume(model, varargin)
            gp@StateFunction(model, varargin{:});
            if isfield(model.fluid, 'pvMultR')
                gp = gp.dependsOn({'pressure, temperature'}, 'state');
            end
        end
        function pv = evaluateOnDomain(prop, model, state)
            % Get effective pore-volume, accounting for possible
            % rock-compressibility-dilatation as a function of p or T
            f = model.fluid;
            pv = model.operators.pv;
            if isfield(f, 'pvMultR')
                p = model.getProp(state, 'pressure');
                T = model.getProp(state, 'temperature');
                pvMult = prop.evaluateFunctionOnDomainWithArguments(f.pvMultR, p, T);
                pv = pv.*pvMult;
            end
        end
    end
end