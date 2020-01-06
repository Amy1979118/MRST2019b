classdef ThermalFluxDiscretization < FluxDiscretization
    
    properties
        HeatTransRock
        HeatTransFluid
        HeatFluxRock
        HeatFluxFluid
    end
    methods
        function props = ThermalFluxDiscretization(model)
            props = props@FluxDiscretization(model);
            
            props.HeatTransRock    = HeatTransRock(model);
            props.HeatTransFluid   = HeatTransFluid(model);
            props.HeatFluxRock     = HeatFluxRock(model);
            props.HeatFluxFluid    = HeatFluxFluid(model);
            props.Transmissibility = pTTransmissibility(model);

        end
    end  
end