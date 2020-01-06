classdef ThermalFlowPropertyFunctions < FlowPropertyFunctions
    properties
        ThermalInternalEnergy
        ThermalEnthalpy
        ThermalInternalEnergyRock
        
    end
    methods
        function props = ThermalFlowPropertyFunctions(model)
            props = props@FlowPropertyFunctions(model);
            
            props.ThermalInternalEnergy     = ThermalInternalEnergy(model);
            props.ThermalEnthalpy           = ThermalEnthalpy(model);
            props.ThermalInternalEnergyRock = ThermalInternalEnergyRock(model);
            props.Density                   = ThermalDensity(model);
            props.Viscosity                 = ThermalViscosity(model);
            props.PoreVolume                = pTMultipliedPoreVolume(model);
        end        
    end    
end