classdef GeothermalWaterModel < GeothermalModel
    % DESCRIPTION:
    % Model for single liquid phase (pure water) with temperature
    % Equations:
    % {eq1} conservation of mass
    % MAcc + Div(DFlux) - q = 0;
    %
    % DFlux = -k* (kr*rho)/mu (grad(P) - rho*g) : darcy flow, k: permea,
    % kr: rel. perm and mu: viscosity.
    %
    % {eq2} conservation of energy
    % HAcc + div(HFlux) - q = 0;
    %
    % HFlux = -lambda_r * grad(T) + hf*DFlux 
    % heat flux with hf = fluid enthalpy
    % 
    % SYNOPSIS:
    %   model = GeothermalWaterModel(G, rock, fluid)
    %
    % REQUIRED PARAMETERS:
    %   G     - Simulation grid.
    %
    %   rock  - Valid rock used for the model.
    %
    %   fluid - Fluid model used for the model.
    %
    %
    % OPTIONAL PARAMETERS:
    %   'property' - Set property to the specified value.
    %
    % RETURNS:
    %   Class instance.
    %
    % SEE ALSO:
    %   `GeothermalModel`, `GeothermalWaterNaClModel`
    

    properties
    end
    
    methods
         function model = GeothermalWaterModel(G, rock, fluid, varargin)
            model       =   model@GeothermalModel(G, rock, fluid);
            model       =   merge_options(model, varargin{:});
         end 
         
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsGeothermalWater( state0, state, model, dt, drivingForces, varargin{:});
        end
    end
    
end

%{
Copyright 2009-2020 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
