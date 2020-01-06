classdef GeothermalWaterNaClModel < GeothermalModel
    % Model for single phase liquid (brine) with temperature
    % Equations:
    % {eq1} conservation of mass (water)
    % MAcc + Div(DFlux) - q = 0;
    %
    % DFlux = -k* (kr*rho)/mu (grad(P) - rho*g) : darcy flow, k: permea,
    % kr: rel. perm and mu: viscosity.
    %
    % {eq2} conservation of mass (concentration NaCl in kg/kg)
    % MAcc + Div(NaClFlux) - q = 0; 
    % 
    % NaClFlux = c.*DFlux - phi.*Tau.*dNaCl.*rho_b.* grad(c);
    % molecular diffusion (dispersion is neglected)
    % 
    % {eq3} conservation of energy
    % HAcc + div(HFlux) - q = 0;
    %
    % HFlux = -lambda_r * grad(T) + hf*DFlux 
    % heat flux with hf = fluid enthalpy
    %
    % SYNOPSIS:
    %   model = GeothermalWaterNaClModel(G, rock, fluid)
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
    %   `GeothermalModel`, `GeothermalWaterModel`
      
    properties
        NaCl            %salt present
    end
    
    methods
        function model =  GeothermalWaterNaClModel(G, rock, fluid, varargin)
            model      =  model@GeothermalModel(G, rock, fluid);
            model.NaCl =  true;
            poro       =  rock.poro;
            LS         =  poro.*rock.tau.*fluid.dNaCl;
            fluid_salt =  struct('perm',LS);
            T_salt     =  computeTrans(G,fluid_salt);
            
            % harmonic average of transmissibilities
            cf                            =  G.cells.faces(:,1);
            nf                            =  G.faces.num;
            T_salt                        =  1 ./ accumarray(cf, 1./T_salt, [nf, 1]);
            model.operators.T_salt_all    = T_salt;
            intInx                        = all(G.faces.neighbors ~= 0, 2);
            model.operators.T_salt        = T_salt(intInx);
            model                         =  merge_options(model, varargin{:});
        end
        
        % Define the equation for the problem
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsGeothermalWaterNaCl( state0, state, model, dt, drivingForces, varargin{:});
        end
        
       function [fn, index] = getVariableField(model, name, varargin)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch((name))
            case 'NaCl'
                d     = model.getComponentNames(); 
                index = find(strcmpi(d, 'NaCl')); 
                fn    = 'c';                 
            otherwise
                [fn, index] = getVariableField@GeothermalModel(model, name, varargin{:});
            end 
       end
       
       function names = getComponentNames(model)
            names = getComponentNames@GeothermalModel(model);
            if model.NaCl
                names{end+1} = 'NaCl';
            end
       end
       
       function [state, report] = updateState(model, state, problem, dx, drivingForces)
            state = model.updateStateFromIncrement(state, dx, problem, 'NaCl', Inf, Inf);
            [problem.primaryVariables, removed] = model.stripVars(problem.primaryVariables, 'NaCl');
            dx = dx(~removed);

            [state, report] = updateState@GeothermalModel(model, state, problem, dx, drivingForces);
       end
       
       function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
        %% here the salt is treated as a component
        %  For a given component conservation equation, compute and add in
        %  source terms for a specific source/bc where the fluxes have
        %  already been computed.
        %
        % PARAMETERS:
        %   model  - (Base class, automatic)
        %
        %   cname  - Name of the component. Must be a property known to the
        %            model itself through `getProp` and `getVariableField`.
        %
        %   eq     - Equation where the source terms are to be added. Should
        %            be one value per cell in the simulation grid (model.G)
        %            so that the src.sourceCells is meaningful.
        %
        %   component - Cell-wise values of the component in question. Used
        %               for outflow source terms only.
        %
        %   src    - Source struct containing fields for fluxes etc. Should
        %            be constructed from force and the current reservoir
        %            state by `computeSourcesAndBoundaryConditionsAD`.
        %
        %   force  - Force struct used to produce src. Should contain the
        %            field defining the component in question, so that the
        %            inflow of the component through the boundary condition
        %            or source terms can accurately by estimated.
            
            if isempty(force)
                return
            end
            c = model.getProp(force, cname);
            cells = src.sourceCells;
            switch (cname)
              case 'NaCl'
                % Water based EOR, multiply by water flux divided by
                % density and add into corresponding equation
                qW = src.phaseMass{1}./model.fluid.rhoWS;
                isInj = qW > 0;
                qC = (isInj.*c + ~isInj.*component(cells)).*qW;
              otherwise
                error(['Unknown component ''', cname, '''. BC not implemented.']);
            end
            eq(cells) = eq(cells) - qC;
            src.components{end+1} = qC;
        end

       function [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions@GeothermalModel(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);
            if model.NaCl
                assert(model.water, 'NaCl injection requires a water phase.');
                f = model.fluid;
                if well.isInjector
                    concWell = model.getProp(well.W, 'NaCl');
                else
                    pix = strcmpi(model.getComponentNames(), 'NaCl');
                    concWell = packed.components{pix};
                end
             
                % Water is always first
                wix = 1;
                cqWs = qMass{wix}./f.rhoWS; % connection volume flux at surface condition
                cqP = concWell.*cqWs;

                compSrc{end+1} = cqP;
            end
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
