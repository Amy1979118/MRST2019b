classdef FluxDiscretization < StateFunctionGrouping
    % Function grouping for Darcy-type flux discretization. The defaults
    % gives a industry-standard single-point upwind scheme with a two-point
    % flux discretization which emphasizes robustness and efficiency.
    properties
        PermeabilityPotentialGradient % K * (grad(p) + rho g dz)
        PressureGradient % Gradient of phase pressures
        GravityPotentialDifference % rho * g * dz term
        PhasePotentialDifference % (grad p_alpha + dpdz)
        PhaseFlux % Phase volumetric fluxes
        FaceComponentMobility % Composition * mobility on face
        PhaseUpwindFlag % Upwind flag for each phase
        ComponentTotalFlux % Total mass flux for each component
        ComponentPhaseFlux % Phase fluxes for each component
        Transmissibility % Face-based transmissibility
    end

    properties (Access = protected)
        FlowStateBuilder
    end
    
    methods
        function props = FluxDiscretization(model)
            upstr = UpwindFunctionWrapperDiscretization(model);
            tpfa = TwoPointFluxApproximation(model);

            props@StateFunctionGrouping();
            % Darcy flux
            props = props.setStateFunction('Transmissibility', Transmissibility(model));
            props = props.setStateFunction('PermeabilityPotentialGradient', PermeabilityPotentialGradient(model, tpfa));
            if ~isempty(model.inputdata) && isfield(model.inputdata.SOLUTION, 'THPRES')
                pgrad = PressureGradientWithThresholdPressure(model, model);
            else
                pgrad = PressureGradient(model);
            end
            props = props.setStateFunction('PressureGradient', pgrad);
            props = props.setStateFunction('GravityPotentialDifference', GravityPotentialDifference(model));
            
            % Phase flux
            props = props.setStateFunction('PhasePotentialDifference', PhasePotentialDifference(model));
            props = props.setStateFunction('PhaseUpwindFlag', PhaseUpwindFlag(model));
            
            % Face values - typically upwinded
            props = props.setStateFunction('FaceComponentMobility', FaceComponentMobility(model, upstr));
            % Phase mobility on face
            props = props.setStateFunction('FaceMobility', FaceMobility(model, upstr));
            % 
            props = props.setStateFunction('ComponentPhaseFlux', ComponentPhaseFlux(model));
            props = props.setStateFunction('ComponentTotalFlux', ComponentTotalFlux(model));
            props = props.setStateFunction('PhaseFlux', PhaseFlux(model));

            % Flow discretizer
            props.FlowStateBuilder = ImplicitFlowStateBuilder();
            
            % Define storage
            props.structName = 'FluxProps';
        end

        function fd = setFlowStateBuilder(fd, fb)
            fd.FlowStateBuilder = fb;
        end
        
        function fb = getFlowStateBuilder(fd)
            fb = fd.FlowStateBuilder;
        end
        
        function [acc, v, names, types] = componentConservationEquations(fd, model, state, state0, dt)
            % Compute discretized conservation equations in the interior of the domain.
            % REQUIRED PARAMETERS:
            %   model  - ReservoirModel derived class
            %   state  - State corresponding to the given model at the end
            %            of the current time-integration interval.
            %   state0 - State at previous time-step.
            %   dt     - Time-step length.
            % RETURNS:
            %   acc    - Cell array of component accumulation terms
            %            discretized with a first order finite-difference
            %            scheme.
            %   v      - Cell array of component fluxes in the interior of
            %            the domain.
            %   names  - The names of each equation (corresponds to
            %            component names.
            %   types  - Types of each equation (defaults to 'cell', to
            %            indicate cell-wise equations on the whole domain)
            ncomp = model.getNumberOfComponents;
            [acc, types] = deal(cell(1, ncomp));
            names = model.getComponentNames();
            [types{:}] = deal('cell');
            mass = model.getProps(state, 'ComponentTotalMass');
            mass0 = model.getProps(state0, 'ComponentTotalMass');
            flowState = fd.buildFlowState(model, state, state0, dt);
            v = model.getProps(flowState, 'ComponentTotalFlux');
            for c = 1:ncomp
                acc{c} = (mass{c} - mass0{c})./dt;
            end
        end

        function flowState = buildFlowState(fd, model, state, state0, dt)
            flowState = fd.FlowStateBuilder.build(fd, model, state, state0, dt);
        end
        
        function dt = getMaximumTimestep(fd, model, state, state0, dt, forces)
            dt = fd.FlowStateBuilder.getMaximumTimestep(fd, model, state, state0, dt, forces);
        end
        
        function [fd, state] = prepareTimestep(fd, model, state, state0, dt, drivingForces)
            % Called before each time-step solve
            [fd.FlowStateBuilder, state] = fd.FlowStateBuilder.prepareTimestep(fd, model, state, state0, dt, drivingForces);
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
