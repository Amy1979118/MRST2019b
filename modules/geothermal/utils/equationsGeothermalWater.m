function [ problem, state ] = equationsGeothermalWater( state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for the geothermal pure water model.
%
% SYNOPSIS:
%   [problem, state] = equationsGeothermalWater(state0, state, model, dt, drivingForces)
%
% DESCRIPTION:
%   This is the core function of the geothermal pure water solver. This
%   function assembles the residual equations for the conservation of water
%   and energy as well as required well equations. By default, Jacobians are
%   also provided by the use of automatic differentiation.
%
% REQUIRED PARAMETERS:
%   state0    - Reservoir state at the previous timestep. Assumed to have
%               physically reasonable values.
%
%   state     - State at the current nonlinear iteration. The values do not
%               need to be physically reasonable.
%
%   model     - GeothermalModel-derived class. Typically,
%               equationsGeothermalWater will be called from the class
%               getEquations member function.
%
%   dt        - Scalar timestep in seconds.
%
%   drivingForces - Struct with fields:
%                   * W for wells. Can be empty for no wells.
%                   * bc for boundary conditions. Can be empty for no bc.
%                   * src for source terms. Can be empty for no sources.
%
% OPTIONAL PARAMETERS:
%   'Verbose'    -  Extra output if requested.
%
%   'reverseMode'- Boolean indicating if we are in reverse mode, i.e.
%                  solving the adjoint equations. Defaults to false.
%
%   'resOnly'    - Only assemble residual equations, do not assemble the
%                  Jacobians. Can save some assembly time if only the
%                  values are required.
%
%   'iterations' - Nonlinear iteration number. Special logic happens in the
%                  wells if it is the first iteration.
% RETURNS:
%   problem - LinearizedProblemAD class instance, containing the water
%             and energy conservation equations, as well as well equations
%             specified by the WellModel class.
%
%   state   - Updated state. Primarily returned to handle changing well
%             controls from the well model.
%
% SEE ALSO:
%   'equationsGeothermalWaterNaCl'
% 
% Equations: 
    % {eq1} conservation of mass
    % phi * dP/dt + Div(perm/mu(Grad(P)-rho_f*g)) = 0
    % {eq2} conservation of heat 
    % rho_b*c_b * dT/dt + rho_f*c_f * Div(uT) - Div(lambda_r*Grad(T)) = 0

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false, ...
             'resOnly', false, ... 
             'iteration', -1);

         
opt = merge_options(opt, varargin{:});
assert(isempty(drivingForces.src))

op  = model.operators; 
f   = model.fluid; 
r   = model.rock; 

% get properties for current and previous timesteps
[p, T, wellSol]     = model.getProps(state, 'pressure','temperature', 'wellSol');
[p0 , T0, wellSol0] = model.getProps(state0, 'pressure','temperature','wellSol');

% Initialization of independent variables
[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

if ~opt.resOnly,
    if ~opt.reverseMode,
        % define primary varible x and initialize
        [p, T, wellVars{:}]    =  model.AutoDiffBackend.initVariablesAD(p, T, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, T0, wellVars0{:}] =  model.AutoDiffBackend.initVariablesAD(p0, T0, wellVars0{:}); 
    end
end
primaryVars = {'pressure','temperature', wellVarNames{:}};

% Update state with AD-variables
state = model.setProps(state, {'pressure', 'temperature'}, {p, T});
state0 = model.setProps(state0, {'pressure', 'temperature'}, {p0, T0});

% Water properties and fluxes 
[rho, mob]          = model.getProps(state, 'Density', 'Mobility');
[phaseFlux, flags]  = model.getProps(state, 'PhaseFlux', 'PhaseUpwindFlag');
pv                  = model.getProp(state, 'PoreVolume');
pv0                 = model.getProp(state0, 'PoreVolume');
rho0                = model.getProp(state0, 'Density');

rhoW  = rho{1};
rhoW0 = rho0{1};
mobW  = mob{1};
vW    = phaseFlux{1};
upcw  = flags{1};

% Thermal properties and fluxes 
[u, h, uR] = model.getProps(state, 'ThermalInternalEnergy', 'ThermalEnthalpy', 'ThermalInternalEnergyRock');
[u0, uR0]  = model.getProps(state0, 'ThermalInternalEnergy', 'ThermalInternalEnergyRock');
[vHr, vHf] = model.getProps(state,'HeatFluxRock','HeatFluxFluid');

uW  = u{1};
uW0 = u0{1};
hW  = h{1};

% Output if needed
if model.outputFluxes
    state = model.storeFluxes(state,vW,[],[]);
end 
 

if model.extraStateOutput
    state = model.storeDensity(state, rhoW, [], []);
    state = model.storeMobilities(state, mobW, [], []);
    state = model.storeUpstreamIndices(state, upcw, [], []);
end

% Equations new implementations
% fluxes
rhoWvW   = op.faceUpstr(upcw,rhoW).*vW;
rhoWvWhW = op.faceUpstr(upcw, hW).*rhoWvW;

vol    =  model.G.cells.volumes; 

% Water accumulation term
water  = (1/dt).*(pv.*rhoW - pv0.*rhoW0);

% Energy accumulation term (rock)
energy = (1./dt).* ((vol-pv).*uR - (vol-pv0).*uR0);
% Energy accumulation term (fluid)
energy = energy + (1/dt).*(pv.*rhoW.*uW - pv0.*rhoW0.*uW0);

% correction equation for compatibility with BCs and wells
water = water./f.rhoWS;
eqs     =   {water, energy};
names   =   {'water','temperature'};
types   =   {'cell','cell'};

sW      =   ones(model.G.cells.num,1);

% add bc and src terms (mass equation)
[eqs, state, src] = model.addBoundaryConditionsAndSources(eqs, names, types, state,{p}, {sW}, {mobW}, {rhoW}, {}, {}, drivingForces);

% update thermal BCs
[ eqs ] = model.addThermalBoundaryConditionsGeothermal(eqs, names, p, T, drivingForces );
if (~isempty(src.bc.sourceCells))
    [ eqs ] = addEnthalpyContributionGeothermal(model, state, eqs, names, src,p, T, 0, drivingForces);
end 

% well equations and scaling 
components = {};
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, ...
                                                                p, {mobW}, {rhoW}, T, components, dt, opt); 

tol = 1e-14;
if ~isempty(drivingForces.W)
    W = drivingForces.W;
    for wNo = 1:numel(W)
        isInj = state.wellSol(wNo).qTs > tol;
        if isInj
            state.wellSol(wNo).T = W(wNo).T;
        else
            wc = W(wNo).cells;
            state.wellSol(wNo).T = mean(value(state.T(wc)));
        end
    end
end

% add flux terms to equations 
eqs{1} = eqs{1} + op.Div(rhoWvW)./f.rhoWS;
eqs{2} = eqs{2} + op.Div(rhoWvWhW) + op.Div(vHf) + op.Div(vHr);

% add thermal source term to equations
[ eqs, srcT ] = model.addThermalSourcesGeothermal(eqs, names, T, drivingForces );

% equation scaling 
scaleMass = (dt./op.pv);
eqs{1} = eqs{1}.*scaleMass;

uRRef = mean(value(uR));
scaleEnergy = (dt./op.pv)./(uRRef);
if~isempty(srcT)
    scaleEnergy(srcT.cell) = 1./273.15; 
end
eqs{2} = eqs{2}.*scaleEnergy;
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
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