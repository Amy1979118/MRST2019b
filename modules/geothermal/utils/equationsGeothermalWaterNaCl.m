function [problem, state] = equationsGeothermalWaterNaCl( state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for the geothermal brine model.
%
% SYNOPSIS:
%   [problem, state] = equationsGeothermalWaterNaCl(state0, state, model, dt, drivingForces)
%
% DESCRIPTION:
%   This is the core function of the geothermal brine solver. This
%   function assembles the residual equations for the conservation of
%   water, NaCl and energy as well as required well equations. By default, 
%   Jacobians are also provided by the use of automatic differentiation.
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
%   'equationsGeothermalWater'

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false, ...
             'resOnly', false, ... 
             'iteration', -1);

         
opt = merge_options(opt, varargin{:});
assert(isempty(drivingForces.src))

op  = model.operators; 
G   = model.G; 
f   = model.fluid; 
r   = model.rock; 

% properties at current and previous timesteps
[p, c, T, wellSol]      = model.getProps(state, 'pressure','NaCl','temperature','wellSol');
[p0 , c0, T0, wellSol0] = model.getProps(state0, 'pressure','NaCl','temperature','wellSol');

% Initialization of independent variables
[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

if ~opt.resOnly,
    if ~opt.reverseMode,
        % define primary varible x and initialize
        [p, c, T, wellVars{:}]    = model.AutoDiffBackend.initVariablesAD(p, c, T, wellVars{:});       
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, c0, T0, wellVars0{:}] = model.AutoDiffBackend.initVariablesAD(p0, c0, T0, wellVars0{:}); 
       
    end
end

primaryVars = {'pressure','NaCl','temperature', wellVarNames{:}};

% Update state with AD-variables
state = model.setProps(state, {'pressure','NaCl', 'temperature'}, {p, c, T});
state0 = model.setProps(state0, {'pressure', 'NaCl', 'temperature'}, {p0, c0, T0});

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

% Salt Flux
rhoWvc = -op.faceUpstr(upcw,rhoW).*op.T_salt.*op.Grad(c);

% Thermal properties and fluxes 
[u, h, uR] = model.getProps(state, 'ThermalInternalEnergy', 'ThermalEnthalpy', 'ThermalInternalEnergyRock');
[u0, uR0]  = model.getProps(state0, 'ThermalInternalEnergy', 'ThermalInternalEnergyRock');
[vHr, vHf] = model.getProps(state,'HeatFluxRock','HeatFluxFluid');

uW  = u{1};
uW0 = u0{1};
hW  = h{1};

% output if needed
if model.outputFluxes
    state = model.storeFluxes(state,vW,[],[]);
end

if model.extraStateOutput
    state = model.storeDensity(state, rhoW, [], []);
    state = model.storeMobilities(state, mobW, [], []);
    state = model.storeUpstreamIndices(state, upcw, [], []);
end


% Equations 
% flux term 
rhoWvW   = op.faceUpstr(upcw,rhoW).*vW; 
rhoWvWhW = op.faceUpstr(upcw, hW).*rhoWvW;
crhoWvW  = op.faceUpstr(upcw,c).*rhoWvW; 

vol    =  model.G.cells.volumes; 

% Water accumulation term
water = (1/dt).*(pv.*rhoW - pv0.*rhoW0);  

% Salt accumulation term
conc  = (1/dt).*(pv.*rhoW.*c - pv0.*rhoW0.*c0);  

% Energy accumulation term (rock) 
energy = (1./dt).* ((vol-pv).*uR - (vol-pv0).*uR0);
% Energy accumulation term (fluid)
energy = energy + (1/dt).*(pv.*rhoW.*uW - pv0.*rhoW0.*uW0);

% correction equation for compatibility with BCs and wells
water = water./f.rhoWS;
conc  = conc./f.rhoWS; 

eqs         =  {water, conc, energy}; 
names       =  {'water','NaCl','temperature'};
types       =  {'cell','cell','cell'};
components  =  {c};

sW                 =  ones(model.G.cells.num,1);

% add bc and src terms (mass equation)
[eqs, state, src]  =  model.addBoundaryConditionsAndSources( eqs, names, types, state,{p}, {sW}, {mobW}, {rhoW}, {}, components, drivingForces);

% update thermal BCs
[ eqs ] = model.addThermalBoundaryConditionsGeothermal(eqs, names, p, T, drivingForces );
if (~isempty(src.bc.sourceCells))
    [ eqs ] = addEnthalpyContributionGeothermal(model, state, eqs, names, src,p, T, c, drivingForces);
end 

% well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, ...
                                                                p, {mobW}, {rhoW}, T, components, dt, opt); 

if ~isempty(drivingForces.W)
    W = drivingForces.W;
    for wNo = 1:numel(W)
        isInj = state.wellSol(wNo).qTs > 0;
        if isInj
            state.wellSol(wNo).T = W.T;
        else
            wc = W(wNo).cells;
            state.wellSol(wNo).T = mean(state.T(wc));
        end
        
    end
end

% add flux terms to equations 
eqs{1} = eqs{1} + op.Div(rhoWvW)./f.rhoWS;
eqs{2} = eqs{2} + (op.Div(rhoWvc) + op.Div(crhoWvW))./f.rhoWS;
eqs{3} = eqs{3} + op.Div(rhoWvWhW) + op.Div(vHf) + op.Div(vHr);

% add thermal source term to equations
[ eqs, srcT ] = model.addThermalSourcesGeothermal(eqs, names, T, drivingForces );

% equation scaling 
scaleMass = (dt./op.pv);
eqs{1} = eqs{1}.*scaleMass;

scaleConc = (dt./op.pv)./0.3;
eqs{2}    = eqs{2}.*scaleConc; 

uRRef     = mean(value(uR));
scaleEnergy = (dt./op.pv)./(uRRef);
if~isempty(srcT)
    scaleEnergy(srcT.cell) = 1./273.15; 
end
eqs{3} = eqs{3}.*scaleEnergy;

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
