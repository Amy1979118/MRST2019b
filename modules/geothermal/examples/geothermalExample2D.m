%% Simple 2D case with gravity and variable fluid properties 
% Simulate mass and heat transfer in a reservoir. The model considers
% injection at the left boundary and outflow at the right boundary
% (pressure set at reservoir condition). The top and bottom layer have
% no-flow boundary conditions. The initial reservoir is set at 100 bars and
% 10 C. Pure water is injected at 50 C.

%% Add necessary MRST modules
mrstModule add ad-props  ad-core ad-blackoil geothermal mrst-gui

%% Set plot options
setAxProps = @(ax) set(ax, 'View'              , [-140,20]    , ...
                           'PlotBoxAspectRatio', [2,2,1]      , ...
                           'Projection'        , 'Perspective', ...
                           'Box'               , 'on'         , ...
                           'XLimSpec'          , 'tight'      , ...
                           'YLimSpec'          , 'tight'      , ...
                           'ZLimSpec'          , 'tight'      );

%% Set up the grid
physdim  = [100 50];               % domain size in x, y directions
celldim  = [50  50];               % number of cells in x, y directions
G = cartGrid(celldim, physdim);    % Cartesian grid
G = computeGeometry(G);            % grid geometry and connectivity

%% Make fluid structure properties
rhoWS          =   1000;
% Define fluid structure
fluid          =   initSimpleADIFluid('mu', 1.0e-3, 'rho', rhoWS, ...
                                        'phases', 'W');         
% Add thermal properties. The module supports using an EOS from Spivey et
% al (2004) with p/T-dependent density and viscosity. For convenience, this
% can be added to the fluid by setting the optional argument 'useEOS' to
% true in 'addThermalFluidProps'.
fluid          =   addThermalFluidProps(fluid, 'cp'     , 4.2e3, ...
                                               'lambdaF', 0.6  , ...
                                               'useEOS' , true , ...
                                               'brine'  , false);

%% Inspect fluid model
% We plot the fluid density and viscosity as a function of pressure and
% temperature. We also indicate the normal operational range by a black
% triangle
K0   = 273.15*Kelvin;
pMin = 1e6*Pascal;      % Minimum pressure
pMax = 200e6*Pascal;    % Maximum pressure
TMin = K0;              % Minimum temperature 
TMax = K0 + 275*Kelvin; % Maximum temperature 
n    = 100;
% Get pressure/temperature grid
[p,T] = ndgrid(linspace(pMin, pMax, n), linspace(TMin, TMax, n));
% "Normal" operational range
p1    = linspace(pMin, 500*barsa, n)';
T1    = ones(n, 1)*(K0 + 100*Kelvin);
p2    = ones(n, 1)*500*barsa;
T2    = linspace(TMin + 3*Kelvin, K0 + 100*Kelvin, n)';
% Plot density
rhoW = reshape(fluid.rhoW(p(:),T(:)), n, n);
figure(), surf(p/barsa, T-K0, rhoW, 'EdgeColor', 'none'); setAxProps(gca);
xlabel('Pressure (bar)'), ylabel('Temperature (C)'), zlabel('Density, (kg/m^3)');
hold on
rhoW1 = fluid.rhoW(p1, T1);
plot3(p1/barsa, T1-K0, rhoW1, 'k');
rhoW2 = fluid.rhoW(p2, T2);
plot3(p2/barsa, T2-K0, rhoW2, 'k');
hold off
% Plot viscosity
muW = reshape(fluid.muW(p(:),T(:)), n, n);
figure(), surf(p/barsa, T-K0, muW/(centi*poise), 'EdgeColor', 'none'); setAxProps(gca);
xlabel('Pressure (bar)'), ylabel('Temperature (C)'), zlabel('Viscosity, (cp)');
hold on
muW1 = fluid.muW(p1, T1);
plot3(p1/barsa, T1-K0, muW1/(centi*poise), 'k');
muW2 = fluid.muW(p2, T2);
plot3(p2/barsa, T2-K0, muW2/(centi*poise), 'k');
hold off

%% Make rock structure
perm            = 1e-14;
poro            = 0.1;  
% Define rock structure
rock            = makeRock(G, perm, poro);   
% Add thermal props
rock            = addThermalRockProps(rock, 'lambdaR', 2, 'rhoR', 2700, ...
                                     'CpR', 1000);     

%% Inititate state transient flow
state0   = initResSol(G, 100*barsa, 1);       % Pressure and saturation 
state0.T = ones(G.cells.num,1).*(273.15+10);  % Temperature

%% Define the problem - GeothermalWaterModel for pure water problem
% Define the gravity vector for a 2D x,y case using MRST convention
gravity reset on
gravity([0 -9.81]); % gravity is downward 
wModel  =  GeothermalWaterModel(G, rock, fluid, 'extraStateOutput', true);
% Define limits of temperature validity for EOS
wModel.minimumTemperature = TMin; 
wModel.maximumTemperature = TMax;
wModel.maximumPressure    = pMax;

%% Set up boundary conditions
pRes = 100*barsa;
pInj = pRes + 100*barsa;
Tinj = K0+50; 
TRes = K0+10;
[bc, bcT] = deal([]);
% Flow BCs
bc   = pside(bc, G, 'left', pInj, 'sat', 1);
bc   = pside(bc, G, 'right', pRes, 'sat', 1);
% Thermal BCs
bcT  = setThermalBC(bcT, G, 'left', 'temperature', Tinj);
bcT  = setThermalBC(bcT, G, 'right', 'temperature', TRes);

%% Set up schedule
timesteps =  rampupTimesteps(1*year, 30*day, 8); 
schedule  =  simpleSchedule(timesteps, 'bc', bc);
% Add thermal BCs to schedule
schedule  =  addThermalSchedule(schedule,'bcT', bcT);  

%% Run simulation
[~, states]  =  simulateScheduleAD(state0, wModel, schedule);

%% Visualization
figure(), plotToolbar(G, states); axis tight equal, colormap(hot)

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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