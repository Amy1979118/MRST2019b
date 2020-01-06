%% Benchmark with TOUGH2
% In this example we simulate injection of pure water at 50C in a
% reservoir of brine at 10C. The same simulation was done with TOUGH2
% The results are compared to validate the implementation of equations in
% the geothermal module. 
%
% Setup: 2D grid without gravity -> equivalent to a top view 
% Pressure applied on the left side (injection), on two faces only
% Reservoir with constant pressure, temperature and NaCl mass fraction
% conditions (100 bars, 10C, 0.1)
% 
% Pressure, temperature and NaCl mass fraction are evaluated at two cells
% in the model: near the injection and in the model centre.

%% Add necessary MRST modules
mrstModule add ad-props  ad-core ad-blackoil geothermal mrst-gui

%% Set up the grid
physdim  = [100 200];              % domain size in x,y directions 
celldim  = [50 100];               % number of cells in x, y dims 
G = cartGrid(celldim, physdim);    % Cartesian Grid 
G = computeGeometry(G);            % grid geometry and connectivity 

%% Find injection faces and monitoring cell indices
% monitoring point near injection
Minj    = find(G.cells.centroids(:,1)==1 & G.cells.centroids(:,2)==101);   
% monitoring point in the center
Mcen    = find(G.cells.centroids(:,1)==51 & G.cells.centroids(:,2)==101);   

Finj1   = find(G.faces.centroids(:,1)==0 & G.faces.centroids(:,2)==99);
Finj2   = find(G.faces.centroids(:,1)==0 & G.faces.centroids(:,2)==101);
faces   = [Finj1; Finj2];

%% Plot the grid with the monitoring cells
plotting_grid = true; % set to true to plot the grid
if plotting_grid
figure()
plotGrid(G, 'faceColor', 'none')
hold on 
plotGrid(G, Minj, 'faceColor', 'green')
plotGrid(G, Mcen, 'faceColor', 'blue')
axis equal 
axis tight
title('Grid and location of monitoring cells')
txt1 = 'Injection';
txt2 = 'Centre';
text(5, 95, txt1,'FontSize',12)
text(55, 105, txt2,'FontSize',12)
end 

%% Set fluid structure properties
rhoWS               =   1000;
% define fluid structure
fluidBrine          =   initSimpleADIFluid('mu', 1.0e-3, ...
                                           'rho', rhoWS, ...
                                           'phases', 'W'); 
% add up thermal props                                       
fluidBrine          =   addThermalFluidProps(fluidBrine,'cp', 4.2e3, ...
                                             'lambdaF', 0.6, ...
                                             'useEOS', true, ...
                                             'brine', true, ...
                                             'dNaCl', 1e-6);

%% Make rock structure
perm            = 1e-14;
poro            = 0.1;  
% define rock structure
rock            = makeRock(G, perm, poro);    
% add thermal props
rock            = addThermalRockProps(rock, 'lambdaR', 2, 'rhoR', 2700, ...
                                      'CpR', 1000, 'tau', 1);     
%% Inititate state transient flow
gravity off 
clear wModelBrine
clear nonlinear
clear states; 

stateBrine0   = initResSol(G, 100*barsa, 1);  % pressure and saturation 
stateBrine0.T = ones(G.cells.num,1).*283.15;  % temperature
stateBrine0.c = ones(G.cells.num,1).*0.1;     % NaCl mass fraction    

%% Define the problem - GeothermalWaterNaclModel for brine problem
wModelBrine  = GeothermalWaterNaClModel(G, rock, fluidBrine, ...
                                        'extraStateOutput', true);

%% Set Boundary conditions and schedule

% bc: structure for the mass boundary conditions
% setPressureBC: function that merges pside and fluxeside. 
% set the pressure or flux condition on a given side of the model 
bc = []; 
bc = setPressureBC( bc, G, 'Right', ...         % Right side, all faces
                   'pressure', 100*barsa, ...
                   'sat',1); 
               
bc = addBC(bc, faces, ...                       % Left side, 2 faces
          'pressure',250*barsa,...
          'sat', 1);  
      
bc.c = [ones(100,1).*0.1; zeros(2,1)]; % salt contribution for BCs                                 
                                       % include right side (0.1) and 
                                       % injection (0): same order as for 
                                       % pressure bc
                                       
% bcT: structure for thermal boundary conditions 
% setThermalBC: analog to setPressureBC but for temperature
% set the temperature or heat flux condition on a given side of the model 
bcT = [];
bcT = setThermalBC(bcT, G, 'Right', ...         % Right side
                    'temperature', 283.15);  
               
bcT = addThermalBC(bcT, faces, ...              % Left side, 2 faces
                    'temperature', 323.15);                     

% Define the same timesteps as for the simulation with Tough2
data = load('tough2sim.mat'); tough2 = data.tough2;
time      = tough2(:,2);
timesteps = diff([0;time]);

scheduleBrine              = simpleSchedule(timesteps, 'bc', bc);
% add thermal BCs to schedule
scheduleBrine              = addThermalSchedule(scheduleBrine,'bcT', bcT);  

%% Run simulation
[~, states_brine]  = simulateScheduleAD(stateBrine0, wModelBrine, ...
                                        scheduleBrine);

%% Post processing
% Allocate memory for the different variables at the two monitoring cells 
it        = numel(states_brine);
Minj_b_p  = zeros(it,1);      % pressure at injection
Minj_b_t  = zeros(it,1);      % temperature at injection
Minj_b_c  = zeros(it,1);      % NaCl mass fraction at injection
Mcen_b_p  = zeros(it,1);      % pressure in the centre
Mcen_b_t  = zeros(it,1);      % temperature in the centre
Mcen_b_c  = zeros(it,1);      % NaCl mass fraction in the centre    

for i = 1:it     
    % Extract data from states
    Minj_b_p(i)  = states_brine{i}.pressure(Minj);
    Minj_b_t(i)  = states_brine{i}.T(Minj);
    Minj_b_c(i)  = states_brine{i}.c(Minj);
    Mcen_b_p(i)  = states_brine{i}.pressure(Mcen); 
    Mcen_b_t(i)  = states_brine{i}.T(Mcen); 
    Mcen_b_c(i)  = states_brine{i}.c(Mcen); 
end 

%% plot benchmark tough vs. MRST

tough2Ix    = [4,9,7,12,5,10];
tough2Data = tough2(:, tough2Ix);
tough2Data(:,1:2) = tough2Data(:,1:2)./barsa;
tough2Data(:,3:4) = tough2Data(:,3:4).*0.1;

mrstData   = [Minj_b_p/barsa , Mcen_b_p/barsa  , ...
              Minj_b_c       , Mcen_b_c        , ...
              Minj_b_t-273.15, Mcen_b_t-273.15];

plotTime   = cumsum(scheduleBrine.step.val)/day;

yl = {'Pressure [bar]', 'NaCl mass fraction', 'Temperature [C]'};
tl = {'Injection monitoring point', 'Centre monitoring point'};

figure('Position', [0,0,1200,800])
for i = 1:numel(tough2Ix)
    subplot(3,2,i), hold on
    plot(plotTime, tough2Data(:,i), '-' , 'LineWidth', 2);
    plot(plotTime, mrstData(:,i)  , '--', 'LineWidth', 5);
    hold off, box on
    ylabel(yl{ceil(i/2)});
    xlabel('Time [day]');
    title(tl{mod(i-1,2)+1})
    legend({'TOUGH2', 'MRST'},'Location','northeast')
end

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
