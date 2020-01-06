%% Subset of SEP10 Model 2
% In this example, we simulate injection of hot water into Layer 10 of
% SPE10 Model 2, once with thermally insulated boundary conditions, and
% once with fixed-temperature boundary conditions.

%% Add modules
mrstModule add geothermal upr vem vemmech spe10 mrst-gui ...
    ad-props ad-core ad-blackoil

%% Set up fine-scale model
% We extract the lower half of layer 13 of SPE10 2
[state0Ref, modelRef, scheduleRef] = setupSPE10_AD('layers', 10);
GRef    = modelRef.G;
rockRef = modelRef.rock;
WRef    = scheduleRef.control(1).W;

%% Make PEBI grid
% We use the upr module to construct a PEBI grid with refinement around the
% wells
rng(2019) % For reproducibility
n = 40;   % Approx number of cells in each direction
% Get well coordinates
l = max(GRef.nodes.coords(:,1:2));
wellLines = mat2cell(GRef.cells.centroids(vertcat(WRef.cells),1:2), ...
                                                    ones(numel(WRef),1), 2)';
% Construct PEBI grid
G = pebiGrid(max(l)/n, l, 'wellLines'     , wellLines, ... % Well coords
                          'wellRefinement', true     , ... % Refine
                          'wellGridFactor', 0.4      );
G = computeGeometry(G);       % Compute geometry

%% Sample rock properties
% We assign rock properties in the PEBI grid cells by sampling from the
% fine grid using sampleFromBox
poro = sampleFromBox(G, reshape(rockRef.poro, GRef.cartDims));
perm = zeros(G.cells.num,G.griddim);
for i = 1:G.griddim
    perm(:,i) = sampleFromBox(G, reshape(rockRef.perm(:,i), GRef.cartDims));
end
rock = makeRock(G, perm, poro);
rock = addThermalRockProps(rock);

%% Set up schedule and initial state
W        = WRef;
schedule = scheduleRef;
x  = G.cells.centroids;
xwR = GRef.cells.centroids(vertcat(WRef.cells),1:2);
% Slick oneliner to find corresponding cells in the new grid
[~, c] = min(sum(bsxfun(@minus, reshape(xwR, [], 1 , G.griddim), ...
                             reshape(x , 1 , [], G.griddim)).^2,3), [], 2);
K0 = 273.15*Kelvin;
Tinj = K0 + 100*Kelvin;
for i = 1:numel(W)
    W(i).cells = c(i);
    W(i).T = Tinj;
    if strcmpi(W(i).type, 'rate')
        W(i).val = 5*W(i).val; % Increase injection rate
    end
end
xw = G.cells.centroids(vertcat(W.cells),1:2);
schedule.control.W = W;
% Set initial state
state0   = initResSol(G, state0Ref.pressure(1), 1);
Tres     = (273.15 + 30)*Kelvin;
state0.T = repmat(Tres, G.cells.num, 1);

%% Inspect the geological model
figure('Position', [0,0,500,800])
K = convertTo(rock.perm(:,1),milli*darcy);
plotCellData(G,log10(K),'edgeAlpha',.1); axis tight;
set(gca,'FontSize',12)
mrstColorbar(K,'South',true); axis equal tight
hold on; plot(xw(:,1),xw(:,2),'.r','MarkerSize',18); hold off, drawnow

%% Create model
% We use a single-phase geothermal model. Viscosity and density are
% p/T-dependent, and will be set later
fluid = initSimpleADIFluid('phases', 'W', 'n', 1, 'mu', 1, 'rho', 1);
% Assign thermal properties of the fluid, with equation of state from
% Spivey et. al (2004)
fluid = addThermalFluidProps(fluid, 'useEOS', true);
% Add thermal rock properties, with high thermal conductivity and low
% specific heat capacity
Watt  = joule/second;
rock  = addThermalRockProps(rock, 'lamdaR', 100*Watt/(meter*Kelvin), ...
                                  'CpR'   , 250*joule/(kilo*gram*Kelvin));
model = GeothermalWaterModel(G, rock, fluid); % Make model
% The EOS is valid for pressure/temperature within a given range. We
% provide these to the model so that pressure/temperature are within these
% during the nonlinear solution step
model.maximumPressure    = 200e6;           % Maximum pressure
model.minimumTemperature = K0;              % Minimum temperature 
model.maximumTemperature = K0 + 275*Kelvin; % Maximum temperature
model.extraStateOutput   = true;     % Output density and mobility to state

%% Simulate
% No boundary conditions means a completely insulated, closed flow
% compartment
[wellSols, states, reports] = simulateScheduleAD(state0, model, schedule);

%% Interactive plot of the results
figure('Position', [0,0,500,800], 'Name', 'Insulated BCs')
plotToolbar(model.G, states); colormap(hot); axis equal tight
set(gca,'FontSize',12);
plotWellSols(wellSols, schedule.step.val);

%% Simulate with fixed-temperature boundary conditions
bcT = addThermalBC([], boundaryFaces(G), 'temperature', Tres);
schedule.control(1).bcT = bcT;
[wellSolsBC, statesBC, reportsBC] = simulateScheduleAD(state0, model, schedule);

%% Compare with insulated reservoir results
figure('Position', [0,0,500,800], 'Name', 'Fixed-temperature BCs')
plotToolbar(model.G, statesBC); colormap(hot); axis equal tight
set(gca,'FontSize',12)
plotWellSols({wellSols, wellSolsBC}, schedule.step.val, ...
    'dataSetNames', {'Insulated', 'Fixed-temp'});