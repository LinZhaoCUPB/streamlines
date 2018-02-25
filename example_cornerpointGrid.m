%% Streamline Example for Corner-Point Grid
clear
mrstModule add ad-core ad-blackoil ad-props mrst-gui
close all
%% Reservoir geometry and petrophysical properties
% Generate a corner-point grid without faults
grdecl = simple_makeModel3([50, 50, 3], [1000, 1000, 12]*meter);
G = processGRDECL(grdecl);
G = computeGeometry(G);
G.type = [G.type, 'corner-point grid'];

px = 500*milli*darcy*ones(G.cells.num, 1);

perm = [px, px, px];
rock = makeRock(G, perm, 0.3);
%% Define wells and simulation schedule
simTime = 10*year;
nstep   = 25;
refine  = 5;

% Producers
W = verticalWell([], G, rock, 10, 10, 1,...
    'Name', 'P1', 'sign', -1, 'comp_i', [1, 1], 'Val', 250*barsa, 'Type', 'bhp', 'dir','z');
W = verticalWell(W,  G, rock, 10, 25, 2,...
    'Name', 'P2', 'sign', -1, 'comp_i', [1, 1], 'Val', 250*barsa, 'Type', 'bhp', 'dir','z');
W = verticalWell(W,  G, rock, 10, 40, 3,...
    'Name', 'P3', 'sign', -1, 'comp_i', [1, 1], 'Val', 250*barsa, 'Type', 'bhp', 'dir','z');

% Injectors, horizontal well
pv      = poreVolume(G, rock);
injRate = 1*sum(pv)/simTime;

indi  = 40;
indj  = 10:30;
indk  = 3;

[nx, ny] = deal(G.cartDims(1),G.cartDims(2));

wc = nx*ny*(indk-1) + nx*(indj-1) + indi;
wc = find(ismember(G.cells.indexMap, wc));

W = addWell(W, G, rock, wc, 'Name', 'I1',...
    'sign', 1,  'comp_i', [1, 0], 'Val', injRate, 'Type', 'rate', 'dir','y');

% Injectors, vertical well
W = verticalWell(W, G, rock, 40, 40, 2,...
    'Name', 'I2', 'sign', 1,  'comp_i', [1, 0], 'Val', injRate, 'Type', 'rate', 'dir','z');
        

figure
plotGrid(G,'facecol', 'c', 'edgecol', 'k', 'faceAlpha', 0.2)
plotGrid(G, vertcat(W.cells))
plotWell(G, W([W.sign] == -1), 'color', 'r')
plotWell(G, W([W.sign] ==  1), 'color', 'b')
view([30 40])
axis off

% Define the timesteps
timesteps = ones(10,1)*10*day;

% Set up the schedule containing both the wells and the timesteps
schedule = simpleSchedule(timesteps, 'W', W);
%% Set up simulation model
fluid = initSimpleADIFluid('mu',    [1, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);
                       
% Constant oil compressibility
c        = 0.001/barsa;
p_ref    = 300*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);

gravity reset on
model = TwoPhaseOilWaterModel(G, rock, fluid);
%% Define initial state
sW = ones(G.cells.num, 1);
sW(G.cells.centroids(:, 3) < 5) = 0;

sat = [sW, 1 - sW];

g = model.gravity(3);
% Compute initial pressure
p_res = p_ref + g*G.cells.centroids(:, 3).*...
   (sW.*model.fluid.rhoWS + (1 - sW).*model.fluid.rhoOS);
state0 = initResSol(G, p_res, sat);
%% Simulate base case
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule);
%% Generate streamlines
flux = states{end}.flux;
% Total face flux qt = qW + qO;
qt = sum(flux,2);
% Cell porosity
poro = pv./G.cells.volumes;
% Flux on each streamline
q0 = 5e-4;

[sc, sl, tof] = generate_streamlines(G, W, qt, q0, poro);
%%
% Plot streamlines
figure('color','w');hold on
plotWell(G, W([W.sign] == -1), 'color', 'r')
plotWell(G, W([W.sign] ==  1), 'color', 'b')

ind1 = ismember(sc(:,1), W(4).cells);
cellfun(@(x)plot3(x(:,1), x(:,2), x(:,3), 'color', 'm'), sl(ind1))

ind2 = ismember(sc(:,1), W(5).cells);
cellfun(@(x)plot3(x(:,1), x(:,2), x(:,3), 'color', 'g'), sl(ind2))

plotGrid(G, vertcat(W.cells))

view([30 40])
axis off
%%
