mrstModule add incomp

%% Define the model
% To set up a model, we need: a grid, rock properties (permeability), a
% fluid object with density and viscosity, and boundary conditions.
% gravity reset on
G          = cartGrid([10, 10], [1, 1]);
G          = computeGeometry(G);
rock       = makeRock(G, 1, 1);
fluid      = initSingleFluid('mu' ,    1,  'rho', 1);
bc  = pside([], G, 'West', 1);
bc  = pside(bc, G, 'East', 0);

xcell = G.cells.centroids(:,1);
ycell = G.cells.centroids(:,2);

%% Assemble and solve the linear system
% To solve the flow problem, we use the standard two-point
% flux-approximation method (TPFA), which for a Cartesian grid is the same
% as a classical seven-point finite-difference scheme for Poisson's
% equation. This is done in two steps: first we compute the
% transmissibilities and then we assemble and solve the corresponding
% discrete system.
T   = computeTrans(G, rock);
sol = incompTPFA(initResSol(G, 0.0), G, T, fluid, 'bc', bc);
pressureexact = 1-xcell;  %exc 1
pressureerror = abs(pressureexact - sol.pressure);
er = sqrt(sum(pressureerror.^2.*G.cells.volumes))
norm(pressureerror)
en = sqrt(sum(pressureexact.^2.*G.cells.volumes));
er/en
norm(pressureerror)/norm(pressureexact)

%% Plot the face pressures
% % clf
figure()
plotCellData(G,sol.pressure);
% plotFaces(G, 1:G.faces.num, convertTo(sol.facePressure, barsa()));
% set(gca, 'ZDir', 'reverse'), title('Pressure [bar]')
view(2), colorbar
% set(gca,'DataAspect',[1 1 10]);
