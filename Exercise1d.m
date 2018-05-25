clear all;
mrstModule add incomp

%% Define the model
% To set up a model, we need: a grid, rock properties (permeability), a
% fluid object with density and viscosity, and boundary conditions.
% gravity reset on
Nx = 40; Ny = 40;
G          = cartGrid([Nx, Ny], [1, 0.5]);
G          = computeGeometry(G);
rock       = makeRock(G, [1 1], 1);
fluid      = initSingleFluid('mu' , 1, 'rho', 1);


bc  = pside([], G, 'South', 1);
% bc  = pside(bc, G, 'North', 0);

xcell = G.cells.centroids(:,1);
ycell = G.cells.centroids(:,2);

% Add source term
f = 1*G.cells.volumes;;
src = addSource([],1:G.cells.num,f);

% Change the permeability tensor
eta = [0.1, 0.05, 0.025, 0.0125];
alpha = 0;
beta = 2;
% i=1;
for i=length(eta):-1:1
rock.perm = ones(G.cells.num,2);
reg = find(ycell>0.5-eta(i) & ycell<0.5+eta(i));
rock.perm(reg,1) = eta(i)^alpha;
rock.perm(reg,2) = eta(i)^beta;

%% Assemble and solve the linear system
% To solve the flow problem, we use the standard two-point
% flux-approximation method (TPFA), which for a Cartesian grid is the same
% as a classical seven-point finite-difference scheme for Poisson's
% equation. This is done in two steps: first we compute the
% transmissibilities and then we assemble and solve the corresponding
% discrete system.
T   = computeTrans(G, rock);
sol = incompTPFA(initResSol(G, 0.0), G, T, fluid, 'bc', bc, 'src', src);
press(i,:) = sol.pressure(1:Nx:Nx*Ny);
end
%% Plot the face pressures
% % clf
figure()
subplot(1,2,1)
plotCellData(G,sol.pressure);
% plotFaces(G, 1:G.faces.num, convertTo(sol.facePressure, barsa()));
% set(gca, 'ZDir', 'reverse'), title('Pressure [bar]')
view(2), colorbar
% set(gca,'DataAspect',[1 1 10]);
subplot(1,2,2)
% figure()
plot(ycell(1:Ny:Ny*Nx),press);