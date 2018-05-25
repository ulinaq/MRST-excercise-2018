clear all;
mrstModule add incomp

%% Define the model
% To set up a model, we need: a grid, rock properties (permeability), a
% fluid object with density and viscosity, and boundary conditions.
% gravity reset on
Nx = 80; Ny = 80;
G          = cartGrid([Nx, Ny], [1, 1]);
G          = computeGeometry(G);
rock       = makeRock(G, 1, 1);
fluid      = initSingleFluid('mu' , 1, 'rho', 1);

n = 10;
% bc  = pside([], G, 'South', 1);
% bc  = pside(bc, G, 'North', 0);
bc  = pside([], G, 'West', 1);
bc  = pside(bc, G, 'East', 0);

xcell = G.cells.centroids(:,1);
ycell = G.cells.centroids(:,2);

% Change the permeability tensor
K = [10 1];
Upreg = zeros(Nx*Ny,1);
Loreg = zeros(Nx*Ny,1);
for i=1:n
Upreg = Upreg + ycell>=(i-1)/n & ycell<(i-1)/n+0.25/n;
Loreg = Loreg + ycell>=(i-1)/n+0.25/n & ycell<i/n;
end
rock.perm(Upreg) = K(1);
rock.perm(Loreg) = K(2);

%% Assemble and solve the linear system
% To solve the flow problem, we use the standard two-point
% flux-approximation method (TPFA), which for a Cartesian grid is the same
% as a classical seven-point finite-difference scheme for Poisson's
% equation. This is done in two steps: first we compute the
% transmissibilities and then we assemble and solve the corresponding
% discrete system.
T   = computeTrans(G, rock);
sol = incompTPFA(initResSol(G, 0.0), G, T, fluid, 'bc', bc);

%% Plot the face pressures
% % clf
figure()
% plotCellData(G,sol.pressure);
plotCellData(G,rock.perm);
% % plotFaceData(G,sol.flux)
% plotFaces(G, 1:G.faces.num, convertTo(sol.facePressure, barsa()));
% % set(gca, 'ZDir', 'reverse'), title('Pressure [bar]')
view(2), colorbar
% % set(gca,'DataAspect',[1 1 10]);
% figure()
% subplot(1,2,1)
% plot(xcell(1:Ny),sol.pressure(1:Nx:Nx*Ny));
% subplot(1,2,2)
% plot(xcell(1:Ny),sol.pressure(1:Nx));
fx=reshape(sol.flux(1:(Nx+1)*Ny),Nx+1,Ny);
fy=reshape(sol.flux((Nx+1)*Ny+1:end),Nx,Ny+1);
as=sum(fy,1);
ax=sum(fx,2);