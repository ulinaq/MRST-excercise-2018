mrstModule add incomp

%% Define the model
% To set up a model, we need: a grid, rock properties (permeability), a
% fluid object with density and viscosity, and boundary conditions.
% gravity reset on
Nx=125;
Ny=Nx;
G          = cartGrid([Nx, Ny], [1, 1]);
G          = computeGeometry(G);
rock       = makeRock(G, 1, 1);
fluid      = initSingleFluid('mu' ,    1, ...
                             'rho', 1);
bc  = pside([], G, 'West', 1);
bc  = pside(bc, G, 'East', 0);


xcell = G.cells.centroids(:,1);
ycell = G.cells.centroids(:,2);
% Change the permeability tensor

delta = 0.05;
BcellInd = find(xcell<=0.5+delta & xcell>=0.5-delta & ycell>=0.5-delta & ycell<=0.5+delta);
perm = [0.001 0.01 0.1 1 10 100 1000];
fy=find(ycell>0.5-0.5/Ny & ycell<0.5+0.5/Ny);
% i=7;
for i = 1:length(perm)
rock.perm(BcellInd) = perm(i);
%% Assemble and solve the linear system
% To solve the flow problem, we use the standard two-point
% flux-approximation method (TPFA), which for a Cartesian grid is the same
% as a classical seven-point finite-difference scheme for Poisson's
% equation. This is done in two steps: first we compute the
% transmissibilities and then we assemble and solve the corresponding
% discrete system.
T   = computeTrans(G, rock);
sol = incompTPFA(initResSol(G, 0.0), G, T, fluid, 'bc', bc);
pressurenum{i} = sol.pressure(fy);
end
% pre = reshape(sol.pressure,Nx,Ny);
fx=xcell(fy);
f1=find(ycell>0.1-0.5/Ny & ycell<0.1+0.5/Ny);
%% Plot the face pressures
% % clf
% figure()
% subplot(1,2,1)
% plotCellData(G,as); title('Pressure (K_{\Omega_b}=0.01)');
% view(2), colorbar
% subplot(1,2,2)
% plotCellData(G,sol.pressure); title('Pressure (K_{\Omega_b}=100)');
% plotFaces(G, 1:G.faces.num, convertTo(sol.facePressure, barsa()));
% set(gca, 'ZDir', 'reverse'), title('Pressure [bar]')
% view(2), colorbar
% set(gca,'DataAspect',[1 1 10]);
figure()
% plot(fx,pre(:,[1 7 9 10 11 13]))
% plot(fx,pre(:,[1 20 40 70 80 90 100 113]))
% plot(fx,sol.pressure(f1),fx,sol.pressure(f3),fx,sol.pressure(fy));
plot(fx,pressurenum{1},fx,pressurenum{2});
hold on
for i=3:length(perm)
    plot(fx,pressurenum{i})
end