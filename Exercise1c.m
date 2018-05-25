mrstModule add incomp

%% Define the model
% To set up a model, we need: a grid, rock properties (permeability), a
% fluid object with density and viscosity, and boundary conditions.
% gravity reset on
for i=1:5
N = 5*2^(i-1);
G          = cartGrid([N, N], [1, 1]);
G          = computeGeometry(G);
rock       = makeRock(G, [1 1], 1);
fluid      = initSingleFluid('mu' , 1, 'rho', 1);
bc  = pside([], G, 'West', 1);
bc  = pside(bc, G, 'East', 1);
bc  = pside(bc, G, 'North', 1);
bc  = pside(bc, G, 'South', 1);

xcell = G.cells.centroids(:,1);
ycell = G.cells.centroids(:,2);
% Change the permeability tensor
% rock.perm = [1+xcell.^2,ones(G.cells.num,1)];
rock.perm = [1+xcell.^2,ycell];
% Add source term
% f = [2*ycell.*(1-ycell)+2*xcell.*(1-xcell)].*G.cells.volumes;
% f = [(6*ycell-2).*xcell.^2.*(1-xcell)+(12*xcell.^3-6*xcell.^2+6*xcell-2).*ycell.^2.*(1-ycell)].*G.cells.volumes;
f = [(9*ycell.^2-4*ycell).*xcell.^2.*(1-xcell)+(12*xcell.^3-6*xcell.^2+6*xcell-2).*ycell.^2.*(1-ycell)].*G.cells.volumes;
src = addSource([],1:G.cells.num,f);

%% Assemble and solve the linear system
% To solve the flow problem, we use the standard two-point
% flux-approximation method (TPFA), which for a Cartesian grid is the same
% as a classical seven-point finite-difference scheme for Poisson's
% equation. This is done in two steps: first we compute the
% transmissibilities and then we assemble and solve the corresponding
% discrete system.
T   = computeTrans(G, rock);
sol = incompTPFA(initResSol(G, 0.0), G, T, fluid, 'bc', bc, 'src', src);

% pressureexact = ycell.*(1-ycell).*xcell.*(1-xcell);  %exc 2
pressureexact = ycell.^2.*(1-ycell).*xcell.^2.*(1-xcell)+1;
pressureerror = abs(pressureexact - sol.pressure);
er(i) = sqrt(sum(pressureerror.^2.*G.cells.volumes));
en(i) = sqrt(sum(pressureexact.^2.*G.cells.volumes));
ts(i)=1/N;
end
re=er./en;
da = log(re(1:end-1)./re(2:end))./log(ts(1:end-1)./ts(2:end))
%% Plot the face pressures
% clf
figure()
plotCellData(G,sol.pressure);
% plotFaces(G, 1:G.faces.num, convertTo(sol.facePressure, barsa()));
% set(gca, 'ZDir', 'reverse'), title('Pressure [bar]')
view(2), colorbar
% set(gca,'DataAspect',[1 1 10]);
