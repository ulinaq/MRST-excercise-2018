mrstModule add mimetic incomp

nx = 40; ny = 40; nz = 1;

G         = cartGrid([nx ny nz]);
G         = computeGeometry(G);
rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]);
rock.poro = repmat(0.3            , [G.cells.num, 1]);

x = linspace(0, 1, 11).';
y = linspace(1, 0, 11).';
pc_form = 'nonwetting';
cap_scale = 10;
[kr, pc]  = tabulatedSatFunc([x, x.^2, y.^2, y.*cap_scale*barsa]);
props = constantProperties([   1,  10] .* centi*poise, ...
                           [1000, 700] .* kilogram/meter^3);
fluid = struct('properties', props                  , ...
               'saturation', @(x, varargin)    x.s  , ...
               'relperm'   , kr);
fluid_pc = struct('properties', props                  , ...
                  'saturation', @(x, varargin)    x.s  , ...
                  'relperm'   , kr                     , ...
                  'pc'        , @(x, varargin) pc(x.s));
xDummy   = initState(G, [], [0, 1]);
xDummy.s = linspace(0, 1, numel(xDummy.s))'; ...
pc = convertTo(fluid_pc.pc(xDummy), barsa);

clf
plot(xDummy.s, pc);
xlabel('s_w'); ylabel('pc [bar]');
title('Capillary pressure curve')

rate = 0.5*meter^3/day;
bhp  = 1*barsa;

W = verticalWell([], G, rock, 1, 1, 1:nz,          ...
                 'Type', 'rate', 'Val', rate, ...
                 'Radius', .1, 'Name', 'I', 'Comp_i', [1 0]);
W = verticalWell(W, G, rock, nx, ny, 1:nz,     ...
                 'Type','bhp', 'Val', bhp, ...
                 'Radius', .1, 'Dir', 'x', 'Name', 'P', 'Comp_i', [0 1]);
rSol    = initState(G, W, 0, [0.2, 0.8]);
rSol_pc = initState(G, W, 0, [0.2, 0.8]);

gravity off
verbose = false;

S  = computeMimeticIP(G, rock, 'Verbose', verbose,'InnerProduct','ip_tpf');
psolve  = @(state, fluid) incompMimetic(state, G, S, fluid, 'wells', W);
tsolve  = @(state, dT, fluid) implicitTransport(state, G, dT, rock, ...
                                                fluid, 'wells', W, ...
                                                'verbose', verbose);
% tsolve = @(state, dT) explicitTransport(state, G, dT, rock, fluid, ...
%                                        'wells', W, 'verbose', verbose);
rSol    = psolve(rSol, fluid);
rSol_pc = psolve(rSol_pc, fluid_pc);
T      = 300*day();
dT     = T/15;
dTplot = 100*day();  % plot only every 100th day
N      = fix(T/dTplot);
pv     = poreVolume(G,rock);
t  = 0; plotNo = 1;
h1 = 'No pc - '; H2 = 'Linear pc - ';
e = []; p_org = []; p_pc = [];
figure;

while t < T,
   % TRANSPORT SOLVE
   rSol    = tsolve(rSol, dT, fluid);
   rSol_pc = tsolve(rSol_pc, dT, fluid_pc);

   % Check for inconsistent saturations
   s = [rSol.s(:,1); rSol_pc.s(:,1)];
   assert(max(s) < 1+eps && min(s) > -eps);

   % Update solution of pressure equation.
   rSol    = psolve(rSol,    fluid);
   rSol_pc = psolve(rSol_pc, fluid_pc);

   % Measure water saturation in production cells in saturation
   e = [e; sum(abs(rSol.s(:,1) - rSol_pc.s(:,1)).*pv)/sum(pv)]; %#ok
   p_org = [p_org; rSol.s(W(2).cells,1)' ];                 %#ok
   p_pc  = [p_pc; rSol_pc.s(W(2).cells,1)'];                 %#ok

   % Increase time and continue if we do not want to plot saturations
   t = t + dT;
   if ( t < plotNo*dTplot && t<T) continue, end

   % Plot saturation
   heading = [num2str(convertTo(t,day)),  ' days'];
   r = 0.01;
   subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.48]), cla
   plotCellData(G, rSol.s(:,1));
   caxis([0 1]), view(60,50), axis equal off, title([h1 heading])

   subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.48]), cla
   plotCellData(G, rSol_pc.s(:,1));
   caxis([0 1]), view(60,50), axis equal off, title([H2 heading])

   plotNo = plotNo+1;

end
figure()
n = numel(p_org(:,1));
plot(1:n,p_org(:,1),'-o',1:n,p_pc(:,1),'--*')
legend('No capillary pressure','Linear capillary pressure','Location','Best');
title('Water breakthrough at heel');