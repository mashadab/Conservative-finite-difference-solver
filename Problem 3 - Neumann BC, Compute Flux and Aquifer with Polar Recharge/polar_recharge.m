% Steady confined aquifer with polar recharge
% date: 03/03/2021
% author: Mohammad Afzal Shadab
h_ana = @(x) 1-x;
q_ana = @(x) 1 + x*0;
x_ana = linspace(0,1,1e2);

%% Build grid and ops
Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 20;
Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid);
L = - D*G;
fs = zeros(Grid.N,1);

%% Build BCs
BC.dof_dir   = Grid.dof_xmax;
BC.dof_f_dir = Grid.dof_f_xmax;
BC.g         = h_ana(Grid.xc(BC.dof_dir));

BC.dof_neu   = Grid.dof_xmin;
BC.dof_f_neu = Grid.dof_f_xmin;
BC.qb        = 1.0;
[B,N,fn] = build_bnd(BC,Grid,I);

%% Solve for ha and q
h = solve_lbvp(L,fs+fn,B,BC.g,N);
q = comp_flux(D,1,G,h,fs,Grid,BC);

%% Plot solution
subplot 121
plot(x_ana,h_ana(x_ana),'-'), hold on
plot(Grid.xc,h,'o','markersize',6,'markerfacecolor','w')
xlabel 'x'' '
ylabel 'h'' '
ylim([0 1.2])
pbaspect([1 .8 1])

subplot 122
plot(x_ana,q_ana(x_ana),'-'), hold on
plot(Grid.xf,q,'o','markersize',6,'markerfacecolor','w')
xlabel 'x'' '
ylabel 'q'' '
ylim([0 1.2])
pbaspect([1 .8 1])