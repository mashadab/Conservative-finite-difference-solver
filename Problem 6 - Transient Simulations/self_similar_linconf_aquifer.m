% Time stepping
tDmax = 1;
Nt = 100;
dtD = tDmax/Nt
theta_dt = 0; % theta parameter for the time integration

%% Grid and operators
Grid.xmin = 0; Grid.xmax=10; Grid.Nx=100; 
Grid.geom = 'cartesian';
Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid);

L = -D * G; 
fs = zeros(Grid.N,1);
IM = I + dtD * L * (1-theta_dt);
EX = I - dtD * theta_dt * L;

%% Boundary conditions
BC.dof_dir   = Grid.dof_xmin;
BC.dof_f_dir = Grid.dof_f_xmin;
BC.g         = 0.0;
BC.dof_neu   = [];
BC.dof_f_neu = [];
BC.qb        = [];
[B,N,fn] = build_bnd(BC,Grid,I);

%% Initial condition
hD_num = ones(Grid.N,1);

%% Time stepping
for i=1:Nt
    hD_num = solve_lbvp(IM,EX*hD_num + dtD*(fs+fn),B,BC.g,N); % this is just one line!
end

%% Compute the self-similar solutions
tD      = tDmax;
eta_num = Grid.xc/sqrt(4 * tD);

eta_ana = linspace(0,3,1e2);
hD_ana = erf(eta_ana);

subplot 121
plot(Grid.xc,hD_num,'linewidth',1.5);
xlabel 'x'' ', ylabel 'h'' '
pbaspect([1 .8 1])
subplot 122
plot(eta_ana,hD_ana,'r-'), hold on
plot(eta_num,hD_num,'bo','linewidth',1.5,'markersize',6,'markerfacecolor','w')
xlabel '\eta ', ylabel 'h'' '
xlim([0 3])
pbaspect([1 .8 1])
legend('analytic','numeric')
