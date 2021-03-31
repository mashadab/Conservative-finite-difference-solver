% Self-similar recharge in a linear confined aquifer
% Author: Mohammad Afzal Shadab
% Date: 03/31/2021
% Email: mashadab@utexas.edu

% Time stepping
tDmax = 1;
Nt = 100;
dtD = tDmax/Nt;
theta_dt = 0;

%% Grid and operators
Grid.xmin = 0; Grid.xmax = 10; Grid.Nx = 1e2;
Grid.geom = 'cartesian';
Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid);

L = - D * G;
fs = zeros(Grid.N,1);
IM = I + dtD * L * (1-theta_dt);
EX = I - dtD * theta_dt * L;

%% Boundary conditions
BC.dof_dir   = Grid.dof_xmax;
BC.dof_f_dir = Grid.dof_f_xmax;
BC.g         = 0.0;
BC.dof_neu   = Grid.dof_xmin;
BC.dof_f_neu = Grid.dof_f_xmin;
BC.qb        = 1.0;
[B,N,fn] = build_bnd(BC,Grid,I);
%% Initial condition
hD_num = zeros(Grid.N,1)

%% Time stepping
for i=1:Nt
    hD_num = solve_lbvp(IM,EX*hD_num + dtD*(fs+fn),B,BC.g,N);
end

%% Self-similar solution
tD = tDmax;
eta = Grid.xc / sqrt(4*tD);
Pi = hD_num / sqrt(4*tD);

eta_ana = linspace(0,2,1e2);
Pi_ana = @(eta) exp(-eta.^2)/sqrt(pi) - eta.*erfc(eta);

%% Plotting
subplot 121
plot(Grid.xc,hD_num,'linewidth',1.5);
xlabel 'x'' ', ylabel 'h'' '
pbaspect([1 .8 1])
xlim([0 5]), ylim([0 1.2])
set(gca,'ytick',[0:.4:1.2])

subplot 122
plot(eta_ana,Pi_ana(eta_ana),'r-'), hold on
plot(eta,Pi,'bo','linewidth',1.5,'markersize',6,'markerfacecolor','w')
xlabel '\eta ', ylabel '\Pi '
xlim([0 2]), ylim([0 .8])
pbaspect([1 .8 1])
set(gca,'ytick',[0:.2:.8])
legend('analytic','numeric')
