%% Aquifer with a dike
% Parameters
Ka = 1e-11;     % [m^2] aquifer permeability 
Kd = 1e-14;     % [m^2] dike permeability

% Grid and ops
Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 100;
Grid = build_grid(Grid)
[D,G,I] = build_ops(Grid)

% Permeability
K_hom = ones(Grid.N,1);
K_het = K_hom; K_het(51) = Kd/Ka;
Kd_hom = comp_mean_matrix(K_hom,-1,Grid,G);
Kd_het = comp_mean_matrix(K_het,-1,Grid,G);
L_hom = -D*Kd_hom*G;    % heterogeneus operator
L_het = -D*Kd_het*G;    % homogeneous operator
fs = ones(Grid.N,1);

% BC's
BC.dof_dir   = Grid.dof_xmax;
BC.dof_f_dir = Grid.dof_f_xmax;
BC.g = 0.0;

BC.dof_neu = [];
BC.dof_f_neu = [];
BC.qb = [];
[B,N,fn] = build_bnd(BC,Grid,I);

% Solve systems
hD_hom = solve_lbvp(L_hom,fs+fn,B,BC.g,N);
qD_hom = comp_flux(D,Kd_hom,G,hD_hom,fs,Grid,BC);

hD_het = solve_lbvp(L_het,fs+fn,B,BC.g,N);
qD_het = comp_flux(D,Kd_het,G,hD_het,fs,Grid,BC);

% 
subplot 121
plot(Grid.xc,hD_hom,'-',Grid.xc,hD_het,'--')
xlabel 'x'' ', ylabel 'h'' ', pbaspect([1 .8 1])
legend('homogeneous','heterogeneous','location','northeast')

subplot 122
plot(Grid.xf,qD_hom,'-',Grid.xf,qD_het,'--')
xlabel 'x'' ', ylabel 'q'' ', pbaspect([1 .8 1])
legend('homogeneous','heterogeneous','location','northwest')