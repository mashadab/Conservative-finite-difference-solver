% Grid and ops
Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 30;
Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid);
L = -D*G;
fs_con = ones(Grid.N,1); % constant source term
fs_var = zeros(Grid.N,1); fs_var(Grid.xc>=0.5) = 2+2*cos(2*pi*Grid.xc(Grid.xc>=0.5)); % variable source term

% BC's
BC.dof_dir   = Grid.dof_xmax;
BC.dof_f_dir = Grid.dof_f_xmax;
BC.g = 0;

BC.dof_neu = [];
BC.dof_f_neu = [];
BC.qb = [];
[B,N,fn] = build_bnd(BC,Grid,I);

% Solve systems
hD_con = solve_lbvp(L,fs_con+fn,B,BC.g,N); % soln with constant source term
qD_con = comp_flux(D,1.0,G,hD_con,fs_con,Grid,BC); % flux with constant source term

hD_var = solve_lbvp(L,fs_var+fn,B,BC.g,N); % soln with variable source term
qD_var = comp_flux(D,1.0,G,hD_var,fs_var,Grid,BC); % flux with variable source term

%% Check mass concervation
err_con = sum(Grid.dx*fs_con)-qD_con(Grid.dof_f_xmax)
err_var = sum(Grid.dx*fs_var)-qD_var(Grid.dof_f_xmax)

subplot 311
plot(Grid.xc,fs_con,'-',Grid.xc,fs_var,'-')
xlabel 'x'' ', ylabel 'q_p'' ', pbaspect([1 .2 1])

subplot 312
plot(Grid.xf,qD_con,'-',Grid.xf,qD_var,'-')
xlabel 'x'' ', ylabel 'q'' ', pbaspect([1 .2 1])

subplot 313
plot(Grid.xc,hD_con,'-',Grid.xc,hD_var,'-')
xlabel 'x'' ', ylabel 'h'' ', pbaspect([1 .2 1])