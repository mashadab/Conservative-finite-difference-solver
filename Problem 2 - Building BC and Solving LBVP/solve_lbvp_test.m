%% Build grid and operators
Grid.xmin = -1; Grid.xmax = 1; Grid.Nx = 100;
Grid = build_grid(Grid);
[D,G,I]=build_ops(Grid);
L = -D*G;                   % Laplacian operator
fs = spalloc(Grid.N,1,0);   % r.h.s. (zero)

%% Setboundary conditions
Param.dof_dir   = [Grid.dof_xmin;Grid.dof_xmax];     % identify cells on Dirichlet bnd
Param.dof_f_dir = [Grid.dof_f_xmin;Grid.dof_f_xmax]; % identify faces on Dirichlet bnd
Param.dof_neu   = [];     % identify cells on Neumann bnd
Param.dof_f_neu = [];     % identify faces on Neumann bnd
Param.g  = [1;2];         % set bnd temperature
Param.qb = [];            % set bnd heat flux
[B,N,fn] = build_bnd(Param,Grid,I);  % Build constraint matrix and basis for its nullspace

%% Solve linear boundary value problem & plot solution
u = solve_lbvp(L,fs+fn,B,Param.g,N);
plot(Grid.xc,u)
xlabel 'x', ylabel 'u'
ylim([0 2])