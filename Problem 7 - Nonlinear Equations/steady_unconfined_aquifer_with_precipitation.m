% Steady unconfined aquifer with precipitation
% Author: Mohammad Afzal Shadab
% Date: 04/11/2021
% Email: mashadab@utexas.edu

%% Parameters
K0 = 1;
qp = 1;
l = 1;
hb = .3;
n_exp = 5.0;
hc = @(n) (qp*l^2/K0)^(1/(n+2));
Pi = @(n) hb/(qp*l^2/K0)^(1/(n+2));

%% Analytic solution
xD_ana = linspace(0,1,5e2);
hD_ana = @(x,n) (Pi(n).^(n+2) + 0.5*(n+1)*(n+2)*(1-x.^2)).^(1/(n+2));

%% Grid and operators
Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 100;
Grid=build_grid(Grid);
[D,G,I] = build_ops(Grid);
fs = ones(Grid.Nx,1);
M = 0.5*Grid.dx*abs(G);

%% Residual and Jacobian
% Function in flux term
f  = @(u) u.^(n_exp+1)/(n_exp+1);
df = @(u) u.^n_exp;

F  = @(u) spdiags(M*u.^(n_exp+1)/(n_exp+1),0,Grid.Nx+1,Grid.Nx+1);
dF = @(u) spdiags(u.^n_exp,0,Grid.Nx,Grid.Nx);

% Other 'function matrices'
GU = @(u) spdiags(G*u,0,Grid.Nx+1,Grid.Nx+1);

% Residual and Jacobian
res = @(u) D*F(u)*G*u + fs;
Jac = @(u) D*(GU(u)*M*dF(u)+F(u)*G);

%% Boundary conditions - for update!
BC.dof_dir = Grid.dof_xmax; 
BC.dof_f_dir = Grid.dof_f_xmax;
BC.g= hD_ana(Grid.xc(Grid.dof_xmax),n_exp);
BC.dof_neu= [];
BC.dof_f_neu= [];
BC.qb = [];
[B,N,fn] = build_bnd(BC,Grid,I);

%% Newton iteration
tol = 1e-9; % convergence tolerance
kmax = 20;  % maximum number of iterations

% Initial guess
hD = ones(Grid.N,1);

nres = norm(res(hD)); ndhD = 1; k = 0;
while (nres > tol || ndhD > tol) && k < kmax
    dhD = solve_lbvp(Jac(hD),-res(hD),B,BC.g,N);
    hD  = hD + dhD;
    
    nres = norm(N'*res(hD)); ndhD = norm(N'*dhD);
    k = k+1;
    fprintf('it = %d: nres = %3.2e  ndhD = %3.2e\n',k,nres,ndhD)
    if k == 1; ndhD = 0; end % to allow exit on first iteration
end

plot(xD_ana,hD_ana(xD_ana,n_exp),'r-'), hold on
plot(Grid.xc,hD,'b--','markerfacecolor','w','markersize',6)
legend('analytic','numeric')
ylim([0 2])
xlabel 'x'' ', ylabel 'h'''
pbaspect([1 .8 1])
