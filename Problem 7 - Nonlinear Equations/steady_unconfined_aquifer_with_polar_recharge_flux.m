% Steady unconfined aquifer with polar recharge
% Author: Mohammad Afzal Shadab
% Date: 04/11/2021
% Email: mashadab@utexas.edu

%% Parameters
n_exp = 3.0;
Pi = 0.1; 

%% Analytic solution
xD_ana = linspace(0,1,5e2);
hD_ana = @(xD,Pi,n) (Pi.^(n+2) + (n+1)*(n+2)*(1-xD)).^(1/(n+2));


%% Grid and operators
Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 50;
Grid=build_grid(Grid);
[D,G,I] = build_ops(Grid);
fs = zeros(Grid.Nx,1);
M = 0.5*Grid.dx*abs(G);

%% Boundary conditions:
% Dir => for update!
BC.dof_dir = Grid.dof_xmax;
BC.dof_f_dir = Grid.dof_f_xmax;
BC.g = 0;

% Neu => residual
BC.dof_neu = Grid.dof_xmin;
BC.dof_f_neu = Grid.dof_f_xmin;
BC.qb = 1.0;
[B,N,fn] = build_bnd(BC,Grid,I);

%% Residual and Jacobian
% Function in flux term
f  = @(u) u.^(n_exp+1)/(n_exp+1);
df = @(u) u.^n_exp;

F  = @(u) spdiags(M*u.^(n_exp+1)/(n_exp+1),0,Grid.Nx+1,Grid.Nx+1);
dF = @(u) spdiags(u.^n_exp,0,Grid.Nx,Grid.Nx);

% Other 'function matrices'
GU = @(u) spdiags(G*u,0,Grid.Nx+1,Grid.Nx+1);
    
% Residual and Jacobian
res = @(u) D*F(u)*G*u + fs + fn;
Jac = @(u) D*(GU(u)*M*dF(u)+F(u)*G);

%% Newton iteration
tol = 1e-9; % convergence tolerance
kmax = 20;  % maximum number of iterations

% Initial guess
hD = ones(Grid.N,1);
hD(BC.dof_dir) = hD_ana(Grid.xc(Grid.dof_xmax),Pi,n_exp);

%% Newton-Raphson loop
nres = norm(res(hD)); ndhD = 1; k = 0;

while (nres > tol || ndhD > tol) && k < kmax
    dhD = solve_lbvp(Jac(hD),-res(hD),B,BC.g,N);
    hD  = hD + dhD;
    hD(BC.dof_dir) = hD_ana(Grid.xc(Grid.dof_xmax),Pi,n_exp);
    nres = norm(N'*res(hD)); ndhD = norm(N'*dhD);
    k = k+1;
    fprintf('it = %d: nres = %3.2e  ndhD = %3.2e\n',k,nres,ndhD)
    if k == 1; ndhD = 0; end % to allow exit on first iteration
end

%% Plotting
plot(xD_ana,hD_ana(xD_ana,Pi,n_exp),'r-'), hold on
plot(Grid.xc,hD,'b--','markerfacecolor','w','markersize',6)
legend('analytic','numeric')
ylim([0 2])
xlabel 'x'' ', ylabel 'h'''
pbaspect([1 .8 1])
disp 'done'