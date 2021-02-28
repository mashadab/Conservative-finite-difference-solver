% file: SteadyLinearConfinedAquifer.m
% date: 1894
% author: Percival Lowel

% Description: Solve for the steady linear confined aquifer in the southern
% highlands of Mars, without the symmetry assumption at the south pole.
% See problem statement in Matlab Grader for details.

%% Problem parameters
% primary
Param.yr2s = 60^2*24*365.25;    % [-] seconds per year
Param.theta_bnd = acos(1/3);    % [rad] co-lattitude of dichotomy b
Param.R = 3389508;              % [m] Mars' mean radius 
Param.grav = 3.711;             % [m/s^2] grav. acceleration on Mars
Param.rho = 1e3;                % [kg/m^3] density of water 
Param.grav = 3.711;             % [m/s^2] grav. acceleration on Mars
Param.k = 1e-11;                % [m^2] permeability (Hanna & Phillips 2005)
Param.mu = 1e-3;                % [Pa s] water viscosity
Param.b = 5e3;                  % [m] aquifer thickness
Param.fp = 0.13;                % [mm/yr] mean precipitation
Param.ho = -500;                % [m] sealevel (for now)

% secondary
Param.l = Param.R * (pi - Param.theta_bnd);  % [m] distance from south pole to dichotomy boundary
Param.K = Param.k * Param.rho * Param.grav / Param.mu;  % [m/s] hydraulic conductivity
Param.fs = Param.fp * 1e-3 / Param.yr2s;  % [m/s] mean precipitation

% Characterisitc scales
Scales.x_c = Param.l;
Scales.h_c = Param.fs * Param.l^2 / (Param.b * Param.K);

%% Generate grid and operators
Grid.xmin = -1;
Grid.xmax = 1;
Grid.Nx = 100;
Grid = build_grid(Grid);

[D,G,I] = build_ops(Grid);
L = -D*G;
fs = 1.0 * ones(Grid.N,1);

%% Set boundary conditions
BC.dof_dir = [Grid.dof_xmin;Grid.dof_xmax];
BC.g = [0.0;0.0]; 

[B,N,fn] = build_bnd(BC,Grid,I);

%% Solve problem
% dimensionless solution
h = solve_lbvp(L,fs+fn,B,BC.g,N);

% dimensional solution
x_dim = Scales.x_c * Grid.xc;
h_dim = Param.ho + Scales.h_c * h;


%% Plot results
subplot 121
plot(Grid.xc,h)
xlabel 'x''', ylabel 'h''', title 'dimensionless'
pbaspect([1 .8 1])

subplot 122
plot(x_dim/1e3,h_dim)
xlabel 'x [km]', ylabel 'h [m]', title 'dimensional'
pbaspect([1 .8 1])