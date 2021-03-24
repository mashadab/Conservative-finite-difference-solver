%% Parameters
yr2s = 60^2*24*365.25;  % second per year
R_mars = 3389508;    % [m] Mars' mean radius
grav = 3.711;           % [m/s^2] grav. acceleration on Mars
rho = 1e3;              % [kg/m^3] desity of water 
k = 1e-11;              % [m^2] permeability (Hanna & Phillips 2005)
mu = 1e-3;              % [Pa s] water viscosity
h_b = -500;             % [m] sea level
b = 5e3;                % [m] aquifer thickness
Qi = 1e9/yr2s           % [m^3/s] polar rechrge rate Q = 1km^3/yr from Clifford 1993

% derived values
K = k * rho * grav / mu;      % [m/s] hydraulic conductivity

theta_p = deg2rad(1);
theta_b = pi-acos(1/3);
theta_ana = linspace(theta_p,theta_b,100);

%% Characteristic scales
h_c = @(theta_p) Qi./(2*pi*b*K*sin(theta_p));
q_c = @(theta_p) K.*h_c(theta_p)/R_mars;

%% Analytic solutions
% Dimensionless
hD_ana = @(theta,theta_p) sin(theta_p).*log((csc(theta)+cot(theta))./(csc(theta_b)+cot(theta_b)));
qD_ana = @(theta,theta_p) sin(theta_p)./sin(theta);
% Dimensional
h_ana = @(theta,theta_p) h_b + h_c(theta_p)*hD_ana(theta,theta_p);
q_ana = @(theta,theta_p) q_c(theta_p).*qD_ana(theta,theta_p);


%% Grid and ops
Grid.xmin = theta_p;
Grid.xmax = theta_b;
Grid.Nx   = 20;
Grid.geom = 'spherical_shell';
Grid.R_shell = 1.0;
Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid);
L = - D * G;
fs = zeros(Grid.N,1);

%% BC's
BC.dof_dir = Grid.dof_xmax;
BC.dof_f_dir = Grid.dof_f_xmax;
BC.g = hD_ana(Grid.xc(BC.dof_dir),theta_p);
BC.dof_neu = Grid.dof_xmin;
BC.dof_f_neu = Grid.dof_f_xmin;
BC.qb = 1.0;
[B,N,fn] = build_bnd(BC,Grid,I);

%% Compute hD and qD
% Dimensionless solution
hD = solve_lbvp(L,fs+fn,B,BC.g,N);
qD = comp_flux(D,1.0,G,hD,fs,Grid,BC);

% Dimensional solution
h = h_b + h_c(theta_p)*hD;
q = q_c(theta_p)*qD;

%% Plotting
subplot(2,2,1)
plot(rad2deg(theta_ana),h_ana(theta_ana,theta_p),'-','linewidth',1.5), hold on
plot(rad2deg(Grid.xc),h,'o','MarkerFaceColor','w')
xlabel('\theta [\circ]'), ylabel('h [m]')
xlim([0 120])
legend('analytical','numerical','location','northeast')

subplot(2,2,3)
plot(rad2deg(theta_ana),qD_ana(theta_ana,theta_p),'-','linewidth',1.5), hold on
plot(rad2deg(Grid.xf),qD,'o','MarkerFaceColor','w')
xlabel('\theta [\circ]'), ylabel('q [m/s]')
xlim([0 120])
legend('analytical','numerical','location','northeast')

subplot(2,2,[2 4])
theta_sphere = linspace(0,pi,5e2);
x_base = sin(theta_sphere);
z_base = cos(theta_sphere);
plot(x_base,z_base,'k','linewidth',1.5), hold on
h_scale_plot = 2;
x_h = (1+hD_ana(theta_ana,theta_p)*h_scale_plot).*sin(theta_ana);
z_h = (1+hD_ana(theta_ana,theta_p)*h_scale_plot).*cos(theta_ana);
plot(x_h,-z_h,'r-','linewidth',1.5)
axis equal
xlim([0 1.5]), ylim(1.5*[-1 1])
