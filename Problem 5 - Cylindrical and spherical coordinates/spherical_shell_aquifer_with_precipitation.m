%% Problem 3: Aquifer on spherical cap with precipitation

theta_p = deg2rad(5);
theta_b = pi-acos(1/3);
theta_ana = linspace(theta_p,theta_b,100);

hD_ana = @(theta,theta_p) log(sin(theta)/sin(theta_b)) + cos(theta_p).*log((csc(theta)+cot(theta))./(csc(theta_b)+cot(theta_b)));
qD_ana = @(theta,theta_p) cos(theta_p)*csc(theta)-cot(theta);

%% Grid and ops
Grid.xmin = theta_p; 
Grid.xmax = theta_b; 
Grid.Nx = 20;
Grid.geom = 'spherical_shell';
Grid.R_shell = 1.0;
Grid = build_grid(Grid);

% Operators
[D,G,I] = build_ops(Grid);
L = - D * G;
fs = ones(Grid.N,1);

%% Boundary conditions
BC.dof_dir = Grid.dof_xmax;
BC.dof_f_dir = Grid.dof_f_xmax; %At theta_b, we have zero head
BC.g = hD_ana(Grid.xc(BC.dof_dir),theta_p)
BC.dof_f_neu = []; %At theta_b, we have zero head
BC.dof_neu = [];
BC.qb = [];
[B,N,fn] = build_bnd(BC,Grid,I);

%% Solve for hD and qD
hD = solve_lbvp(L,fs+fn,B,BC.g,N);
qD = comp_flux(D,1.0,G,hD,fs,Grid,BC);

subplot(2,2,1)
plot(rad2deg(theta_ana),hD_ana(theta_ana,theta_p),'-','linewidth',1.5), hold on
plot(rad2deg(Grid.xc),hD,'o','MarkerFaceColor','w')
xlabel('\theta'), ylabel('h''')%, pbaspect([1 .2 1])
xlim([0 120])
legend('analytical','numerical','location','southwest')

subplot(2,2,3)
plot(rad2deg(theta_ana),qD_ana(theta_ana,theta_p),'-','linewidth',1.5), hold on
plot(rad2deg(Grid.xf),qD,'o','MarkerFaceColor','w')
xlabel('\theta'), ylabel('q''')% pbaspect([1 .2 1])
xlim([0 120])
legend('analytical','numerical','location','northwest')

subplot(2,2,[2 4])
h_scale_plot = 0.3;
x_h = (1+hD*h_scale_plot).*sin(Grid.xc);
z_h = (1+hD*h_scale_plot).*cos(Grid.xc);
plot(x_h,-z_h,'-','linewidth',1.5), hold on
theta_sphere = linspace(0,pi,5e2);
x_base = sin(theta_sphere);
z_base = cos(theta_sphere);
plot(x_base,z_base,'k','linewidth',1.5)
axis equal