% author: Mohammad Afzal Shadab
% date: 26 April 2021
% Description:
% Testing to find the dofs of the faces associated with a given vector of cell dof's
% and separates them into interior and exterior faces. The information 
% which faces belong to which cell is in the divergence matrix.

theta_b = pi-acos(1/3);

%% Physical paramters
R_Mars = 3389508; % [m] Mars' mean radius 
grav = 3.711;   % [m/s^2] grav. acceleration on Mars

% Hellas geometry
% 42° 24? 0? S, 70° 30? 0? E? (-42.4°, 70.5°)
theta_Hellas  = deg2rad(41+24/60);
phi_Hellas  = deg2rad(73+30/60);
R_Hellas = 2200e3/2;
R_theta_Hellas = R_Hellas/R_Mars;

%% Grid and ops
Grid.xmin = 0; Grid.xmax = theta_b; Grid.Nx = 25;
Grid.ymin = 0; Grid.ymax = 2*pi; Grid.Ny = 75;
Grid = build_grid(Grid);
[Dref,G,I] = build_ops(Grid); % Need a cartesian D for the complex domains!

Grid.geom = 'spherical_shell_theta_phi';
Grid.R_shell = 1;
Grid = build_grid(Grid);
[Theta,Phi] = meshgrid(Grid.xc,Grid.yc);

% Operators
[D,G,I] = build_ops(Grid);
L = -D*G;
fs = ones(Grid.N,1);

%% Hellas geometry
% 1. find Hellas dof's
[dof_Hellas,dof_act] = find_crater_dofs(theta_Hellas,phi_Hellas,R_theta_Hellas,Grid,Theta,Phi);
dof_f_Hellas = find_faces(dof_Hellas,Dref,Grid);


% Plotting
[X_f_Hellas,Y_f_Hellas] = comp_face_coords(dof_f_Hellas,Grid);

plot(rad2deg([Grid.xf';Grid.xf']),rad2deg([Grid.ymin*ones(1,Grid.Nx+1);Grid.ymax*ones(1,Grid.Nx+1)]),'k','linewidth',1), hold on
plot(rad2deg([Grid.xmin*ones(1,Grid.Ny+1);Grid.xmax*ones(1,Grid.Ny+1)]),rad2deg([Grid.yf';Grid.yf']),'k','linewidth',1)
plot(rad2deg(X_f_Hellas),rad2deg(Y_f_Hellas),'-','color','r','linewidth',2)
plot(Grid.xdom,Grid.ydom,'k-')
axis equal tight
xlabel '\theta [\circ]', ylabel '\phi [\circ]'

