% author: Mohammad Afzal Shadab
% date: 26 April 2021
% Description:
% Testing to find the dofs of the faces associated with a given vector of cell dof's
% and separates them into interior and exterior faces. The information 
% which faces belong to which cell is in the divergence matrix.

theta_b = pi-acos(1/3);
Pi_ocean = 0.0;
Pi_hellas = 0.0;


%% Physical parameters
R_Mars = 3389508; % [m] Mars' mean radius 
grav = 3.711;   % [m/s^2] grav. acceleration on Mars

% Hellas geometry
% 42° 24? 0? S, 70° 30? 0? E? (-42.4°, 70.5°)
theta_Hellas  = deg2rad(42+24/60);
phi_Hellas  = deg2rad(70+30/60)
R_Hellas = 1150e3
R_theta_Hellas = R_Hellas/R_Mars;

%% Grid and ops
Grid.xmin = 0; Grid.xmax = theta_b; Grid.Nx = 25;
Grid.ymin = 0; Grid.ymax = 2*pi; Grid.Ny = 75;
Grid = build_grid(Grid);
[Dref,G,I] = build_ops(Grid); % Need a cartesian D for the complex domains!

Grid.geom = 'spherical_shell_theta_phi';
Grid.R_shell = 1.0;
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
[dof_bnd_act,dof_bnd_inact] = find_bnd_cells(dof_act,dof_Hellas,dof_f_Hellas,D,Grid);
N_Hellas = length(dof_bnd_act);


%% Boundary condtions
G(dof_f_Hellas,:) = 0; % Natural boundary conditions at crater
BC.dof_dir   = [Grid.dof_xmax;... % Dichotomy boundary
                dof_bnd_act;...   % Hellas bnd
                dof_Hellas];      % Eliminate inactive cells
BC.dof_f_dir = [Grid.dof_f_xmax;... % Dichotomy boundary
                dof_f_Hellas];   % Hellas bnd
                % Note we don't need the equivalent of the "inactive cells" in this array
BC.g = [Pi_ocean*ones(Grid.Ny,1);...     % Dichotomy boundary
        Pi_hellas*ones(N_Hellas,1);...     % Hellas bnd
       -NaN(size(dof_Hellas))];        % Eliminate inactive cells
[B,N,fn] = build_bnd(BC,Grid,I);

%% Solve linear system
hD = solve_lbvp(L,fs+fn,B,BC.g,N);

% Plotting
[X_f_Hellas,Y_f_Hellas] = comp_face_coords(dof_f_Hellas,Grid);

figure('position',[10 10 500 1000])
subplot 121
plot(rad2deg([Grid.xf';Grid.xf']),rad2deg([Grid.ymin*ones(1,Grid.Nx+1);Grid.ymax*ones(1,Grid.Nx+1)]),'k','linewidth',1), hold on
plot(rad2deg([Grid.xmin*ones(1,Grid.Ny+1);Grid.xmax*ones(1,Grid.Ny+1)]),rad2deg([Grid.yf';Grid.yf']),'k','linewidth',1)
plot(rad2deg(X_f_Hellas),rad2deg(Y_f_Hellas),'-','color','r','linewidth',1.5)
plot(rad2deg(Theta(dof_bnd_act)),rad2deg(Phi(dof_bnd_act)),'go','markersize',4)
plot(Grid.xdom,Grid.ydom,'k-')
axis equal tight
xlabel '\theta [\circ]', ylabel '\phi [\circ]'

subplot 122
contourf(rad2deg(Theta),rad2deg(Phi),reshape(hD,Grid.Ny,Grid.Nx),30)
axis equal tight
xlabel '\theta [\circ]', ylabel '\phi [\circ]'

figure

subplot 121
hellas = ones(Grid.N,1); hellas(dof_Hellas)=0.1;
plot_spherical_shell(hellas,1,Grid)
title 'Hellas planita'
shading interp
view(rad2deg(phi_Hellas+pi/2)-30,rad2deg(theta_Hellas))

subplot 122
plot_spherical_shell(hD,.5,Grid)
title 'Hellas planita'
shading interp
view(rad2deg(phi_Hellas+pi/2)-30,rad2deg(theta_Hellas))


function [] = plot_spherical_shell(hD,scale,Grid)
% Plot sphere
[Theta_s,Phi_s] = meshgrid(linspace(0,pi,100),linspace(0,2*pi,100));
R = 1;
Xs = R*sin(Theta_s).*cos(Phi_s);
Ys = R*sin(Theta_s).*sin(Phi_s);
Zs = R*cos(Theta_s);

% solution
[Theta,Phi] = meshgrid(Grid.xc,Grid.yc);
Hd = reshape(hD,Grid.Ny,Grid.Nx);
Xh = (R+scale*Hd).*sin(Theta).*cos(Phi);
Yh = (R+scale*Hd).*sin(Theta).*sin(Phi);
Zh = (R+scale*Hd).*cos(Theta);

s = surf(Xh,Yh,Zh);
s.CData = Hd;
shading interp
axis equal
end


function [dof_in,dof_out] = find_crater_dofs(theta0,phi0,r,Grid,Theta,Phi)
% We find the cells in the crater by the "great circle distance'
theta = Theta(:); phi = Phi(:);
dof_in =[];
for i = 1:length(r)
    dist_great_circle = 2*asin(sqrt((sin((theta-theta0)/2)).^2 + sin(theta).*sin(theta0).*(sin((phi-phi0)/2)).^2  ));
    dof_in = [dof_in;Grid.dof(dist_great_circle <= r(i))];  
end
dof_in = sort(dof_in,1,'ascend');
dof_out = setdiff(Grid.dof,dof_in);
end
