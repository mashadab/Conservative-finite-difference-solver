function [D,G,I] = build_ops(Grid)
% author: Marc Hesse
% date: 09/08/2014, 09/23/2016, 12/31/2017
% description:
% This function computes the discrete divergence and gradient matrices on a
% regular staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions.
% Input:
% Grid = structure containing all pertinent information about the grid.
% Output:
% D = discrete divergence matrix
% G = discrete gradient matrix
% I = identity matrix
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 4;
% >> Grid.ymin = 0; Grid.ymax = 1; Grid.Ny = 3;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);

% this will help
Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; N = Grid.N;

if (Nx>1) && (Ny>1)  % 2D case
    % 1D divergence matrices
    Dx = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[0 1],Nx,Nx+1); % 1D div-matrix in x-dir
    Dy = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[0 1],Ny,Ny+1); % 1D div-matrix in y-dir
    Ix = speye(Nx); Iy = speye(Ny);  % 1D identities in x and y dirs

    % 2D Tensor-product divergence matrices
    Dx = kron(Dx,Iy);  % 2D div-matrix in x-dir
    Dy = kron(Ix,Dy);  % 2D div-matrix in y-dir

    % Complete 2D divergence
    D = [Dx Dy];
    dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax;... % boundary faces
                 Grid.dof_f_ymin; Grid.dof_f_ymax];
elseif (Nx>1) && (Ny==1)
    D = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[0 1],Nx,Nx+1);
    dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax];   % boundary faces
elseif (Nx==1) && (Ny>1)
    D = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[0 1],Ny,Ny+1);
    dof_f_bnd = [Grid.dof_f_ymin; Grid.dof_f_ymax];   % boundary faces
end
%% Gradient
% adjoint relation
G = -D';
% natural bc's
G(dof_f_bnd,:) = 0;

%% Identity
I = speye(Grid.N);

%% Other geometries
if strcmp(Grid.geom,'cylindrical_r')
    fprintf('Operators built for 1D cylindrical geometry\n.')
    Rf = spdiags(Grid.xf,0,Nx+1,Nx+1);
    Rcinv = spdiags(1./Grid.xc,0,Nx,Nx);
    D = Rcinv*D*Rf;
elseif strcmp(Grid.geom,'spherical_shell')
    fprintf('Operators built for 1D spherical shell geometry (R = %3.2e)n.',Grid.R_shell)
    Sin_f = spdiags(sin(Grid.xf),0,Nx+1,Nx+1);
    R_sin_c_inv = spdiags(1./(Grid.R_shell*sin(Grid.xc)),0,Nx,Nx);
    D = R_sin_c_inv*D*Sin_f;
    G = G/Grid.R_shell;
elseif strcmp(Grid.geom,'spherical_shell_theta_phi')
    % assumes x is polar angle and y is the azimutal angle
    Dx = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[0 1],Nx,Nx+1); % 1D div-matrix in x-dir
    Gx = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[-1 0],Nx+1,Nx); % 1D grad-matrix in x-dir
    Dy = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[0 1],Ny,Ny+1);  % 1D div-matrix in y-dir
    Gy = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[-1 0],Ny+1,Ny); % 1D grad-matrix in y-dir
    Gy(1,Ny) = -1/Grid.dy; Gy(Ny+1,1)=1/Grid.dy;   % periodic BC's in y-dir
    
    % polar angle
    Sin_f = spdiags(sin(Grid.xf),0,Nx+1,Nx+1);
    R_sin_c_inv = spdiags(1./(Grid.R_shell.*sin(Grid.xc)),0,Nx,Nx);
    Dx = R_sin_c_inv*Dx*Sin_f;  % 1D polar D
    Gx = Gx./Grid.R_shell;      % 1D polar G
    Dx = kron(Dx,speye(Ny));    % 2D polar D
    Gx = kron(Gx,speye(Ny));    % 2D polar G
    
    % azimuthal angle
    % Note: The the R_sin_c_inv matrix is at cell centers for both Dy and Gy
    %       because it is the cell center in the x-dir!
    R_sin_c_inv = spdiags(1./(Grid.R_shell.*sin(Grid.xc)),0,Nx,Nx);
    Dy = kron(R_sin_c_inv,Dy);    % 2D azimutal D
    Gy = kron(R_sin_c_inv,Gy);    % 2D azimutal G
    D  = [Dx,Dy];
    G  = [Gx;Gy];
    
    dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax]; % Natural BC's x-dir
    G(dof_f_bnd,:) = 0;
elseif strcmp(Grid.geom,'cartesian')
    % nothing needs to be done
else
    error('Unknown geometry.')
end
