function [D,G,I]=build_ops(Grid)
% author: Mohammad Afzal Shadab
% date: 04/19/2021
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
    Dx = spdiags([-ones(Grid.Nx,1) ones(Grid.Nx,1)]/Grid.dx,[0 1],Grid.Nx,Grid.Nx+1); % 1D div-matrix in x-dir
    Dy = spdiags([-ones(Grid.Ny,1) ones(Grid.Ny,1)]/Grid.dy,[0 1],Grid.Ny,Grid.Ny+1); % 1D div-matrix in y-dir
    Ix = eye(Grid.Ny); % 1D Nx-identity
    Iy = eye(Grid.Nx); % 1D Ny-identity

    % 2D Tensor-product divergence matrices
    Dx = kron(Dx,Ix);  % 2D div-matrix in x-dir
    Dy = kron(Iy,Dy);  % 2D div-matrix in y-dir

    % Complete 2D divergence
    D = [Dx,Dy];
    % Boundary faces
    dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax...
                ;Grid.dof_f_ymin; Grid.dof_f_ymax];
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
