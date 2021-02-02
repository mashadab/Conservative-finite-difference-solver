function [D,G,I]=build_ops(Grid)
% author: Mohammad Afzal Shadab
% date: 02/02/2021
% description:
% This function computes the discrete divergence and gradient matrices on a
% regular staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions.
% Input:
% Grid = structure containing all pertinent information about the grid.
% Output:
% D = Nx by Nx+1 discrete divergence matrix 
% G = Nx+1 by Nx discrete gradient matrix
% I = Nx by Nx identity matrix
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);

% Making a vector of all ones of length Nx.1
e = ones(Grid.Nx,1);

% 1) Build sparse Divergence operator
D = spdiags([-e e]/Grid.dx,0:1,Grid.Nx,Grid.Nfx);
% 2) Obtain sparse Gradient operator in interior from D
G = -D';
% 3) Set natural (homogeneous Neumann) boundary conditions
dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax] ; % all dof's on boundary
G(dof_f_bnd,:) = 0;

% Identity
I = speye(Grid.Nx);