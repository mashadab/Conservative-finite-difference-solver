function [D,G,I]=build_ops(Grid)
% author: Mohammad Afzal Shadab
% date: -3/17/2021
% description:
% This function computes the discrete divergence and gradient matrices on a
% regular staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous Neumann boundary conditions.
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

Nx = Grid.Nx; 
% 1) Build sparse Divergence operator
D = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[0 1],Nx,Nx+1);
% 2) Obtain sparse Gradient operator in interior
G = -D';
% 3) Set natural (homogeneous Neumann) boundary conditions
dof_f_bnd = [Grid.dof_f_xmin,Grid.dof_f_xmax]; % all dof's on boundary
G(dof_f_bnd,:) = 0;

% Identity
I = speye(Grid.Nx);

if strcmp(Grid.geom,'cylindrical_r')
    fprintf('Operators built for 1D cylindrical geometry\n.')
    Rf = spdiags(Grid.xf,0,Grid.Nx+1,Grid.Nx+1);
    Rcinv = spdiags(1./Grid.xc,0,Grid.Nx,Grid.Nx);
    D = Rcinv*D*Rf;
elseif strcmp(Grid.geom,'spherical_shell')
    fprintf('Operators built for 1D spherical shell geometry (R = %3.2e)n.',Grid.R_shell)
    SINf = spdiags(sin(Grid.xf),0,Grid.Nx+1,Grid.Nx+1);
    SINcinv = spdiags(1./sin(Grid.xc),0,Grid.Nx,Grid.Nx);
    D = 1/Grid.R_shell*SINcinv*D*SINf;
    G = G/Grid.R_shell;
elseif strcmp(Grid.geom,'cartesian')
    % nothing needs to be done
else
    error('Unknown geometry.')
end
