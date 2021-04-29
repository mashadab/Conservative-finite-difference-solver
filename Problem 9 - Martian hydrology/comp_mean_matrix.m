function [Kd] = comp_mean_matrix(K,p,Grid) % repo MDOT 
% author: Mohammad Afzal Shadab	
% date: 9 April 2021
% Description:
% Takes coefficient field, K, defined at the cell centers and computes the
% mean specified by the power, p and returns it in a sparse diagonal
% matrix, Kd.
%
% Input:
% K = N by 1 column vector of cell centered values
% p = power of the generalized mean
%       1 (arithmetic mean)
%      -1 (harmonic mean)
% Grid = structure containing information about the grid.
%
% Output:
% Kd = Nf by Nf diagonal matrix of power means at the cell faces.
%
% Example call:
% K = @(x) 1+x.^3;
% Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% Grid = build_grid(Grid);
% Kd = comp_mean(K(Grid.xc),1,Grid);
Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; N = Grid.N;

if (Nx>1) && (Ny>1)  % 2D case
    Mx = spdiags([ones(Nx,1) ones(Nx,1)]/2,[-1 0],Nx+1,Nx); % 1D mean-matrix in x-dir
    My = spdiags([ones(Ny,1) ones(Ny,1)]/2,[-1 0],Ny+1,Ny); % 1D mean-matrix in y-dir
    Mx = kron(Mx,speye(Ny));                                 % 2D mean-matrix in x-dir
    My = kron(speye(Nx),My);                                 % 2D mean-matrix in y-dir
    M  = [Mx;My];                                            % 2D mean-matrix
    dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax;...
                 Grid.dof_f_ymin; Grid.dof_f_ymax];          % 0 on bnd's
    M(dof_f_bnd) = 0;
    
elseif (Nx>1) && (Ny==1) % 1D x-direction
    M = spdiags([ones(Nx,1) ones(Nx,1)]/2,[-1 0],Nx+1,Nx); % 1D mean-matrix in x-dir
    dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax];   % boundary faces
elseif (Nx==1) && (Ny>1) % 1D y-direction
    M = spdiags([ones(Ny,1) ones(Ny,1)]/2,[-1 0],Ny+1,Ny); % 1D mean-matrix in y-dir
    dof_f_bnd = [Grid.dof_f_ymin; Grid.dof_f_ymax]-Grid.Nfx;   % boundary faces
end

% for spherical shell coordinates we need periodicity in the y-direction
if strcmp(Grid.geom,'spherical_shell_theta_phi')
    Mx = spdiags([ones(Nx,1) ones(Nx,1)]/2,[-1 0],Nx+1,Nx); % 1D mean-matrix in x-dir
    My = spdiags([ones(Ny,1) ones(Ny,1)]/2,[-1 0],Ny+1,Ny); % 1D mean-matrix in y-dir
    My(1,Ny) = 1/2; My(Ny+1,1) = 1/2;        % periodicity in y-dir
    Mx = kron(Mx,speye(Ny));                 % 2D mean-matrix in x-dir
    My = kron(speye(Nx),My);                 % 2D mean-matrix in y-dir
    M  = [Mx;My];                            % 2D mean-matrix
    dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax];   % 0 on bnd's
    M(dof_f_bnd,:) = 0;   
end

Kmean = (M*K.^p).^(1/p);                % Compute general power mean
Kmean(dof_f_bnd) = 0;                   % apply zero bnd's
Kd = spdiags(Kmean,0,Grid.Nf,Grid.Nf);  % place on diagonal
