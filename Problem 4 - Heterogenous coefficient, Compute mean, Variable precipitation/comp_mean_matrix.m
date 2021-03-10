function [Kd] = comp_mean_matrix(K,p,Grid,G)
% author: Mohammad Afzal Shadab
% date: 10 March 2021
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
% [D,G,I]=build_ops(Grid);
% Kd = comp_mean_matrix(K(Grid.xc),1,Grid,G);

if (p == -1) | (p == 1)
    M = 0.5*Grid.dx*abs(G); % averaging matrix
    Kmean = (M*K.^p).^(1/p); % generalized mean vector Nf by 1
    dof_f_bnd = [Grid.dof_f_xmin;Grid.dof_f_xmax]; % all boundary dof's
    Kmean(dof_f_bnd) = 0; % Zero out bnd's
    Kd = spdiags([Kmean],0,Grid.Nfx,Grid.Nfx); % Diagonal matrix
else
    error('This power does not have significance.')
end