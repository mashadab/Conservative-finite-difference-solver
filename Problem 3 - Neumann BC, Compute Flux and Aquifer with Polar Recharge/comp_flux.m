function [q] = comp_flux(D,Kd,G,u,fs,Grid,BC)
% author: Mohammad Afzal Shadab
% date:  03/03/2021
% Description:
% Computes the conservative fluxes across all boundaries from the 
% residual of the compatability condition over the boundary cells.
% Note: Current implmentation works for all cases where one face 
%       is assigend to each bnd cell. So corner cells must have
%       natual BC's on all but one face.
%
% Input:
% D = Nx by Nfx discrete divergence matrix.
% Kd = Nfx by Nfx conductivity matrix.
% G = Nfx by Nx discrete gradient matrix.
% u = Nx by 1 vector of the potential in cell centers.
% fs = Nx by 1 right hand side vector containing only source terms.
% Grid = structure containing grid information.
% BC = structure contaning information about BC's
%
% Output:
% q = Nfx by 1 vector of fluxes across all cell faces,
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I] = build_ops(Grid);
% >> L = -D*G; fs = ones(Grid.Nx,1);
% >> BC.dof_dir = Grid.dof_xmin;
% >> BC.dof_f_dir = Grid.dof_f_xmin;
% >> BC.g = 0;
% >> BC.dof_neu = []; BC.dof_f_neu =[];
% >> [B,N,fn] = build_bnd(BC,Grid);
% >> h = solve_lbvp(L,fs+fn,B,BC.g,N);
% >> q = comp_flux(D,1,G,h,fs,Grid,BC);

%% Compute interior fluxes
q = - Kd*G*u
    
%% Compute residual everywhere
res   = D*q - fs;        % column vector of residual in all cells

%% Reconstruct boundary fluxes
dof_cells = [BC.dof_dir;BC.dof_neu]; % column vector of all cells on the boundary with non-natural BC's
dof_faces = [BC.dof_f_dir;BC.dof_f_neu]; % column vector of all faces on the boundary with non-natural BC's
res_bnd   = res(dof_cells); % column vector of residual in all dof_bnd cells
the_sign  = ismember(dof_faces,Grid.dof_f_xmin)-ismember(dof_faces,Grid.dof_f_xmax) ; % column vector of flips the sign on xmax, ymax, zmax (in 3D) bnd's

q(dof_faces) =  the_sign.*res_bnd.*Grid.V(dof_cells)./Grid.A(dof_faces);