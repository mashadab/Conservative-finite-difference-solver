function [dof_f_ext, dof_f] = find_faces(cell_dofs,D,Grid)
% author: Mohammad Afzal Shadab
% date: 26 April 2021
% Description:
% Find the dofs of the faces associated with a given vector of cell dof's
% and separates them into interior and exterior faces. The information 
% which faces belong to which cell is in the divergence matrix.
% 
% Input:
% cell_dofs = column vector of cells in the active subdomain
% D         = discrete divergence matrix
% Grid      = structure containing pertinent information about the grid
%
% Output:
% dof_f    = column vector containing the dof's of all faces in and on 
%             the exterior of the active subdomain

% Sample call:
% Grid.xmin = -2; Grid.xmax = 2; Grid.Nx = 50;
% Grid.ymin = -2; Grid.ymax = 2; Grid.Ny = 50;
% theta = linspace(0,2*pi,4);
% x_poly = cos(theta);
% y_poly = sin(theta);
% [Xc,Yc] = meshgrid(Grid.xc,Grid.yc);
% [dof_in,dof_out] = find_poly_cells(x_poly,y_poly,Grid,Xc,Yc);
% [dof_f_bnd,dof_f] = find_faces(dof_in,D,Grid);


%% Select appropriate rows of the divergence matrix
DD = D(cell_dofs,:);

%% Exterior faces
% Sum columns of DD. In the columns corresponding to internal faces the
% entries will cancel. Therefore, only the sum of the columns corresponding
% to exterior faces will be non-zero.
dof_f_ext = Grid.dof_f(abs(sum(DD,1))>1e-10);
dof_f     = Grid.dof_f(abs(sum(abs(DD),1))>1e-10);