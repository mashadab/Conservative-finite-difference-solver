function [B,N,fn] = build_bnd(Param,Grid,I)
% author: Mohammad Afzal Shadab
% date: 02/10/2021
% Description:
% This function computes the operators for Dirichlet boundary conditions. 
%
% Input:
% Param = structure containing all information about the physical problem
%         in particular this function needs the fields
%         Param.dof_dir = Nc by 1 column vector containing 
%                         the dof's of the cells on the Dirichlet boundary.
%         Param.dof_f_dir = Nc by 1 column vector containing on the Dirichlet boundary.
%                          the dof's of the faces
% Grid = structure cntaining all relevant information about the grid
% I = identity matrix in the full space
%
% Output:
% B = Nc by N matrix of the Dirichlet constraints
% N = (N-Nc) by (N-Nc) matrix of the nullspace of B
% fn = N by 1 r.h.s. vector of Neuman contributions
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);
% >> Param.dof_dir   = Grid.dof_xmin;    % identify cells on Dirichlet bnd
% >> Param.dof_f_dir = Grid.dof_f_xmin;  % identify faces on Dirichlet bnd
% >> [B,N,fn] = build_bnd(Param,Grid,I);

%% Check input format
if isrow(Param.dof_dir)   && length(Param.dof_dir)>1;   error('Param.dof_dir is not a column vector'); end
if isrow(Param.dof_f_dir) && length(Param.dof_f_dir)>1; error('Param.dof_f_dir is a not column vector'); end

B = I(Param.dof_dir, :);
N = I;
N(:,Param.dof_dir) = [];

fn = spalloc(Grid.Nx,1,0);
