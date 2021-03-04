function [B,N,fn] = build_bnd(BC,Grid,I)
% author: Mohammad Afzal Shadab
% date: 03/03/2021
% Description:
% This function computes the operators and r.h.s vectors for both Dirichlet
% and Neumann boundary conditions. 
% Note:
% N is not created from I the same way B is created from I, because 
% the vector dof_dir contains the columns that must be eliminated rather
% then the colums that are retained in N. If you wanted to do it this way
% you would have to create a new vector
% dof_non_dir = setdiff(dof,dof_dir)
% I suspect that the set-operators are expensive on large vectors, hence
% we simply eliminate the rows.
%
% Input:
% BC = structure containing all information about the physical problem
%         in particular this function needs the fields
%         BC.dof_dir = Nc by 1 column vector containing 
%                         the dof's of the Dirichlet boundary.
%         BC.dof_neu = N by 1 column vector containing 
%                         the dof's of the Neumann boundary.
%         BC.qb      = column vector of prescribed fluxes on Neuman bnd.
% Grid = structure containing all pertinent information about the grid.
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
% >> BC.dof_dir   = Grid.dof_xmin;    % identify cells on Dirichlet bnd
% >> BC.dof_f_dir = Grid.dof_f_xmin;  % identify faces on Dirichlet bnd
% >> BC.dof_neu   = Grid.dof_xmax;    % identify cells on Neumann bnd
% >> BC.dof_f_neu = Grid.dof_f_xmax;  % identify cells on Neumann bnd
% >> BC.qb = 1;                   % set bnd flux
% >> [B,N,fn] = build_bnd(BC,Grid,I);

%% Check input format
if isrow(BC.dof_dir)   && length(BC.dof_dir)>1;   error('BC.dof_dir is not a column vector'); end
if isrow(BC.dof_neu)   && length(BC.dof_neu)>1;   error('BC.dof_neu is not a column vector'); end
if isrow(BC.dof_f_dir) && length(BC.dof_f_dir)>1; error('BC.dof_f_dir is a not column vector'); end
if isrow(BC.dof_f_neu) && length(BC.dof_f_neu)>1; error('BC.dof_f_neu is a not column vector'); end
if isfield(BC,'qb') && isrow(BC.qb) && length(BC.qb)>1;        error('BC.qb is not a column vector'); end

B = I(BC.dof_dir,:);
N = I; N(:,BC.dof_dir) = [];

%% Neumann boundary conditions
if isempty(BC.dof_neu)
    fn = spalloc(Grid.N,1,0);         % allocate sparse zero vector
else
    fn = spalloc(Grid.N,1,length(BC.dof_neu)); % allocate sparse vector with appropriate non-zero elements
    fn(BC.dof_neu) = BC.qb.*Grid.A(BC.dof_f_neu)./Grid.V(BC.dof_neu);
end
