function [Grid] = build_grid(Grid)
% Author: Mohammad Afzal Shadab
% Date: 02/02/2021
% Description:
% This function takes in minimal definition of the computational
% domain and grid and computes all containing all pertinent information 
% about the grid. 
% Input:
% Grid.xmin = left boundary of the domain
% Grid.xmax = right bondary of the domain
% Grid.Nx   = number of grid cells
% 
% Output: (suggestions)
% Grid.Lx = scalar length of the domain
% Grid.dx = scalar cell width
% Grid.Nfx = number of fluxes in x-direction
% Grid.xc = Nx by 1 column vector of cell center locations
% Grid.xf = Nfx by 1 column vector of cell face locations
% Grid.dof = Nx by 1 column vector from 1 to N containig the degrees of freedom, i.e. cell numbers
% Grid.dof_f = Nfx by 1 column vector from 1 to Nfx containig the degrees of freedom, i.e. face numbers
% Grid.dof_xmin  = scalar cell degree of freedom corrsponding to the left boundary
% Grid.dof_xmax  = scalar cell degree of freedom corrsponding to the right boundary
% Grid.dof_f_xmin = scalar face degree of freedom corrsponding to the left boundary
% Grid.dof_f_xmax = scalar face degree of freedom corrsponding to the right boundary
% + anything else you might find useful
%
% Example call: 
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10; 
% >> Grid = build_grid(Grid);

%% Set up catesian geometry
if ~isfield(Grid,'xmin'); Grid.xmin = 0;  fprintf('Grid.xmin is not defined and has been set to zero.\n');end
if ~isfield(Grid,'xmax'); Grid.xmax = 10; fprintf('Grid.xmax is not defined and has been set to 10.\n'); end
if ~isfield(Grid,'Nx');   Grid.Nx   = 10; fprintf('Grid.Nx is not defined and has been set to 10.\n');end
Grid.Lx = Grid.xmax - Grid.xmin;    % domain length in x
Grid.dx = Grid.Lx / Grid.Nx;        % dx of the gridblocks

%% Number for fluxes
Grid.Nfx = Grid.Nx + 1;

% Set up mesh
% cell centers 'xc' and cell faces 'xf'   
Grid.xc = (linspace(Grid.xmin + 0.5*Grid.dx, Grid.xmax - 0.5*Grid.dx, Grid.Nx))'; % x-coords of gridblock centers
Grid.xf = (linspace(Grid.xmin, Grid.xmax,Grid.Nfx))'; % x-coords of gridblock faces

%% Set up dof vectors
Grid.dof  = (linspace(1, Grid.Nx, Grid.Nx));% cell centered degree of freedom/gridblock number
Grid.dof_f= (linspace(1, Grid.Nfx, Grid.Nfx));% face degree of freedom/face number

%% Boundary dof's
% Boundary cells
Grid.dof_xmin = 1;
Grid.dof_xmax = Grid.Nx;
% Boundary faces
Grid.dof_f_xmin = 1;
Grid.dof_f_xmax = Grid.Nfx;