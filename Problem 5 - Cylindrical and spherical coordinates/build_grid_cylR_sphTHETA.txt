function [Grid] = build_grid(Grid)
% author: Mohammad Afzal Shadab
% date: 17 March 2021
% Description:
% This function computes takes in minimal definition of the computational
% domain and grid and computes all containing all pertinent information 
% about the grid. 
% Input:
% Grid.xmin = left boundary of the domain
% Grid.xmax = right bondary of the domain
% Grid.Nx   = number of grid cells
% Grid.geom = coordinate system

% Output: (suggestions)
% Grid.Lx = scalar length of the domain
% Grid.dx = scalar cell width
% Grid.Nfx = number of fluxes in x-direction
% Grid.xc = Nx by 1 column vector of cell center locations
% Grid.xf = Nfx by 1 column vector of cell face locations
% Grid.dof = Nx by 1 column vector from 1 to N containig the degrees of freedom, i.e. cell numbers
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
Grid.Lx = Grid.xmax-Grid.xmin;    % domain length in x
Grid.dx = Grid.Lx/Grid.Nx;        % dx of the gridblocks

%% Number for fluxes
Grid.Nfx = (Grid.Nx+1);

% Set up mesh
% cell centers 'xc' and cell faces 'xf'   
Grid.xc = [Grid.xmin+Grid.dx/2:Grid.dx:Grid.xmax-Grid.dx/2]'; % x-coords of gridblock centers
Grid.xf = [Grid.xmin:Grid.dx:Grid.xmax]'; % x-coords of gridblock faces

%% Set up dof vectors
Grid.dof = 1:Grid.Nx;              % cell centered degree of freedom/gridblock number
Grid.dof_f = 1:Grid.Nfx;            % face degree of freedom/face number

%% Boundary dof's
% Boundary cells
Grid.dof_xmin = Grid.dof(1);
Grid.dof_xmax = Grid.dof(Grid.Nx);
% Boundary faces
Grid.dof_f_xmin = Grid.dof_f(1);
Grid.dof_f_xmax = Grid.dof_f(Grid.Nfx);

%% Cell volumes and face areas
Grid.dy = 1;
Grid.dz = 1;
switch Grid.geom
    case 'cartesian' 
        % assumes unit crossesctional area (dy=dz=1)
        Grid.A = ones(Grid.Nfx,1);          % cross-sectional area of the cell faces
        Grid.V = Grid.dx*ones(Grid.Nx,1);   % volume of the cells
    case 'cylindrical_r'
        % assumes x = r
        Grid.A = 2*pi*Grid.xf;   % cross-sectional area of the cell faces
        Grid.V = 2*pi*Grid.xc*Grid.dx;% volume of the cells
    case 'spherical_shell'
        % assumes x = theta = co-lattitude
        if ~isfield(Grid,'R_shell')
            fprintf('Grid.R_shell not initialized. Radius of spherical shell initialized to unity.\n ')
            Grid.R_shell = 1;  
        end

        Grid.A = 2*pi*Grid.R_shell.*sin(Grid.xf);   % cross-sectional area of the cell faces
        Grid.V =-2*pi*Grid.R_shell.^2.*(cos(Grid.xc+0.5*Grid.dx)-cos(Grid.xc-0.5*Grid.dx)) ;   % volume of the cells
    otherwise
        error('Unknown grid geometry.')
end