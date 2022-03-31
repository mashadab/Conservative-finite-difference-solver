% Richards plotting
clc, close all, clear
Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 30;
Grid.ymin = 0; Grid.ymax = 1; Grid.Ny = 30;
Grid = build_grid2D(Grid);
[D,G,C,I,M] = build_ops2D(Grid);
[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);

% suppose our saturated region looks like an elephant!
[x,y] = elephant(.5,.5,'l',.5,'k');
ind_sat = inpolygon(Xc(:),Yc(:),x,y);
dof_sat = Grid.dof(ind_sat);

% find faces
[dof_f_bnd,dof_f] = find_faces(dof_sat,D,Grid);
[X_sat,Y_sat] = comp_face_coords(dof_f,Grid);
[X_bnd,Y_bnd] = comp_face_coords(dof_f_bnd,Grid);

subplot 121
plot(x,y), hold on
% plot(Xc(dof_sat),Yc(dof_sat),'r.')
plot([0 1 1 0 0],[0 0 1 1 0],'k-')
axis equal tight
title 'actual domain'

subplot 122
% basic grid
plot([Grid.xf';Grid.xf'],[ones(1,Grid.Nx+1); zeros(1,Grid.Nx+1)],'-','color',0.5*[1 1 1]), hold on
plot([ones(1,Grid.Ny+1); zeros(1,Grid.Ny+1)],[Grid.yf';Grid.yf'],'-','color',0.5*[1 1 1])
plot(X_sat,Y_sat,'r-')
plot(X_bnd,Y_bnd,'r-','linewidth',2)
axis equal tight
title 'something like this'



function [x,y] = elephant(x0,y0,dir,scale,col)
% author: Marc Hesse
% date: 17 Jan 2018
% Description: Plots an elephant with four complex parameters
%              following the paper by Meyer et al. (2010),
%              Am. J. Phys., 78(?6), p. 648-649
% Input: x0, y0 = scalars that define center of mass
%        dir = string that is either 'l' or 'r' and
%              defined the direction the elephant is looking
%        scale = elephant height
%        col = 1 by 3 row vector defining the color
% Output: None
% Example function call:
%
% plot_elephant(0,0,'l',1,[0 0 0])

Ax = [-60  0   0 0   0];
Bx = [-30  8 -10 0   0];
Ay = [  0  0  12 0 -14];
By = [ 50 18   0 0   0];

t = linspace(0,2*pi,1e2);
x = 0*t; y = x;
for k = 1:5
    x = x + Ax(k)*cos(k*t) + Bx(k)*sin(k*t);
    y = y - Ay(k)*cos(k*t) - By(k)*sin(k*t);
end
% rescale 
Dy = max(y) - min(y);
y = scale*y/Dy; x = scale*x/Dy;
% reflect
if strcmp(dir,'l')
    x = -x;
end
% shift
x = x0 + x; y = y0 + y;


end
