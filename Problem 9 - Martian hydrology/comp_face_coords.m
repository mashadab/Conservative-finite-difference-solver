function [X_faces,Y_faces] = comp_face_coords(dof_faces,Grid)

% Coordinates of the face centers
[Xx,Yx] = meshgrid(Grid.xf,Grid.yc);
[Xy,Yy] = meshgrid(Grid.xc,Grid.yf);
Xf = [Xx(:);Xy(:)]; Yf = [Yx(:);Yy(:)];

% Coordinates of x-face centers
Xfx = Xf(1:Grid.Nfx); 
Yfx = Yf(1:Grid.Nfx);

% Coordinates of y-face centers
Xfy = Xf(Grid.Nfx+1:Grid.Nf);
Yfy = Yf(Grid.Nfx+1:Grid.Nf);

% Center coordinates of bounding faces
dof_x_faces = dof_faces(dof_faces<=Grid.Nfx);
dof_y_faces = dof_faces(dof_faces >Grid.Nfx)-Grid.Nfx;

Xbx = Xfx(dof_x_faces); Ybx = Yfx(dof_x_faces);
Xby = Xfy(dof_y_faces); Yby = Yfy(dof_y_faces);

% Create two matrices size 2 by Nb containing the 
% x and y coords of the endpoints of the faces
X_faces = [[Xbx' Xby'+Grid.dx/2];[Xbx' Xby'-Grid.dx/2]];
Y_faces = [[Ybx'+Grid.dy/2 Yby'];[Ybx'-Grid.dy/2 Yby']];
