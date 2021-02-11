Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
Grid = build_grid(Grid);
[D,G,I]=build_ops(Grid);
Param.dof_dir   = Grid.dof_xmin; 
Param.dof_f_dir = Grid.dof_f_xmin;
[B,N,fn] = build_bnd(Param,Grid,I)