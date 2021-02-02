Grid.xmin = 0; Grid.xmax = 3; Grid.Nx = 30;
Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid); L = D*G;

xa = linspace(Grid.xmin,Grid.xmax,3e2);
f = @(x) exp(cos(2*pi*x));
dfdx = @(x) -exp(cos(2*pi*x)).*sin(2*pi*x)*2*pi; 
d2fdx2 = @(x) -2*pi^2*exp(cos(2*pi*x)).*(2*cos(2*pi*x)+cos(4*pi*x)-1);

subplot 311
plot(xa,f(xa),'r',Grid.xc,f(Grid.xc),'bo')
xlabel 'x'
ylabel 'f'
legend('analytical','numerical')

subplot 312
plot(xa,dfdx(xa),'r',Grid.xf,G*f(Grid.xc),'bo')
xlabel 'x'
ylabel 'df/dx'
legend('analytical','numerical')

subplot 313
plot(xa,d2fdx2(xa),'r',Grid.xc,L*f(Grid.xc),'bo')
xlabel 'x'
ylabel 'd^2f/dx^2'
legend('analytical','numerical')