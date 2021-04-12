% Newton-Raphson method for solution of non-linear system of equations
% Author: Mohammad Afzal Shadab
% Date: 04/11/2021
% Email: mashadab@utexas.edu

% Parameters for the Newton iteration
tol = 1e-6; 
kmax = 8;

% Root
U = [1;2;3]; 

% Residual vector
r = @(u) [u(1)^2+3*u(1)*u(2)+5*u(2)*u(3) ;  u(1)*u(2)^2+4*u(1)*u(2)*u(3) ; u(1)^3-5*u(2)^2+3*u(3)^2];
r = @(u) r(u)-r(U); % generates a residual with U as solution
% Jacobian matrix
J = @(u) [2*u(1)+3*u(2), 3*u(1)+5*u(3),5*u(2);...
          u(2)^2+4*u(2)*u(3), 2*u(1)*u(2)+4*u(1)*u(3), 4*u(1)*u(2);...
          3*u(1)^2, -10*u(2), 6*u(3)];
    
% Newton's method with numerical Jacobian
u = [2;2;2]; % initial guess
k = 1; nres = 1; ndu = 1; % initialize Newton
while (nres > tol || ndu > tol) && k <= kmax
    du = -(J(u)\r(u));
    u = u + du; 
    nres = norm(r(u));
    ndu = norm(du);
    fprintf('%d: res = %3.2e, du = %3.2e;\n',k,nres,ndu)
    k = k+1;
end
if k <= kmax+1 && nres <= tol && ndu <= tol
    fprintf('\nNewton iteration converged to tol = %3.2e\n',tol)
else
    fprintf('\nNewton iteration did NOT converge to tol = %3.2e in %d iterations!\n',tol,k)
end
u
