function [f1,f2] = f_stable_manifold_ext(a,par)

N = (length(a)+2)/2; % a = (a1,a2) with ai = (ai_n)_{n=2}^N \in R^{N-1}

a1 = a(1:N-1); a2 = a(N:2*N-2);

scaling = intval(par(1)); 
rho1 = intval(par(2)); 
rho2 = intval(par(3)); 
c = intval(par(4));
c1 = intval(par(5));
c2 = intval(par(6));

lambda = -(rho2/2)*(rho2-rho1); % Stable eigenvalue
tx = [rho2;0]; % Fixed point 
v = scaling*[-rho2^2*(c*rho2+c1);-(3*rho2/2)*(rho2-rho1)]; % Stable eigenvector

a1 = [tx(1);v(1);a1]; a2 = [tx(2);v(2);a2];

n = (2:N)';

a1a1rho1a1rho2 = cauchy_ext([a1 a1-[rho1;zeros(N,1)] a1-[rho2;zeros(N,1)]]);
a1a1a1a2 = cauchy_ext([a1 a1 a1 a2]);
a1a1a2 = cauchy_ext([a1 a1 a2]);
a1a1a2a2 = cauchy_ext([a1 a1 a2 a2]);
a1a1a2a2a2 = cauchy_ext([a1 a1 a2 a2 a2]);

phi1 = [a1a1rho1a1rho2;zeros(N,1)] - c*a1a1a1a2 - c1*[a1a1a2;zeros(N,1)];
phi2 = -(1/2)*[a1a1a2;zeros(2*N,1)] + (1/2)*rho1*rho2*[a2;zeros(4*N,1)] + c*[a1a1a2a2;zeros(N,1)] + c2*a1a1a2a2a2;

f1 = [lambda*n.*a1(3:end);zeros(3*N,1)] - phi1(3:end); 
f2 = [lambda*n.*a2(3:end);zeros(4*N,1)] - phi2(3:end);

end

