function f = f_stable_manifold(a,par)

N = (length(a)+2)/2; % a = (a1,a2) with ai = (ai_n)_{n=2}^N \in R^{N-1}

a1 = a(1:N-1); a2 = a(N:2*N-2);

scaling = par(1); 
rho1 = par(2); 
rho2 = par(3); 
c = par(4);
c1 = par(5);
c2 = par(6);

lambda = -(rho2/2)*(rho2-rho1); % Stable eigenvalue
tx = [rho2;0]; % Fixed point 
v = scaling*[-rho2^2*(c*rho2+c1);-(3*rho2/2)*(rho2-rho1)]; % Stable eigenvector

a1 = [tx(1);v(1);a1]; a2 = [tx(2);v(2);a2];

n = (2:N)';

phi1 = cauchy([a1 a1-[rho1;zeros(N,1)] a1-[rho2;zeros(N,1)]]) - c*cauchy([a1 a1 a1 a2]) - c1*cauchy([a1 a1 a2]);
phi2 = -(1/2)*cauchy([a1 a1 a2]) + (1/2)*rho1*rho2*a2 + c*cauchy([a1 a1 a2 a2]) + c2*cauchy([a1 a1 a2 a2 a2]);

f = [lambda*n.*a1(3:end) - mid(phi1(3:end)) ; lambda*n.*a2(3:end) - mid(phi2(3:end))];

end

