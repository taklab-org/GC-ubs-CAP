function Df = Df_stable_manifold(a,par)

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

a1a1 = cauchy([a1 a1]);
a1a2 = cauchy([a1 a2]);
a1a1a1 = cauchy([a1 a1 a1]);
a1a1a2 = cauchy([a1 a1 a2]);
a1a2a2 = cauchy([a1 a2 a2]);
a1a2a2a2 = cauchy([a1 a2 a2 a2]);
a1a1a2a2 = cauchy([a1 a1 a2 a2]);

dphi_11 = intval(zeros(N-1));
dphi_12 = intval(zeros(N-1));
dphi_21 = intval(zeros(N-1));
dphi_22 = intval(zeros(N-1));

for n = 2:N
    j = (2:n);
    k = n-j+1;
    dphi_11(n-1,j-1) = 3*a1a1(k) - 2*(rho1+rho2)*a1(k) - 3*c*a1a1a2(k) - 2*c1*a1a2(k);
    dphi_12(n-1,j-1) = -c*a1a1a1(k) - c1*a1a1(k);
    dphi_21(n-1,j-1) = -a1a2(k) + 2*c*a1a2a2(k) + 2*c2*a1a2a2a2(k);
    dphi_22(n-1,j-1) = -(1/2)*a1a1(k) + 2*c*a1a1a2(k) + 3*c2*a1a1a2a2(k);
end

dphi_11 = dphi_11 + rho1*rho2*eye(N-1);
dphi_22 = dphi_22 + (1/2)*rho1*rho2*eye(N-1);

n = (2:N)';

Df_11 = lambda*diag(n) - dphi_11 ;
Df_12 = - dphi_12 ;
Df_21 = - dphi_21 ;
Df_22 = lambda*diag(n) - dphi_22 ;

Df = [[Df_11 Df_12];[Df_21 Df_22]];

end

