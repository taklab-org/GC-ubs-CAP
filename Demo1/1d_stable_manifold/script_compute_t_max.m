clear
clc
close all

load sol

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

a = [a1;a2];
nu = 1.01;
r0 = 4.2e-13;

t_max = rigorous_computation_tmax(a,rho1,rho2,r0,nu);
