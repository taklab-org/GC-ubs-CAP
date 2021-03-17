% close all
clear 
clc

scaling = .1; 

par = [scaling;1];

N = 45; 

a1 = .0001*rand(N-1,1); a2 = zeros(N-1,1);

a = [a1;a2];
a = newton(a,par);

N = (length(a)+2)/2; 
%i3 = intval(3); i4 = intval(4); x1 = (1/(13-6*sqrt(i3)))^(1/i4); x2 = sqrt(3-sqrt(i3));

x1 = (1/(13-6*sqrt(3)))^(1/4); x2 = sqrt(3-sqrt(3));

lambda = x2^3*( 1/24 + 5*x1^4*(1/8*x2^4 - x2^2 + 23/24 ) );

tx = [x1;x2]; 
v = scaling*[1;0]; % Stable eigenvector

nu = 1;

a1 = a(1:N-1); a2 = a(N:2*N-2);

plot_manifold([tx(1);v(1);a1;tx(2);v(2);a2],tx,nu)

a1 = [tx(1);v(1);a1]
