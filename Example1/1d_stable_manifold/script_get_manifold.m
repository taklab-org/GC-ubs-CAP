close all
clear 

scaling = .033; 
rho1 = 1; 
rho2 = 2; 

xL = [1.9;0.25];
beta_L = 1.9; 
v_L = 4;
xR = [1.5;0.2];
beta_R = 1.5;
v_R = 5;

c = (v_R*B1(beta_R,rho1,rho2)-v_L*B1(beta_L,rho1,rho2))/(beta_R-beta_L);

c1_L = v_L*B1(beta_L,rho1,rho2)-c*beta_L; c2_L = v_L^2*B2(beta_L,rho1,rho2)-c*v_L;
c1_R = v_R*B1(beta_R,rho1,rho2)-c*beta_R; c2_R = v_R^2*B2(beta_R,rho1,rho2)-c*v_R;

c1 = c1_L;
c2 = c2_L;

par = [scaling;rho1;rho2;c;c1;c2];

N = 50; 

a = .0001*rand(2*N-2,1);
a = newton(a,par);

a = padding(a,100); a = newton(a,par);
a = padding(a,100); a = newton(a,par);
a = padding(a,50); a = newton(a,par);

N = (length(a)+2)/2; a1 = a(1:N-1); a2 = a(N:2*N-2);

lambda = -(rho2/2)*(rho2-rho1); % Stable eigenvalue
tx = [rho2;0]; % Fixed point 
v = scaling*[-rho2^2*(c*rho2+c1);-(3*rho2/2)*(rho2-rho1)]; % Stable eigenvector

nu = 1;

plot_manifold([tx(1);v(1);a1;tx(2);v(2);a2],rho2,nu)

