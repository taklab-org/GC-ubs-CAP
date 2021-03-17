function [I,success] = rad_poly(a,par,nu)

N = (length(a)+2)/2; % a = (a1,a2) with ai = (ai_n)_{n=2}^N \in R^{N-1}

a1 = a(1:N-1); a2 = a(N:2*N-2);

%%% Small comment: the quantities c, c1, c2, rho1, rho2, lambda, tx & v
%%% sould all be provided with rigorous error bounds (i.e. in interval form)
%%% This will be trivial to do, but still needs to be done.

scaling = par(1);
rho1 = par(2);
rho2 = par(3);
c = par(4);
c1 = par(5);
c2 = par(6);

lambda = -(rho2/2)*(rho2-rho1); % Stable eigenvalue
tx = [rho2;0]; % Fixed point
v = scaling*[-rho2^2*(c*rho2+c1);-(3*rho2/2)*(rho2-rho1)]; % Stable eigenvector

I = [-1 -1];
success = 0 ;

[f1,f2] = f_stable_manifold_ext(a,par);

f1F = f1(1:N-1); f1_tail = f1(N:end); % size of f1_tail = 3*N
f2F = f2(1:N-1); f2_tail = f2(N:end); % size of f2_tail = 4*N

disp('Evaluating the jacobian with interval arithmetic ...')
Df = Df_stable_manifold(a,par);
A = inv(mid(Df));

%%%%%%%%%%
%%% Y0 %%%
%%%%%%%%%%

y0_F = A*[f1F;f2F];

Y0_1 = norma([y0_F(1:N-1);1./(lambda*(N+1:4*N)').*f1_tail],nu);
Y0_2 = norma([y0_F(N:2*N-2);1./(lambda*(N+1:5*N)').*f2_tail],nu);

Y0 = max([sup(Y0_1) sup(Y0_2)]);

disp(['Y0 = ',num2str(Y0)])

%%%%%%%%%%
%%% Z0 %%%
%%%%%%%%%%

iA = intval(A);
B = intval(eye(2*(N-1)) - iA*Df);

i1 = (1:N-1); i2 = (N:2*N-2);

B11 = B(i1,i1); B12 = B(i1,i2); B21 = B(i2,i1); B22 = B(i2,i2);

norm_B11 = mnorma(B11,nu);
norm_B12 = mnorma(B12,nu);
norm_B21 = mnorma(B21,nu);
norm_B22 = mnorma(B22,nu);

Z0_1 = norm_B11 + norm_B12;
Z0_2 = norm_B21 + norm_B22;

Z0 = max([sup(Z0_1) sup(Z0_2)]);

disp(['Z0 = ',num2str(Z0)])

%%%%%%%%%%
%%% Z1 %%%
%%%%%%%%%%

a1 = [tx(1);v(1);a1]; a2 = [tx(2);v(2);a2];

a1a1 = cauchy_ext([a1 a1]);
a1a2 = cauchy_ext([a1 a2]);
a1a1a1 = cauchy_ext([a1 a1 a1]);
a1a1a2 = cauchy_ext([a1 a1 a2]);
a1a2a2 = cauchy_ext([a1 a2 a2]);
a1a2a2a2 = cauchy_ext([a1 a2 a2 a2]);
a1a1a2a2 = cauchy_ext([a1 a1 a2 a2]);

z11 = norma(3*[a1a1;zeros(N,1)] - 2*(rho1+rho2)*[a1;zeros(2*N,1)] - 3*c*a1a1a2 - 2*c1*[a1a2;zeros(N,1)],nu);
z12 = norma(c*a1a1a1 + c1*[a1a1;zeros(N,1)],nu);
Z1_1 = (1/(abs(lambda)*(N+1)))*(z11 + rho1*rho2 + z12);

z21 =  norma(-[a1a2;zeros(2*N,1)] + 2*c*[a1a2a2;zeros(N,1)] + 2*c2*a1a2a2a2,nu);
z22 = norma(-(1/2)*[a1a1;zeros(2*N,1)] + 2*c*[a1a1a2;zeros(N,1)] + 3*c2*a1a1a2a2,nu);
Z1_2 = (1/(abs(lambda)*(N+1)))*(z21 + z22 + rho1*rho2/2 );

Z1 = max([sup(Z1_1) sup(Z1_2)]);

disp(['Z1 = ',num2str(Z1)])

%%%%%%%%%%
%%% Z2 %%%
%%%%%%%%%%

A11 = A(i1,i1); A12 = A(i1,i2); A21 = A(i2,i1); A22 = A(i2,i2);

r_star = 1e-6;

b1 = a1 + infsup(-r_star,r_star)*nu.^-(0:N)';
b2 = a2 + infsup(-r_star,r_star)*nu.^-(0:N)';

b1b2 = cauchy([b1 b2]);
b1b1 = cauchy([b1 b1]);
b2b2 = cauchy([b2 b2]);
b2b2b2 = cauchy([b2 b2 b2]);
b1b1b2 = cauchy([b1 b1 b2]);
b1b2b2 = cauchy([b1 b2 b2]);

nb1 = norma(a1,nu) + r_star;
nb2 = norma(a2,nu) + r_star;

f1_x1x1 = 6*b1 - 6*c*b1b2 - 2*c1*b2;
f1_x1x2 = -3*c*b1b1 - 2*c1*b1;
f2_x1x1 = -b2 + 2*c*b2b2 + 2*c2*b2b2b2;
f2_x1x2 = -b1 + 4*c1*b1b2 + 6*c2*b1b2b2;
f2_x2x2 = 2*c*b1b1 + 6*c2*b1b1b2;

nf1_x1x1 = 6*nb1 + 6*c*nb1*nb2 + 2*c1*nb2 + 2*(rho1+rho2);
nf1_x1x2 = 3*c*nb1*nb1 + 2*c1*nb1;
nf2_x1x1 = nb2 + 2*c*nb2*nb2 + 2*c2*nb2*nb2*nb2;
nf2_x1x2 = nb1 + 4*c1*nb1*nb2 + 6*c2*nb1*nb2*nb2;
nf2_x2x2 = 2*c*nb1*nb1 + 6*c2*nb1*nb1*nb2;

Z2_1_F = norma(A11*f1_x1x1(3:N+1) + A12*f2_x1x1(3:N+1),nu) + 2*(rho1+rho2)*mnorma(A11,nu) ...
    + 2*norma(A11*f1_x1x2(3:N+1) + A12*f2_x1x2(3:N+1),nu) ...
    + norma(A12*f2_x2x2(3:N+1),nu);

Z2_1_tail = (1/(lambda*(N+1)))*(nf1_x1x1 + 2*nf1_x1x2 );

Z2_1 = Z2_1_F + Z2_1_tail;

Z2_2_F = norma(A21*f1_x1x1(3:N+1) + A22*f2_x1x1(3:N+1),nu) + 2*(rho1+rho2)*mnorma(A21,nu) ...
    + 2*norma(A21*f1_x1x2(3:N+1) + A22*f2_x1x2(3:N+1),nu) ...
    + norma(A22*f2_x2x2(3:N+1),nu);

Z2_2_tail = (1/(lambda*(N+1)))*(nf2_x1x1 + 2*nf2_x1x2 + nf2_x2x2);

Z2_2 = Z2_2_F + Z2_2_tail;

Z2 = max([sup(Z2_1) sup(Z2_2)]);

disp(['Z2 = ',num2str(Z2)])

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Radii polynomial %%%
%%%%%%%%%%%%%%%%%%%%%%%%

Y0 = intval(Y0); Z0 = intval(Z0); Z1 = intval(Z1); Z2 = intval(Z2);

if inf(1-Z0-Z1)>0 
  if inf((1-Z0-Z1)^2-4*Y0*Z2) > 0  
    rmin=sup(((1-Z0-Z1) - sqrt((1-Z0-Z1)^2-4*Y0*Z2))/(2*Z2));
    rmax=inf(((1-Z0-Z1) + sqrt((1-Z0-Z1)^2-4*Y0*Z2))/(2*Z2));
    if rmin<rmax 
      success=1;
      I=[rmin rmax];
      disp('success')
      plot_manifold([a1;a2],rho2,nu)

    else
      disp('failure: rmin > rmax')
    end
  else
    disp('failure: discriminant is negative')  
  end
else
    disp('failure: linear term is positive')
end


end

