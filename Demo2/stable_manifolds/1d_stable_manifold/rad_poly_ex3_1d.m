function [I,success] = rad_poly_ex3_1d(a,scaling,nu)

N = (length(a)+3)/3; % a = (a0,a1,a2) with ai = (ai_n)_{n=2}^N \in R^{N-1}

a0 = a(1:N-1); a1 = a(N:2*N-2); a2 = a(2*N-1:3*N-3);

a_par = intval('0.3'); c = intval('0.7'); delta = intval(9); w = intval('0.02');

load initial_data_ex3
lambda = tlambda0; v = tv0; tx = ip0;

I = [-1 -1];
success = 0 ;

[f0,f1,f2] = f_1d_stable_manifold_ex3_ext(a,scaling);

f0F = f0(1:N-1); f0_tail = f0(N:end); % size of f0_tail = 4*N
f1F = f1(1:N-1); f1_tail = f1(N:end); % size of f1_tail = 4*N
f2F = f2(1:N-1); f2_tail = f2(N:end); % size of f2_tail = 4*N

disp('Evaluating the jacobian with interval arithmetic ...')
Df = Df_1d_stable_manifold_ex3(a,scaling);
A = inv(mid(Df));

%%%%%%%%%%
%%% Y0 %%%
%%%%%%%%%%

y0_F = A*[f0F;f1F;f2F];

Y0_0 = norma([y0_F(1:N-1);1./(lambda*(N+1:5*N)').*f0_tail],nu);
Y0_1 = norma([y0_F(N:2*N-2);1./(lambda*(N+1:5*N)').*f1_tail],nu);
Y0_2 = norma([y0_F(2*N-1:3*N-3);1./(lambda*(N+1:5*N)').*f2_tail],nu);

Y0 = max([sup(Y0_0) sup(Y0_1) sup(Y0_2)]);

disp(['Y0 = ',num2str(Y0)])

%%%%%%%%%%
%%% Z0 %%%
%%%%%%%%%%

iA = intval(A);
B = intval(eye(3*(N-1)) - iA*Df);

i0 = (1:N-1); i1 = (N:2*N-2); i2 = (2*N-1:3*N-3);

B00 = B(i0,i0); B01 = B(i0,i1); B02 = B(i0,i2); 
B10 = B(i1,i0); B11 = B(i1,i1); B12 = B(i1,i2); 
B20 = B(i2,i0); B21 = B(i2,i1); B22 = B(i2,i2);

norm_B00 = mnorma(B00,nu);
norm_B01 = mnorma(B01,nu);
norm_B02 = mnorma(B02,nu);
norm_B10 = mnorma(B10,nu);
norm_B11 = mnorma(B11,nu);
norm_B12 = mnorma(B12,nu);
norm_B20 = mnorma(B20,nu);
norm_B21 = mnorma(B21,nu);
norm_B22 = mnorma(B22,nu);

Z0_0 = norm_B00 + norm_B01 + norm_B02;
Z0_1 = norm_B10 + norm_B11 + norm_B12;
Z0_2 = norm_B20 + norm_B21 + norm_B22;

Z0 = max([sup(Z0_0) sup(Z0_1) sup(Z0_2)]);

disp(['Z0 = ',num2str(Z0)])

%%%%%%%%%%
%%% Z1 %%%
%%%%%%%%%%

v = scaling*v;

a0 = [tx(1);v(1);a0];  a1 = [tx(2);v(2);a1]; a2 = [tx(3);v(3);a2];

a0_2 = cauchy_ext([a0 a0],4*N+1);
a1_2 = cauchy_ext([a1 a1],4*N+1);
a2_2 = cauchy_ext([a2 a2],4*N+1);
a0a1 = cauchy_ext([a0 a1],4*N+1);
a0a2 = cauchy_ext([a0 a2],4*N+1);

a0_4 = cauchy_ext([a0 a0 a0 a0],4*N+1);
a0_2a1_2 = cauchy_ext([a0 a0 a1 a1],4*N+1);
a0_2a2_2 = cauchy_ext([a0 a0 a2 a2],4*N+1);
a0_2a1a2 = cauchy_ext([a0 a0 a1 a2],4*N+1);
a1_3a2 = cauchy_ext([a1 a1 a1 a2],4*N+1);
a0a1_2a2 = cauchy_ext([a0 a1 a1 a2],4*N+1);
a0_3a2 = cauchy_ext([a0 a0 a0 a2],4*N+1);
a0_3a1 = cauchy_ext([a0 a0 a0 a1],4*N+1);
a0a1_3 = cauchy_ext([a0 a1 a1 a1],4*N+1);
a0a1a2_2 = cauchy_ext([a0 a1 a2 a2],4*N+1);
a1_4 = cauchy_ext([a1 a1 a1 a1],4*N+1);
a0a2_3 = cauchy_ext([a0 a2 a2 a2],4*N+1);
a1_2a2_2 = cauchy_ext([a1 a1 a2 a2],4*N+1);

one = [1;zeros(4*N,1)];

c00 = 9*a0_2 - one + a1_2 + a2_2 - 10*a0_4 - 6*a0_2a1_2 ...
    - 3*(c/delta+2)*a0_2a2_2 - 3*(1+a_par/delta)*a0_2a1a2 - (1/delta)*a1_3a2 ...
    + 2*((a_par+1)/delta)*a0a1_2a2 - 4*(w/delta)*a0_3a2;

c01 = 2*a0a1 - 4*a0_3a1 ...
    - (1+a_par/delta)*a0_3a2 - 3*(1/delta)*a0a1_2a2 ...
    + 2*((a_par+1)/delta)*a0_2a1a2 ;

c02 = 2*a0a2 - 2*(c/delta+2)*a0_3a2 - (1+a_par/delta)*a0_3a1 - (1/delta)*a0a1_3 ...
    + ((a_par+1)/delta)*a0_2a1_2 - (w/delta)*a0_4;

c10 = 4*a0a1 + 2*a0a2 - 8*a0_3a1 - 4*a0a1_3 - 2*(c/delta+2)*a0a1a2_2 ...
    - 2*(1+a_par/delta)*a0a1_2a2 + ((a_par+1)/delta)*a1_3a2 ...
    - 3*(w/delta)*a0_2a1a2;

c11 = 2*a0_2 - 2*a0_4 - 6*a0_2a1_2 - (c/delta+2)*a0_2a2_2 ...
    - 2*(1+a_par/delta)*a0_2a1a2 - 4*(1/delta)*a1_3a2 + 3*((a_par+1)/delta)*a0a1_2a2 ...
    - (w/delta)*a0_3a2;

c12 = a0_2 - 2*(c/delta+2)*a0_2a1a2 ...
    - (1+a_par/delta)*a0_2a1_2 - (1/delta)*a1_4 + ((a_par+1)/delta)*a0a1_3 ...
    - (w/delta)*a0_3a1;

c20 = 2*(2+c/delta)*a0a2 - ((a_par+1)/delta)*a1_2 + 2*(a_par/delta)*a0a1 ...
    + 3*(w/delta)*a0_2 - 8*a0_3a2 - 4*a0a1_2a2 - 2*(c/delta+2)*a0a2_3 ...
    - 2*(1+a_par/delta)*a0a1a2_2 + ((a_par+1)/delta)*a1_2a2_2 ...
    - 3*(w/delta)*a0_2a2_2;

c21 = 3*(1/delta)*a1_2 - 2*((a_par+1)/delta)*a0a1 + (a_par/delta)*a0_2 - 4*a0_2a1a2 ...
    - (1+a_par/delta)*a0_2a2_2 - 3*(1/delta)*a1_2a2_2 + 2*((a_par+1)/delta)*a0a1a2_2;

c22 = (2+c/delta)*a0_2 - 2*a0_4 - 2*a0_2a1_2 - 3*(c/delta+2)*a0_2a2_2 ...
    - 2*(1+a_par/delta)*a0_2a1a2 - 2*(1/delta)*a1_3a2 + 2*((a_par+1)/delta)*a0a1_2a2 ...
    - 2*(w/delta)*a0_3a2;

Z1_0 = (1/(abs(lambda)*(N+1)))*(norma(c00,nu)+norma(c01,nu)+norma(c02,nu));
Z1_1 = (1/(abs(lambda)*(N+1)))*(norma(c10,nu)+norma(c11,nu)+norma(c12,nu));
Z1_2 = (1/(abs(lambda)*(N+1)))*(norma(c20,nu)+norma(c21,nu)+norma(c22,nu));

Z1 = max([sup(Z1_0) sup(Z1_1) sup(Z1_2)]);

disp(['Z1 = ',num2str(Z1)])

%%%%%%%%%%
%%% Z2 %%%
%%%%%%%%%%

A00 = A(i0,i0); A01 = A(i0,i1); A02 = A(i0,i2); 
A10 = A(i1,i0); A11 = A(i1,i1); A12 = A(i1,i2); 
A20 = A(i2,i0); A21 = A(i2,i1); A22 = A(i2,i2);

r_star = 1e-6;

h = infsup(-r_star,r_star)*nu.^-(0:N)';

b0 = a0 + h;
b1 = a1 + h;
b2 = a2 + h;

a = a_par;

b0_3 = cauchy([b0 b0 b0]);
b1_3 = cauchy([b1 b1 b1]);
b2_3 = cauchy([b2 b2 b2]);
b0b1_2 = cauchy([b0 b1 b1]);
b0b1b2 = cauchy([b0 b1 b2]);
b0b2_2 = cauchy([b0 b2 b2]);
b1_2b2 = cauchy([b1 b1 b2]);
b0_2b1 = cauchy([b0 b0 b1]);
b0_2b2 = cauchy([b0 b0 b2]);
b1b2_2 = cauchy([b1 b2 b2]);

nb0 = norma(a0,nu) + r_star;
nb1 = norma(a1,nu) + r_star;
nb2 = norma(a2,nu) + r_star;

% The following computations only require data of size N-1. Hence, Cauchy
% is used and not cauchy_ext for this part. 

g0_x0x0 = 18*b0 - 40*b0_3 - 12*b0b1_2 - 6*b0b1b2 - 12*b0b2_2 - (6*a*b0b1b2)/delta + (2*b1_2b2)/delta + (2*a*b1_2b2)/delta - (6*b0b2_2*c)/delta - (12*b0_2b2*w)/delta;
g0_x0x1 = 2*b1 - 12*b0_2b1 - 3*b0_2b2 - (3*a*b0_2b2)/delta + (4*b0b1b2)/delta + (4*a*b0b1b2)/delta - (3*b1_2b2)/delta;
g0_x0x2 = -3*b0_2b1 + 2*b2 - 12*b0_2b2 - (3*a*b0_2b1)/delta + (2*b0b1_2)/delta + (2*a*b0b1_2)/delta - b1_3/delta - (6*b0_2b2*c)/delta - (4*b0_3*w)/delta;
g0_x1x1 = 2*b0 - 4*b0_3 + (2*b0_2b2)/delta + (2*a*b0_2b2)/delta - (6*b0b1b2)/delta;
g0_x1x2 = -b0_3 - (a*b0_3)/delta + (2*b0_2b1)/delta + (2*a*b0_2b1)/delta - (3*b0b1_2)/delta;
g0_x2x2 = 2*b0 - 4*b0_3 - (2*b0_3*c)/delta;

g1_x0x0 = 4*b1 - 24*b0_2b1 - 4*b1_3 + 2*b2 - 2*b1_2b2 - 4*b1b2_2 - (2*a*b1_2b2)/delta - (2*b1b2_2*c)/delta - (6*b0b1b2*w)/delta;
g1_x0x1 = 4*b0 - 8*b0_3 - 12*b0b1_2 - 4*b0b1b2 - 4*b0b2_2 - (4*a*b0b1b2)/delta + (3*b1_2b2)/delta + (3*a*b1_2b2)/delta - (2*b0b2_2*c)/delta - (3*b0_2b2*w)/delta;
g1_x0x2 = 2*b0 - 2*b0b1_2 - 8*b0b1b2 - (2*a*b0b1_2)/delta + b1_3/delta + (a*b1_3)/delta - (4*b0b1b2*c)/delta - (3*b0_2b1*w)/delta;
g1_x1x1 = -12*b0_2b1 - 2*b0_2b2 - (2*a*b0_2b2)/delta + (6*b0b1b2)/delta + (6*a*b0b1b2)/delta - (12*b1_2b2)/delta;
g1_x1x2 = -2*b0_2b1 - 4*b0_2b2 - (2*a*b0_2b1)/delta + (3*b0b1_2)/delta + (3*a*b0b1_2)/delta - (4*b1_3)/delta - (2*b0_2b2*c)/delta - (b0_3*w)/delta;
g1_x2x2 = -4*b0_2b1 - (2*b0_2b1*c)/delta;

g2_x0x0 = 4*b2 - 24*b0_2b2 - 4*b1_2b2 - 2*b1b2_2 - 4*b2_3 + (2*a*b1)/delta - (2*a*b1b2_2)/delta + (2*b2*c)/delta - (2*b2_3*c)/delta + (6*b0*w)/delta - (6*b0b2_2*w)/delta;
g2_x0x1 = -8*b0b1b2 - 2*b0b2_2 + (2*a*b0)/delta - (2*b1)/delta - (2*a*b1)/delta - (2*a*b0b2_2)/delta + (2*b1b2_2)/delta + (2*a*b1b2_2)/delta;
g2_x0x2 = 4*b0 - 8*b0_3 - 4*b0b1_2 - 4*b0b1b2 - 12*b0b2_2 - (4*a*b0b1b2)/delta + (2*b1_2b2)/delta + (2*a*b1_2b2)/delta + (2*b0*c)/delta - (6*b0b2_2*c)/delta - (6*b0_2b2*w)/delta;
g2_x1x1 = -4*b0_2b2 - (2*b0)/delta - (2*a*b0)/delta + (6*b1)/delta + (2*b0b2_2)/delta + (2*a*b0b2_2)/delta - (6*b1b2_2)/delta;
g2_x1x2 = -4*b0_2b1 - 2*b0_2b2 - (2*a*b0_2b2)/delta + (4*b0b1b2)/delta + (4*a*b0b1b2)/delta - (6*b1_2b2)/delta;
g2_x2x2 = -2*b0_2b1 - 12*b0_2b2 - (2*a*b0_2b1)/delta + (2*b0b1_2)/delta + (2*a*b0b1_2)/delta - (2*b1_3)/delta - (6*b0_2b2*c)/delta - (2*b0_3*w)/delta;

ng0_x0x0 = 18*nb0 + 40*nb0^3 + 12*nb0*nb1^2 + 6*nb0*nb1*nb2 + 12*nb0*nb2^2 + (6*a*nb0*nb1*nb2)/delta + (2*nb1^2*nb2)/delta + (2*a*nb1^2*nb2)/delta + (6*nb0*nb2^2*c)/delta + (12*nb0^2*nb2*w)/delta;
ng0_x0x1 = 2*nb1 + 12*nb0^2*nb1 + 3*nb0^2*nb2 + (3*a*nb0^2*nb2)/delta + (4*nb0*nb1*nb2)/delta + (4*a*nb0*nb1*nb2)/delta + (3*nb1^2*nb2)/delta;
ng0_x0x2 = 3*nb0^2*nb1 + 2*nb2 + 12*nb0^2*nb2 + (3*a*nb0^2*nb1)/delta + (2*nb0*nb1^2)/delta + (2*a*nb0*nb1^2)/delta + nb1^3/delta + (6*nb0^2*nb2*c)/delta + (4*nb0^3*w)/delta;
ng0_x1x1 = 2*nb0 + 4*nb0^3 + (2*nb0^2*nb2)/delta + (2*a*nb0^2*nb2)/delta + (6*nb0*nb1*nb2)/delta;
ng0_x1x2 = nb0^3 + (a*nb0^3)/delta + (2*nb0^2*nb1)/delta + (2*a*nb0^2*nb1)/delta + (3*nb0*nb1^2)/delta;
ng0_x2x2 = 2*nb0 + 4*nb0^3 + (2*nb0^3*c)/delta;

ng1_x0x0 = 4*nb1 + 24*nb0^2*nb1 + 4*nb1^3 + 2*nb2 + 2*nb1^2*nb2 + 4*nb1*nb2^2 + (2*a*nb1^2*nb2)/delta + (2*nb1*nb2^2*c)/delta + (6*nb0*nb1*nb2*w)/delta;
ng1_x0x1 = 4*nb0 + 8*nb0^3 + 12*nb0*nb1^2 + 4*nb0*nb1*nb2 + 4*nb0*nb2^2 + (4*a*nb0*nb1*nb2)/delta + (3*nb1^2*nb2)/delta + (3*a*nb1^2*nb2)/delta + (2*nb0*nb2^2*c)/delta + (3*nb0^2*nb2*w)/delta;
ng1_x0x2 = 2*nb0 + 2*nb0*nb1^2 + 8*nb0*nb1*nb2 + (2*a*nb0*nb1^2)/delta + nb1^3/delta + (a*nb1^3)/delta + (4*nb0*nb1*nb2*c)/delta + (3*nb0^2*nb1*w)/delta;
ng1_x1x1 = 12*nb0^2*nb1 + 2*nb0^2*nb2 + (2*a*nb0^2*nb2)/delta + (6*nb0*nb1*nb2)/delta + (6*a*nb0*nb1*nb2)/delta + (12*nb1^2*nb2)/delta;
ng1_x1x2 = 2*nb0^2*nb1 + 4*nb0^2*nb2 + (2*a*nb0^2*nb1)/delta + (3*nb0*nb1^2)/delta + (3*a*nb0*nb1^2)/delta + (4*nb1^3)/delta + (2*nb0^2*nb2*c)/delta + (nb0^3*w)/delta;
ng1_x2x2 = 4*nb0^2*nb1 + (2*nb0^2*nb1*c)/delta;

ng2_x0x0 = 4*nb2 + 24*nb0^2*nb2 + 4*nb1^2*nb2 + 2*nb1*nb2^2 + 4*nb2^3 + (2*a*nb1)/delta + (2*a*nb1*nb2^2)/delta + (2*nb2*c)/delta + (2*nb2^3*c)/delta + (6*nb0*w)/delta + (6*nb0*nb2^2*w)/delta;
ng2_x0x1 = 8*nb0*nb1*nb2 + 2*nb0*nb2^2 + (2*a*nb0)/delta + (2*nb1)/delta + (2*a*nb1)/delta + (2*a*nb0*nb2^2)/delta + (2*nb1*nb2^2)/delta + (2*a*nb1*nb2^2)/delta;
ng2_x0x2 = 4*nb0 + 8*nb0^3 + 4*nb0*nb1^2 + 4*nb0*nb1*nb2 + 12*nb0*nb2^2 + (4*a*nb0*nb1*nb2)/delta + (2*nb1^2*nb2)/delta + (2*a*nb1^2*nb2)/delta + (2*nb0*c)/delta + (6*nb0*nb2^2*c)/delta + (6*nb0^2*nb2*w)/delta;
ng2_x1x1 = 4*nb0^2*nb2 + (2*nb0)/delta + (2*a*nb0)/delta + (6*nb1)/delta + (2*nb0*nb2^2)/delta + (2*a*nb0*nb2^2)/delta + (6*nb1*nb2^2)/delta;
ng2_x1x2 = 4*nb0^2*nb1 + 2*nb0^2*nb2 + (2*a*nb0^2*nb2)/delta + (4*nb0*nb1*nb2)/delta + (4*a*nb0*nb1*nb2)/delta + (6*nb1^2*nb2)/delta;
ng2_x2x2 = 2*nb0^2*nb1 + 12*nb0^2*nb2 + (2*a*nb0^2*nb1)/delta + (2*nb0*nb1^2)/delta + (2*a*nb0*nb1^2)/delta + (2*nb1^3)/delta + (6*nb0^2*nb2*c)/delta + (2*nb0^3*w)/delta;

Z2_0_F = norma(A00*g0_x0x0(3:N+1) + A01*g1_x0x0(3:N+1) + A02*g2_x0x0(3:N+1),nu) ...
     + 2*norma(A00*g0_x0x1(3:N+1) + A01*g1_x0x1(3:N+1) + A02*g2_x0x1(3:N+1),nu) ...
     + 2*norma(A00*g0_x0x2(3:N+1) + A01*g1_x0x2(3:N+1) + A02*g2_x0x2(3:N+1),nu) ...
     +   norma(A00*g0_x1x1(3:N+1) + A01*g1_x1x1(3:N+1) + A02*g2_x1x1(3:N+1),nu) ...
     + 2*norma(A00*g0_x1x2(3:N+1) + A01*g1_x1x2(3:N+1) + A02*g2_x1x2(3:N+1),nu) ...
     +   norma(A00*g0_x2x2(3:N+1) + A01*g1_x2x2(3:N+1) + A02*g2_x2x2(3:N+1),nu); 
 
Z2_0_tail = (1/(abs(lambda)*(N+1)))*(ng0_x0x0 + 2*ng0_x0x1 + 2*ng0_x0x2 + ng0_x1x1 + 2*ng0_x1x2 + ng0_x2x2);

Z2_0 = Z2_0_F + Z2_0_tail;

Z2_1_F = norma(A10*g0_x0x0(3:N+1) + A11*g1_x0x0(3:N+1) + A12*g2_x0x0(3:N+1),nu) ...
     + 2*norma(A10*g0_x0x1(3:N+1) + A11*g1_x0x1(3:N+1) + A12*g2_x0x1(3:N+1),nu) ...
     + 2*norma(A10*g0_x0x2(3:N+1) + A11*g1_x0x2(3:N+1) + A12*g2_x0x2(3:N+1),nu) ...
     +   norma(A10*g0_x1x1(3:N+1) + A11*g1_x1x1(3:N+1) + A12*g2_x1x1(3:N+1),nu) ...
     + 2*norma(A10*g0_x1x2(3:N+1) + A11*g1_x1x2(3:N+1) + A12*g2_x1x2(3:N+1),nu) ...
     +   norma(A10*g0_x2x2(3:N+1) + A11*g1_x2x2(3:N+1) + A12*g2_x2x2(3:N+1),nu); 
 
Z2_1_tail = (1/(abs(lambda)*(N+1)))*(ng1_x0x0 + 2*ng1_x0x1 + 2*ng1_x0x2 + ng1_x1x1 + 2*ng1_x1x2 + ng1_x2x2);

Z2_1 = Z2_1_F + Z2_1_tail;

Z2_2_F = norma(A20*g0_x0x0(3:N+1) + A21*g1_x0x0(3:N+1) + A22*g2_x0x0(3:N+1),nu) ...
     + 2*norma(A20*g0_x0x1(3:N+1) + A21*g1_x0x1(3:N+1) + A22*g2_x0x1(3:N+1),nu) ...
     + 2*norma(A20*g0_x0x2(3:N+1) + A21*g1_x0x2(3:N+1) + A22*g2_x0x2(3:N+1),nu) ...
     +   norma(A20*g0_x1x1(3:N+1) + A21*g1_x1x1(3:N+1) + A22*g2_x1x1(3:N+1),nu) ...
     + 2*norma(A20*g0_x1x2(3:N+1) + A21*g1_x1x2(3:N+1) + A22*g2_x1x2(3:N+1),nu) ...
     +   norma(A20*g0_x2x2(3:N+1) + A21*g1_x2x2(3:N+1) + A22*g2_x2x2(3:N+1),nu); 
 
Z2_2_tail = (1/(abs(lambda)*(N+1)))*(ng2_x0x0 + 2*ng2_x0x1 + 2*ng2_x0x2 + ng2_x1x1 + 2*ng2_x1x2 + ng2_x2x2);

Z2_2 = Z2_2_F + Z2_2_tail;

Z2 = max([sup(Z2_0) sup(Z2_1) sup(Z2_2)]);

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
      plot_manifold(mid([a0;a1;a2]),mid(tx),nu)

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