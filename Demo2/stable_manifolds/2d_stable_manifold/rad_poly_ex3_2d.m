function [I,success] = rad_poly_ex3_2d(a,scaling,fp)

N = round((-3+sqrt(25+8/3*length(a)))/2,0); % a = (a0,a1,a2) with ai = (ai_{n})_{|n|=2}^N \in R^{(N+1)*(N+2)/2-3}

i0 = (1:(N+1)*(N+2)/2-3); i1 = i0+(N+1)*(N+2)/2-3; i2 = i1+(N+1)*(N+2)/2-3;
a0 = a(i0); a1 = a(i1); a2 = a(i2);

a_par = intval('0.3'); c = intval('0.7'); delta = intval(9); w = intval('0.02');

I = [-1 -1];
success = 0 ;

[f0,f1,f2] = f_2d_stable_manifold_ex3_ext(a,scaling,fp);

f0F = f0(1:N+1,1:N+1); f0F = reshape_data_2d_manifold_matrix2vec(f0F);
f1F = f1(1:N+1,1:N+1); f1F = reshape_data_2d_manifold_matrix2vec(f1F);
f2F = f2(1:N+1,1:N+1); f2F = reshape_data_2d_manifold_matrix2vec(f2F);

disp('Evaluating the jacobian with interval arithmetic ...')
Df = Df_2d_stable_manifold_ex3(a,scaling,fp);
A = inv(mid(Df));

load initial_data_ex3

if fp == 1   
    lambda = [tlambda1_1;tlambda1_2]; 
    v1 = scaling(1)*tv1_1; 
    v2 = scaling(2)*tv1_2;
    tx = ip1;
end

if fp == 2
    lambda = [tlambda2_1 tlambda2_2]; 
    v1 = scaling(1)*tv2_1; 
    v2 = scaling(2)*tv2_2;
    tx = ip2;
end

%%%%%%%%%%
%%% Y0 %%%
%%%%%%%%%%

y0_F = A*[f0F;f1F;f2F];

y0_0F = y0_F(i0); y0_1F = y0_F(i1); y0_2F = y0_F(i2);
y0_0F_matrix = reshape_data_2d_manifold_vec2matrix(y0_0F,0,0,0);
y0_1F_matrix = reshape_data_2d_manifold_vec2matrix(y0_1F,0,0,0);
y0_2F_matrix = reshape_data_2d_manifold_vec2matrix(y0_2F,0,0,0);

Y0_F = sum(sum(abs(y0_0F_matrix)));
Y1_F = sum(sum(abs(y0_1F_matrix)));
Y2_F = sum(sum(abs(y0_2F_matrix)));

matrix_n1 = repmat((0:5*N)',1,5*N+1);
matrix_n2 = repmat((0:5*N),5*N+1,1);
LAMBDA_EXT = lambda(1)*matrix_n1 + lambda(2)*matrix_n2;
LAMBDA_EXT(1,1)=1; % In order not to divide by zero
f0(1:N+1,1:N+1)=0; f1(1:N+1,1:N+1)=0; f2(1:N+1,1:N+1)=0;
Y0_tail = sum(sum(abs(LAMBDA_EXT.^(-1).*f0)));
Y1_tail = sum(sum(abs(LAMBDA_EXT.^(-1).*f1)));
Y2_tail = sum(sum(abs(LAMBDA_EXT.^(-1).*f2)));
        
Y0_0 = Y0_F + Y0_tail;
Y0_1 = Y1_F + Y1_tail;
Y0_2 = Y2_F + Y2_tail;

Y0 = max([sup(Y0_0) sup(Y0_1) sup(Y0_2)]);

disp(['Y0 = ',num2str(Y0)])

%%%%%%%%%%
%%% Z0 %%%
%%%%%%%%%%

iA = intval(A);
B = intval(eye(length(A)) - iA*Df);

B00 = B(i0,i0); B01 = B(i0,i1); B02 = B(i0,i2); 
B10 = B(i1,i0); B11 = B(i1,i1); B12 = B(i1,i2); 
B20 = B(i2,i0); B21 = B(i2,i1); B22 = B(i2,i2);

norm_B00 = norm(B00,1);
norm_B01 = norm(B01,1);
norm_B02 = norm(B02,1);
norm_B10 = norm(B10,1);
norm_B11 = norm(B11,1);
norm_B12 = norm(B12,1);
norm_B20 = norm(B20,1);
norm_B21 = norm(B21,1);
norm_B22 = norm(B22,1);

Z0_0 = norm_B00 + norm_B01 + norm_B02;
Z0_1 = norm_B10 + norm_B11 + norm_B12;
Z0_2 = norm_B20 + norm_B21 + norm_B22;

Z0 = max([sup(Z0_0) sup(Z0_1) sup(Z0_2)]);

disp(['Z0 = ',num2str(Z0)])

%%%%%%%%%%
%%% Z1 %%%
%%%%%%%%%%

a0 = reshape_data_2d_manifold_vec2matrix(a0,tx(1),v1(1),v2(1));
a1 = reshape_data_2d_manifold_vec2matrix(a1,tx(2),v1(2),v2(2));
a2 = reshape_data_2d_manifold_vec2matrix(a2,tx(3),v1(3),v2(3));

a0_2 = padmat(cauchy2d_ext(a0,a0),2*N);
a1_2 = padmat(cauchy2d_ext(a1,a1),2*N);
a2_2 = padmat(cauchy2d_ext(a2,a2),2*N);
a0a1 = padmat(cauchy2d_ext(a0,a1),2*N);
a0a2 = padmat(cauchy2d_ext(a0,a2),2*N);

a0_4 = cauchy2d_ext(a0,a0,a0,a0);
a0_2a1_2 = cauchy2d_ext(a0,a0,a1,a1);
a0_2a2_2 = cauchy2d_ext(a0,a0,a2,a2);
a0_2a1a2 = cauchy2d_ext(a0,a0,a1,a2);
a1_3a2 = cauchy2d_ext(a1,a1,a1,a2);
a0a1_2a2 = cauchy2d_ext(a0,a1,a1,a2);
a0_3a2 = cauchy2d_ext(a0,a0,a0,a2);
a0_3a1 = cauchy2d_ext(a0,a0,a0,a1);
a0a1_3 = cauchy2d_ext(a0,a1,a1,a1);
a0a1a2_2 = cauchy2d_ext(a0,a1,a2,a2);
a1_4 = cauchy2d_ext(a1,a1,a1,a1);
a0a2_3 = cauchy2d_ext(a0,a2,a2,a2);
a1_2a2_2 = cauchy2d_ext(a1,a1,a2,a2);

one = zeros(4*N+1); one(1,1) = 1;

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

n1 = (0:N+1)';
tail_terms = abs(lambda(2))*(N+1) + n1*(abs(lambda(1))-abs(lambda(2)));

tail = min(tail_terms);

Z1_0 = (1/tail)*(norm1(c00)+norm1(c01)+norm1(c02));
Z1_1 = (1/tail)*(norm1(c10)+norm1(c11)+norm1(c12));
Z1_2 = (1/tail)*(norm1(c20)+norm1(c21)+norm1(c22));

Z1 = max([sup(Z1_0) sup(Z1_1) sup(Z1_2)]);

disp(['Z1 = ',num2str(Z1)])

%%%%%%%%%%
%%% Z2 %%%
%%%%%%%%%%

A00 = A(i0,i0); A01 = A(i0,i1); A02 = A(i0,i2); 
A10 = A(i1,i0); A11 = A(i1,i1); A12 = A(i1,i2); 
A20 = A(i2,i0); A21 = A(i2,i1); A22 = A(i2,i2);

r_star = 1e-6;

h = infsup(-r_star,r_star);

b0 = a0 + h;
b1 = a1 + h;
b2 = a2 + h;

a = a_par;

b0_3 = cauchy2d(b0,b0,b0);
b1_3 = cauchy2d(b1,b1,b1);
b2_3 = cauchy2d(b2,b2,b2);
b0b1_2 = cauchy2d(b0,b1,b1);
b0b1b2 = cauchy2d(b0,b1,b2);
b0b2_2 = cauchy2d(b0,b2,b2);
b1_2b2 = cauchy2d(b1,b1,b2);
b0_2b1 = cauchy2d(b0,b0,b1);
b0_2b2 = cauchy2d(b0,b0,b2);
b1b2_2 = cauchy2d(b1,b2,b2);

nb0 = norm1(a0) + r_star;
nb1 = norm1(a1) + r_star;
nb2 = norm1(a2) + r_star;

% The following computations only require data of size N-1. Hence, Cauchy
% is used and not cauchy2d_ext for this part. 

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

g0_x0x0 = reshape_data_2d_manifold_matrix2vec(g0_x0x0);
g0_x0x1 = reshape_data_2d_manifold_matrix2vec(g0_x0x1);
g0_x0x2 = reshape_data_2d_manifold_matrix2vec(g0_x0x2);
g0_x1x1 = reshape_data_2d_manifold_matrix2vec(g0_x1x1);
g0_x1x2 = reshape_data_2d_manifold_matrix2vec(g0_x1x2);
g0_x2x2 = reshape_data_2d_manifold_matrix2vec(g0_x2x2);

g1_x0x0 = reshape_data_2d_manifold_matrix2vec(g1_x0x0);
g1_x0x1 = reshape_data_2d_manifold_matrix2vec(g1_x0x1);
g1_x0x2 = reshape_data_2d_manifold_matrix2vec(g1_x0x2);
g1_x1x1 = reshape_data_2d_manifold_matrix2vec(g1_x1x1);
g1_x1x2 = reshape_data_2d_manifold_matrix2vec(g1_x1x2);
g1_x2x2 = reshape_data_2d_manifold_matrix2vec(g1_x2x2);

g2_x0x0 = reshape_data_2d_manifold_matrix2vec(g2_x0x0);
g2_x0x1 = reshape_data_2d_manifold_matrix2vec(g2_x0x1);
g2_x0x2 = reshape_data_2d_manifold_matrix2vec(g2_x0x2);
g2_x1x1 = reshape_data_2d_manifold_matrix2vec(g2_x1x1);
g2_x1x2 = reshape_data_2d_manifold_matrix2vec(g2_x1x2);
g2_x2x2 = reshape_data_2d_manifold_matrix2vec(g2_x2x2);

Z2_0_F = norm1(A00*g0_x0x0 + A01*g1_x0x0 + A02*g2_x0x0) ...
     + 2*norm1(A00*g0_x0x1 + A01*g1_x0x1 + A02*g2_x0x1) ...
     + 2*norm1(A00*g0_x0x2 + A01*g1_x0x2 + A02*g2_x0x2) ...
     +   norm1(A00*g0_x1x1 + A01*g1_x1x1 + A02*g2_x1x1) ...
     + 2*norm1(A00*g0_x1x2 + A01*g1_x1x2 + A02*g2_x1x2) ...
     +   norm1(A00*g0_x2x2 + A01*g1_x2x2 + A02*g2_x2x2); 
 
Z2_0_tail = (1/tail)*(ng0_x0x0 + 2*ng0_x0x1 + 2*ng0_x0x2 + ng0_x1x1 + 2*ng0_x1x2 + ng0_x2x2);

Z2_0 = Z2_0_F + Z2_0_tail;

Z2_1_F = norm1(A10*g0_x0x0 + A11*g1_x0x0 + A12*g2_x0x0) ...
     + 2*norm1(A10*g0_x0x1 + A11*g1_x0x1 + A12*g2_x0x1) ...
     + 2*norm1(A10*g0_x0x2 + A11*g1_x0x2 + A12*g2_x0x2) ...
     +   norm1(A10*g0_x1x1 + A11*g1_x1x1 + A12*g2_x1x1) ...
     + 2*norm1(A10*g0_x1x2 + A11*g1_x1x2 + A12*g2_x1x2) ...
     +   norm1(A10*g0_x2x2 + A11*g1_x2x2 + A12*g2_x2x2); 
 
Z2_1_tail = (1/tail)*(ng1_x0x0 + 2*ng1_x0x1 + 2*ng1_x0x2 + ng1_x1x1 + 2*ng1_x1x2 + ng1_x2x2);

Z2_1 = Z2_1_F + Z2_1_tail;

Z2_2_F = norm1(A20*g0_x0x0 + A21*g1_x0x0 + A22*g2_x0x0) ...
     + 2*norm1(A20*g0_x0x1 + A21*g1_x0x1 + A22*g2_x0x1) ...
     + 2*norm1(A20*g0_x0x2 + A21*g1_x0x2 + A22*g2_x0x2) ...
     +   norm1(A20*g0_x1x1 + A21*g1_x1x1 + A22*g2_x1x1) ...
     + 2*norm1(A20*g0_x1x2 + A21*g1_x1x2 + A22*g2_x1x2) ...
     +   norm1(A20*g0_x2x2 + A21*g1_x2x2 + A22*g2_x2x2); 
 
Z2_2_tail = (1/tail)*(ng2_x0x0 + 2*ng2_x0x1 + 2*ng2_x0x2 + ng2_x1x1 + 2*ng2_x1x2 + ng2_x2x2);

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