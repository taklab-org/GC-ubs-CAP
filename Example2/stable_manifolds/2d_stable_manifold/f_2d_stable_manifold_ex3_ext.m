function [f0,f1,f2] = f_2d_stable_manifold_ex3_ext(a,scaling,fp)

N = round((-3+sqrt(25+8/3*length(a)))/2,0); % a = (a0,a1,a2) with ai = (ai_{n})_{|n|=2}^N \in R^{(N+1)*(N+2)/2-3}

i0 = (1:(N+1)*(N+2)/2-3); i1 = i0+(N+1)*(N+2)/2-3; i2 = i1+(N+1)*(N+2)/2-3;
a0 = a(i0); a1 = a(i1); a2 = a(i2);

a_par = intval('0.3'); c = intval('0.7'); delta = intval(9); w = intval('0.02');

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
    
a0 = reshape_data_2d_manifold_vec2matrix(a0,tx(1),v1(1),v2(1));
a1 = reshape_data_2d_manifold_vec2matrix(a1,tx(2),v1(2),v2(2));
a2 = reshape_data_2d_manifold_vec2matrix(a2,tx(3),v1(3),v2(3));

a0_3 = padmat(cauchy2d_ext(a0,a0,a0),2*N);
a0a1_2 = padmat(cauchy2d_ext(a0,a1,a1),2*N);
a0a2_2 = padmat(cauchy2d_ext(a0,a2,a2),2*N);

a0_5 = cauchy2d_ext(a0,a0,a0,a0,a0);
a0_3a1_2 = cauchy2d_ext(a0,a0,a0,a1,a1);
a0_3a2_2 = cauchy2d_ext(a0,a0,a0,a2,a2);
a0_3a1a2 = cauchy2d_ext(a0,a0,a0,a1,a2);
a0a1_3a2 = cauchy2d_ext(a0,a1,a1,a1,a2);
a0_2a1_2a2 = cauchy2d_ext(a0,a0,a1,a1,a2);
a0_4a2 = cauchy2d_ext(a0,a0,a0,a0,a2);

a0_2a1 = padmat(cauchy2d_ext(a0,a0,a1),2*N);
a0_2a2 = padmat(cauchy2d_ext(a0,a0,a2),2*N);

a0_4a1 = cauchy2d_ext(a0,a0,a0,a0,a1);
a0_2a1_3 = cauchy2d_ext(a0,a0,a1,a1,a1);
a0_2a1a2_2 = cauchy2d_ext(a0,a0,a1,a2,a2);
a1_4a2 = cauchy2d_ext(a1,a1,a1,a1,a2);
a1_3 = padmat(cauchy2d_ext(a1,a1,a1),2*N);
a0_2a2_3 = cauchy2d_ext(a0,a0,a2,a2,a2);
a1_3a2_2 = cauchy2d_ext(a1,a1,a1,a2,a2);
a0a1_2a2_2 = cauchy2d_ext(a0,a1,a1,a2,a2);

a0 = padmat(a0,4*N);
a1 = padmat(a1,4*N);
a2 = padmat(a2,4*N);

phi0 = 3*a0_3 - a0 + a0a1_2 + a0a2_2 - 2*a0_5 - 2*a0_3a1_2 ...
    - (c/delta+2)*a0_3a2_2 - (1+a_par/delta)*a0_3a1a2 - (1/delta)*a0a1_3a2 ...
    + ((a_par+1)/delta)*a0_2a1_2a2 - (w/delta)*a0_4a2;

phi1 = 2*a0_2a1 + a0_2a2 - 2*a0_4a1 - 2*a0_2a1_3 - (c/delta+2)*a0_2a1a2_2 ...
    - (1+a_par/delta)*a0_2a1_2a2 - (1/delta)*a1_4a2 + ((a_par+1)/delta)*a0a1_3a2 ...
    - (w/delta)*a0_3a1a2;

phi2 = (2+c/delta)*a0_2a2 + (1/delta)*a1_3 - ((a_par+1)/delta)*a0a1_2 + (a_par/delta)*a0_2a1 ...
    + (w/delta)*a0_3 - 2*a0_4a2 - 2*a0_2a1_2a2 - (c/delta+2)*a0_2a2_3 ...
    - (1+a_par/delta)*a0_2a1a2_2 - (1/delta)*a1_3a2_2 + ((a_par+1)/delta)*a0a1_2a2_2 ...
    - (w/delta)*a0_3a2_2;

matrix_n1 = repmat((0:N)',1,N+1);
matrix_n2 = repmat((0:N),N+1,1);

LAMBDA = intval(zeros(5*N+1));
LAMBDA(1:N+1,1:N+1) = lambda(1)*matrix_n1 + lambda(2)*matrix_n2;

f0 = LAMBDA.*a0 - phi0;
f1 = LAMBDA.*a1 - phi1;
f2 = LAMBDA.*a2 - phi2;

%f0 = reshape_data_2d_manifold_matrix2vec(f0);
%f1 = reshape_data_2d_manifold_matrix2vec(f1);
%f2 = reshape_data_2d_manifold_matrix2vec(f2);
 
% f = [f0;f1;f2];

end