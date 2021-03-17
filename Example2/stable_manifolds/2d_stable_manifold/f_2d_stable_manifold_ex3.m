function f = f_2d_stable_manifold_ex3(a,scaling,fp)

N = round((-3+sqrt(25+8/3*length(a)))/2,0); % a = (a0,a1,a2) with ai = (ai_{n})_{|n|=2}^N \in R^{(N+1)*(N+2)/2-3}

i0 = (1:(N+1)*(N+2)/2-3); i1 = i0+(N+1)*(N+2)/2-3; i2 = i1+(N+1)*(N+2)/2-3;
a0 = a(i0); a1 = a(i1); a2 = a(i2);

a_par = 0.3; c = 0.7; delta = 9; w = 0.02;

if fp == 1
    tx = [0.718092800211620;0.695947361719431;0];
    lambda = [-1.031314539431529;-0.114370869183374];
    v1 = scaling(1)*[0.023249227542643;0.999241674152604;-0.03123379668521];
    v2 = scaling(2)*[0.664951647788332;-0.686110785140320;0.295112345756185];
end

if fp == 2
    tx = [0.998562893876901;-0.053592415248714;0];
    lambda = [-1.994255706055617;-0.187090101867703];
    v1 = scaling(1)*[-0.994185185162860;0.107676008418592;-0.001301850117418];
    v2 = scaling(2)*[0.052670683572220;0.981388690288807;-0.184667370331782];
end
    
a0 = reshape_data_2d_manifold_vec2matrix(a0,tx(1),v1(1),v2(1));
a1 = reshape_data_2d_manifold_vec2matrix(a1,tx(2),v1(2),v2(2));
a2 = reshape_data_2d_manifold_vec2matrix(a2,tx(3),v1(3),v2(3));

a0_3 = cauchy2d(a0,a0,a0);
a0a1_2 = cauchy2d(a0,a1,a1);
a0a2_2 = cauchy2d(a0,a2,a2);
a0_5 = cauchy2d(a0,a0,a0,a0,a0);
a0_3a1_2 = cauchy2d(a0,a0,a0,a1,a1);
a0_3a2_2 = cauchy2d(a0,a0,a0,a2,a2);
a0_3a1a2 = cauchy2d(a0,a0,a0,a1,a2);
a0a1_3a2 = cauchy2d(a0,a1,a1,a1,a2);
a0_2a1_2a2 = cauchy2d(a0,a0,a1,a1,a2);
a0_4a2 = cauchy2d(a0,a0,a0,a0,a2);
a0_2a1 = cauchy2d(a0,a0,a1);
a0_2a2 = cauchy2d(a0,a0,a2);
a0_4a1 = cauchy2d(a0,a0,a0,a0,a1);
a0_2a1_3 = cauchy2d(a0,a0,a1,a1,a1);
a0_2a1a2_2 = cauchy2d(a0,a0,a1,a2,a2);
a1_4a2 = cauchy2d(a1,a1,a1,a1,a2);
a1_3 = cauchy2d(a1,a1,a1);
a0_2a2_3 = cauchy2d(a0,a0,a2,a2,a2);
a1_3a2_2 = cauchy2d(a1,a1,a1,a2,a2);
a0a1_2a2_2 = cauchy2d(a0,a1,a1,a2,a2);

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

f0 = (lambda(1)*matrix_n1 + lambda(2)*matrix_n2).*a0 - phi0;
f1 = (lambda(1)*matrix_n1 + lambda(2)*matrix_n2).*a1 - phi1;
f2 = (lambda(1)*matrix_n1 + lambda(2)*matrix_n2).*a2 - phi2;

f0 = reshape_data_2d_manifold_matrix2vec(f0);
f1 = reshape_data_2d_manifold_matrix2vec(f1);
f2 = reshape_data_2d_manifold_matrix2vec(f2);

f = [f0;f1;f2];

end