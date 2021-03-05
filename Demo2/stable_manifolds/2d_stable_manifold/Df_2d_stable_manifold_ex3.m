function Df = Df_2d_stable_manifold_ex3(a,scaling,fp)

N = round((-3+sqrt(25+8/3*length(a)))/2,0); % a = (a0,a1,a2) with ai = (ai_{n})_{|n|=2}^N \in R^{(N+1)*(N+2)/2-3}

i0 = (1:(N+1)*(N+2)/2-3); i1 = i0+(N+1)*(N+2)/2-3; i2 = i1+(N+1)*(N+2)/2-3;
a0 = a(i0); a1 = a(i1); a2 = a(i2);

if isintval(a(1,1)) == 1
    
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
    
    dphi_00 = intval(zeros((N+1)^2));
    dphi_01 = intval(zeros((N+1)^2));
    dphi_02 = intval(zeros((N+1)^2));
    dphi_10 = intval(zeros((N+1)^2));
    dphi_11 = intval(zeros((N+1)^2));
    dphi_12 = intval(zeros((N+1)^2));
    dphi_20 = intval(zeros((N+1)^2));
    dphi_21 = intval(zeros((N+1)^2));
    dphi_22 = intval(zeros((N+1)^2));
    
else
    
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
    
    dphi_00 = zeros((N+1)^2);
    dphi_01 = zeros((N+1)^2);
    dphi_02 = zeros((N+1)^2);
    dphi_10 = zeros((N+1)^2);
    dphi_11 = zeros((N+1)^2);
    dphi_12 = zeros((N+1)^2);
    dphi_20 = zeros((N+1)^2);
    dphi_21 = zeros((N+1)^2);
    dphi_22 = zeros((N+1)^2);
    
end

a0 = reshape_data_2d_manifold_vec2matrix(a0,tx(1),v1(1),v2(1));
a1 = reshape_data_2d_manifold_vec2matrix(a1,tx(2),v1(2),v2(2));
a2 = reshape_data_2d_manifold_vec2matrix(a2,tx(3),v1(3),v2(3));

a0_2 = cauchy2d(a0,a0);
a1_2 = cauchy2d(a1,a1);
a2_2 = cauchy2d(a2,a2);
a0a1 = cauchy2d(a0,a1);
a0a2 = cauchy2d(a0,a2);

a0_4 = cauchy2d(a0,a0,a0,a0);
a0_2a1_2 = cauchy2d(a0,a0,a1,a1);
a0_2a2_2 = cauchy2d(a0,a0,a2,a2);
a0_2a1a2 = cauchy2d(a0,a0,a1,a2);
a1_3a2 = cauchy2d(a1,a1,a1,a2);
a0a1_2a2 = cauchy2d(a0,a1,a1,a2);
a0_3a2 = cauchy2d(a0,a0,a0,a2);
a0_3a1 = cauchy2d(a0,a0,a0,a1);
a0a1_3 = cauchy2d(a0,a1,a1,a1);
a0a1a2_2 = cauchy2d(a0,a1,a2,a2);
a1_4 = cauchy2d(a1,a1,a1,a1);
a0a2_3 = cauchy2d(a0,a2,a2,a2);
a1_2a2_2 = cauchy2d(a1,a1,a2,a2);

one = zeros(N+1); one(1,1) = 1;

% phi0 = 3*a0_3 - a0 + a0a1_2 + a0a2_2 - 2*a0_5 - 2*a0_3a1_2 ...
%     - (c/delta+2)*a0_3a2_2 - (1+a_par/delta)*a0_3a1a2 - (1/delta)*a0a1_3a2 ...
%     + ((a_par+1)/delta)*a0_2a1_2a2 - (w/delta)*a0_4a2;

c00 = 9*a0_2 - one + a1_2 + a2_2 - 10*a0_4 - 6*a0_2a1_2 ...
    - 3*(c/delta+2)*a0_2a2_2 - 3*(1+a_par/delta)*a0_2a1a2 - (1/delta)*a1_3a2 ...
    + 2*((a_par+1)/delta)*a0a1_2a2 - 4*(w/delta)*a0_3a2;

c01 = 2*a0a1 - 4*a0_3a1 ...
    - (1+a_par/delta)*a0_3a2 - 3*(1/delta)*a0a1_2a2 ...
    + 2*((a_par+1)/delta)*a0_2a1a2 ;

c02 = 2*a0a2 - 2*(c/delta+2)*a0_3a2 - (1+a_par/delta)*a0_3a1 - (1/delta)*a0a1_3 ...
    + ((a_par+1)/delta)*a0_2a1_2 - (w/delta)*a0_4;

% phi1 = 2*a0_2a1 + a0_2a2 - 2*a0_4a1 - 2*a0_2a1_3 - (c/delta+2)*a0_2a1a2_2 ...
%     - (1+a_par/delta)*a0_2a1_2a2 - (1/delta)*a1_4a2 + ((a_par+1)/delta)*a0a1_3a2 ...
%     - (w/delta)*a0_3a1a2;

c10 = 4*a0a1 + 2*a0a2 - 8*a0_3a1 - 4*a0a1_3 - 2*(c/delta+2)*a0a1a2_2 ...
    - 2*(1+a_par/delta)*a0a1_2a2 + ((a_par+1)/delta)*a1_3a2 ...
    - 3*(w/delta)*a0_2a1a2;

c11 = 2*a0_2 - 2*a0_4 - 6*a0_2a1_2 - (c/delta+2)*a0_2a2_2 ...
    - 2*(1+a_par/delta)*a0_2a1a2 - 4*(1/delta)*a1_3a2 + 3*((a_par+1)/delta)*a0a1_2a2 ...
    - (w/delta)*a0_3a2;

c12 = a0_2 - 2*(c/delta+2)*a0_2a1a2 ...
    - (1+a_par/delta)*a0_2a1_2 - (1/delta)*a1_4 + ((a_par+1)/delta)*a0a1_3 ...
    - (w/delta)*a0_3a1;

% phi2 = (2+c/delta)*a0_2a2 + (1/delta)*a1_3 - ((a_par+1)/delta)*a0a1_2 + (a_par/delta)*a0_2a1 ...
%     + (w/delta)*a0_3 - 2*a0_4a2 - 2*a0_2a1_2a2 - (c/delta+2)*a0_2a2_3 ...
%     - (1+a_par/delta)*a0_2a1a2_2 - (1/delta)*a1_3a2_2 + ((a_par+1)/delta)*a0a1_2a2_2 ...
%     - (w/delta)*a0_3a2_2;

c20 = 2*(2+c/delta)*a0a2 - ((a_par+1)/delta)*a1_2 + 2*(a_par/delta)*a0a1 ...
    + 3*(w/delta)*a0_2 - 8*a0_3a2 - 4*a0a1_2a2 - 2*(c/delta+2)*a0a2_3 ...
    - 2*(1+a_par/delta)*a0a1a2_2 + ((a_par+1)/delta)*a1_2a2_2 ...
    - 3*(w/delta)*a0_2a2_2;

c21 = 3*(1/delta)*a1_2 - 2*((a_par+1)/delta)*a0a1 + (a_par/delta)*a0_2 - 4*a0_2a1a2 ...
    - (1+a_par/delta)*a0_2a2_2 - 3*(1/delta)*a1_2a2_2 + 2*((a_par+1)/delta)*a0a1a2_2;

c22 = (2+c/delta)*a0_2 - 2*a0_4 - 2*a0_2a1_2 - 3*(c/delta+2)*a0_2a2_2 ...
    - 2*(1+a_par/delta)*a0_2a1a2 - 2*(1/delta)*a1_3a2 + 2*((a_par+1)/delta)*a0a1_2a2 ...
    - 2*(w/delta)*a0_3a2;

for m = 0:N
    m1 = (0:m);
    for k = 0:N
        for k1 = 0:k
            dphi_00(k+1+m*(N+1),k1+1+m1*(N+1)) = c00(k-k1+1,m-m1+1);
            dphi_01(k+1+m*(N+1),k1+1+m1*(N+1)) = c01(k-k1+1,m-m1+1);
            dphi_02(k+1+m*(N+1),k1+1+m1*(N+1)) = c02(k-k1+1,m-m1+1);
            dphi_10(k+1+m*(N+1),k1+1+m1*(N+1)) = c10(k-k1+1,m-m1+1);
            dphi_11(k+1+m*(N+1),k1+1+m1*(N+1)) = c11(k-k1+1,m-m1+1);
            dphi_12(k+1+m*(N+1),k1+1+m1*(N+1)) = c12(k-k1+1,m-m1+1);
            dphi_20(k+1+m*(N+1),k1+1+m1*(N+1)) = c20(k-k1+1,m-m1+1);
            dphi_21(k+1+m*(N+1),k1+1+m1*(N+1)) = c21(k-k1+1,m-m1+1);
            dphi_22(k+1+m*(N+1),k1+1+m1*(N+1)) = c22(k-k1+1,m-m1+1);
        end
    end
end

LINEAR = linear_operator_L(lambda,N);

Df_00 = LINEAR - dphi_00 ;
Df_01 = - dphi_01 ;
Df_02 = - dphi_02 ;
Df_10 = - dphi_10 ;
Df_11 = LINEAR - dphi_11 ;
Df_12 = - dphi_12 ;
Df_20 = - dphi_20 ;
Df_21 = - dphi_21 ;
Df_22 = LINEAR - dphi_22 ;

indices2remove = [1;2;N+2;2*(N+1)];

for j = 2:N % Column index
    indices2remove = [indices2remove;(N-j+1:N)'+1+j*(N+1)];
end

Df_00(indices2remove,:) = []; Df_00(:,indices2remove) = [];
Df_01(indices2remove,:) = []; Df_01(:,indices2remove) = [];
Df_02(indices2remove,:) = []; Df_02(:,indices2remove) = [];
Df_10(indices2remove,:) = []; Df_10(:,indices2remove) = [];
Df_11(indices2remove,:) = []; Df_11(:,indices2remove) = [];
Df_12(indices2remove,:) = []; Df_12(:,indices2remove) = [];
Df_20(indices2remove,:) = []; Df_20(:,indices2remove) = [];
Df_21(indices2remove,:) = []; Df_21(:,indices2remove) = [];
Df_22(indices2remove,:) = []; Df_22(:,indices2remove) = [];

Df = [[Df_00 Df_01 Df_02];[Df_10 Df_11 Df_12];[Df_20 Df_21 Df_22]];

end
