function Df = Df_1d_stable_manifold_ex3(a,scaling)

N = (length(a)+3)/3; % a = (a0,a1,a2) with ai = (ai_n)_{n=2}^N \in R^{N-1}

a0 = a(1:N-1); a1 = a(N:2*N-2); a2 = a(2*N-1:3*N-3);

a_par = 0.3; c = 0.7; delta = 9; w = 0.02;

tx = [0.933378938375587;0.358892403648868;0]; % Fixed point
lambda = -1.742392485206276;
v = scaling*[0.743908505928788;0.668260964215875;0.005236268906990];

a0 = [tx(1);v(1);a0]; a1 = [tx(2);v(2);a1]; a2 = [tx(3);v(3);a2];

a0_2 = cauchy([a0 a0]);
a1_2 = cauchy([a1 a1]);
a2_2 = cauchy([a2 a2]);
a0a1 = cauchy([a0 a1]);
a0a2 = cauchy([a0 a2]);

a0_4 = cauchy([a0 a0 a0 a0]);
a0_2a1_2 = cauchy([a0 a0 a1 a1]);
a0_2a2_2 = cauchy([a0 a0 a2 a2]);
a0_2a1a2 = cauchy([a0 a0 a1 a2]);
a1_3a2 = cauchy([a1 a1 a1 a2]);
a0a1_2a2 = cauchy([a0 a1 a1 a2]);
a0_3a2 = cauchy([a0 a0 a0 a2]);
a0_3a1 = cauchy([a0 a0 a0 a1]);
a0a1_3 = cauchy([a0 a1 a1 a1]);
a0a1a2_2 = cauchy([a0 a1 a2 a2]);
a1_4 = cauchy([a1 a1 a1 a1]);
a0a2_3 = cauchy([a0 a2 a2 a2]);
a1_2a2_2 = cauchy([a1 a1 a2 a2]);

one = [1;zeros(N,1)];

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

dphi_00 = intval(zeros(N-1));
dphi_01 = intval(zeros(N-1));
dphi_02 = intval(zeros(N-1));
dphi_10 = intval(zeros(N-1));
dphi_11 = intval(zeros(N-1));
dphi_12 = intval(zeros(N-1));
dphi_20 = intval(zeros(N-1));
dphi_21 = intval(zeros(N-1));
dphi_22 = intval(zeros(N-1));

for n = 2:N
    j = (2:n);
    k = n-j+1;
    dphi_00(n-1,j-1) = c00(k);
    dphi_01(n-1,j-1) = c01(k);
    dphi_02(n-1,j-1) = c02(k);
    dphi_10(n-1,j-1) = c10(k);
    dphi_11(n-1,j-1) = c11(k);
    dphi_12(n-1,j-1) = c12(k);
    dphi_20(n-1,j-1) = c20(k);
    dphi_21(n-1,j-1) = c21(k);
    dphi_22(n-1,j-1) = c22(k);
end

n = (2:N)';

Df_00 = lambda*diag(n) - dphi_00 ;
Df_01 = - dphi_01 ;
Df_02 = - dphi_02 ;
Df_10 = - dphi_10 ;
Df_11 = lambda*diag(n) - dphi_11 ;
Df_12 = - dphi_12 ;
Df_20 = - dphi_20 ;
Df_21 = - dphi_21 ;
Df_22 = lambda*diag(n) - dphi_22 ;

Df = [[Df_00 Df_01 Df_02];[Df_10 Df_11 Df_12];[Df_20 Df_21 Df_22]];

end

