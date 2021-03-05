function f = f_1d_stable_manifold_ex3(a,scaling)

N = (length(a)+3)/3; % a = (a0,a1,a2) with ai = (ai_n)_{n=2}^N \in R^{N-1}

a0 = a(1:N-1); a1 = a(N:2*N-2); a2 = a(2*N-1:3*N-3);

a_par = 0.3; c = 0.7; delta = 9; w = 0.02;

tx = [0.933378938375587;0.358892403648868;0]; % Fixed point
lambda = -1.742392485206276;
v = scaling*[0.743908505928788;0.668260964215875;0.005236268906990];

a0 = [tx(1);v(1);a0]; a1 = [tx(2);v(2);a1]; a2 = [tx(3);v(3);a2];

n = (2:N)';

a0_3 = cauchy([a0 a0 a0]);
a0a1_2 = cauchy([a0 a1 a1]);
a0a2_2 = cauchy([a0 a2 a2]);
a0_5 = cauchy([a0 a0 a0 a0 a0]);
a0_3a1_2 = cauchy([a0 a0 a0 a1 a1]);
a0_3a2_2 = cauchy([a0 a0 a0 a2 a2]);
a0_3a1a2 = cauchy([a0 a0 a0 a1 a2]);
a0a1_3a2 = cauchy([a0 a1 a1 a1 a2]);
a0_2a1_2a2 = cauchy([a0 a0 a1 a1 a2]);
a0_4a2 = cauchy([a0 a0 a0 a0 a2]);
a0_2a1 = cauchy([a0 a0 a1]);
a0_2a2 = cauchy([a0 a0 a2]);
a0_4a1 = cauchy([a0 a0 a0 a0 a1]);
a0_2a1_3 = cauchy([a0 a0 a1 a1 a1]);
a0_2a1a2_2 = cauchy([a0 a0 a1 a2 a2]);
a1_4a2 = cauchy([a1 a1 a1 a1 a2]);
a1_3 =  cauchy([a1 a1 a1]);
a0_2a2_3 = cauchy([a0 a0 a2 a2 a2]);
a1_3a2_2 = cauchy([a1 a1 a1 a2 a2]);
a0a1_2a2_2 = cauchy([a0 a1 a1 a2 a2]);

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

f = [lambda*n.*a0(3:end) - mid(phi0(3:end)) ; 
     lambda*n.*a1(3:end) - mid(phi1(3:end)) ; 
     lambda*n.*a2(3:end) - mid(phi2(3:end))];

end

