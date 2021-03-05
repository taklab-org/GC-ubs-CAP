function [f0,f1,f2] = f_1d_stable_manifold_ex3_ext(a,scaling)

N = (length(a)+3)/3; % a = (a0,a1,a2) with ai = (ai_n)_{n=2}^N \in R^{N-1}

a0 = a(1:N-1); a1 = a(N:2*N-2); a2 = a(2*N-1:3*N-3);

a_par = intval('0.3'); c = intval('0.7'); delta = intval(9); w = intval('0.02');

scaling = intval(scaling);

load initial_data_ex3
    
% Stable eigenvalue (includes the rigorous error bound)
tx = ip0;
lambda = tlambda0; 
tv = scaling*tv0; 

a0 = [tx(1);tv(1);a0]; a1 = [tx(2);tv(2);a1]; a2 = [tx(3);tv(3);a2];

n = (2:N)';

a0_3 = cauchy_ext([a0 a0 a0],5*N+1);
a0a1_2 = cauchy_ext([a0 a1 a1],5*N+1);
a0a2_2 = cauchy_ext([a0 a2 a2],5*N+1);
a0_5 = cauchy_ext([a0 a0 a0 a0 a0],5*N+1);
a0_3a1_2 = cauchy_ext([a0 a0 a0 a1 a1],5*N+1);
a0_3a2_2 = cauchy_ext([a0 a0 a0 a2 a2],5*N+1);
a0_3a1a2 = cauchy_ext([a0 a0 a0 a1 a2],5*N+1);
a0a1_3a2 = cauchy_ext([a0 a1 a1 a1 a2],5*N+1);
a0_2a1_2a2 = cauchy_ext([a0 a0 a1 a1 a2],5*N+1);
a0_4a2 = cauchy_ext([a0 a0 a0 a0 a2],5*N+1);
a0_2a1 = cauchy_ext([a0 a0 a1],5*N+1);
a0_2a2 = cauchy_ext([a0 a0 a2],5*N+1);
a0_4a1 = cauchy_ext([a0 a0 a0 a0 a1],5*N+1);
a0_2a1_3 = cauchy_ext([a0 a0 a1 a1 a1],5*N+1);
a0_2a1a2_2 = cauchy_ext([a0 a0 a1 a2 a2],5*N+1);
a1_4a2 = cauchy_ext([a1 a1 a1 a1 a2],5*N+1);
a1_3 =  cauchy_ext([a1 a1 a1],5*N+1);
a0_2a2_3 = cauchy_ext([a0 a0 a2 a2 a2],5*N+1);
a1_3a2_2 = cauchy_ext([a1 a1 a1 a2 a2],5*N+1);
a0a1_2a2_2 = cauchy_ext([a0 a1 a1 a2 a2],5*N+1);

phi0 = 3*a0_3 - [a0;zeros(4*N,1)] + a0a1_2 + a0a2_2 - 2*a0_5 - 2*a0_3a1_2 ...
    - (c/delta+2)*a0_3a2_2 - (1+a_par/delta)*a0_3a1a2 - (1/delta)*a0a1_3a2 ...
    + ((a_par+1)/delta)*a0_2a1_2a2 - (w/delta)*a0_4a2;

phi1 = 2*a0_2a1 + a0_2a2 - 2*a0_4a1 - 2*a0_2a1_3 - (c/delta+2)*a0_2a1a2_2 ...
    - (1+a_par/delta)*a0_2a1_2a2 - (1/delta)*a1_4a2 + ((a_par+1)/delta)*a0a1_3a2 ...
    - (w/delta)*a0_3a1a2;

phi2 = (2+c/delta)*a0_2a2 + (1/delta)*a1_3 - ((a_par+1)/delta)*a0a1_2 + (a_par/delta)*a0_2a1 ...
    + (w/delta)*a0_3 - 2*a0_4a2 - 2*a0_2a1_2a2 - (c/delta+2)*a0_2a2_3 ...
    - (1+a_par/delta)*a0_2a1a2_2 - (1/delta)*a1_3a2_2 + ((a_par+1)/delta)*a0a1_2a2_2 ...
    - (w/delta)*a0_3a2_2;

f0 = [lambda*n.*a0(3:end);zeros(4*N,1)] - phi0(3:end); 
f1 = [lambda*n.*a1(3:end);zeros(4*N,1)] - phi1(3:end); 
f2 = [lambda*n.*a2(3:end);zeros(4*N,1)] - phi2(3:end);

end

