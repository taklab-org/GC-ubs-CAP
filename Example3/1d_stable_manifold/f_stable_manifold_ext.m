function [f1,f2] = f_stable_manifold_ext(a,par)

N = (length(a)+2)/2; % a = (a1,a2) with ai = (ai_n)_{n=2}^N \in R^{N-1}

a1 = a(1:N-1); a2 = a(N:2*N-2);

scaling = intval(par(1)); fp = par(2);

if fp == 0
    
    %%%%%%%%%%%%%%%
    %%%% (0,0) %%%%
    %%%%%%%%%%%%%%%
    
    lambda = -1/4; % Stable eigenvalue
    tx = [0;0]; % Fixed point
    v = scaling*[1;1]; % Stable eigenvector

end

if fp == 1
   
    %%%%%%%%%%%%%%%%%%%%%
    %%%% p_{infty,s} %%%%
    %%%%%%%%%%%%%%%%%%%%%
    
    load initial_data_ex2
    
    % Stable eigenvalue (includes the rigorous error bound)
    lambda = tlambda; %-0.187256681090725;
    
    % Fixed point (includes the rigorous error bound)
    % loads tx = [0.8861081289780320;0.6192579489210105] + error;
    
    % Stable eigenvector (+ rigorous error bound for the proof)
    v = scaling*tv;
    
end

a1 = [tx(1);v(1);a1]; a2 = [tx(2);v(2);a2];

n = (2:N)';

a1_2 = cauchy_ext([a1 a1],11*N+1);
a1_3 = cauchy_ext([a1 a1 a1],13*N+1);
a1_5 = cauchy_ext([a1 a1 a1 a1 a1],13*N+1);
a1_6 = cauchy_ext([a1 a1 a1 a1 a1 a1],11*N+1);
a1_7 = cauchy_ext([a1 a1 a1 a1 a1 a1 a1],13*N+1);
a1_9 = cauchy_ext([a1 a1 a1 a1 a1 a1 a1 a1 a1],13*N+1);
a1_13 = cauchy_ext([a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1],13*N+1);
a1a2_2 = cauchy_ext([a1 a2 a2],13*N+1);
a1a2_4 = cauchy_ext([a1 a2 a2 a2 a2],13*N+1);
a1a2_6 = cauchy_ext([a1 a2 a2 a2 a2 a2 a2],13*N+1);
a1_2a2 = cauchy_ext([a1 a1 a2],11*N+1);
a1_2a2_2 = cauchy_ext([a1 a1 a2 a2],11*N+1);
a1_2a2_3 = cauchy_ext([a1 a1 a2 a2 a2],11*N+1);
a1_2a2_5 = cauchy_ext([a1 a1 a2 a2 a2 a2 a2],11*N+1);
a1_3a2_2 = cauchy_ext([a1 a1 a1 a2 a2],13*N+1);
a1_5a2_2 = cauchy_ext([a1 a1 a1 a1 a1 a2 a2],13*N+1);
a1_4a2 = cauchy_ext([a1 a1 a1 a1 a2],11*N+1);
a1_5a2 = cauchy_ext([a1 a1 a1 a1 a1 a2],13*N+1);
a1_5a2_4 = cauchy_ext([a1 a1 a1 a1 a1 a2 a2 a2 a2],13*N+1);
a1_6a2 = cauchy_ext([a1 a1 a1 a1 a1 a1 a2],11*N+1);
a1_6a2_3 = cauchy_ext([a1 a1 a1 a1 a1 a1 a2 a2 a2],11*N+1);
a1_9a2_2 = cauchy_ext([a1 a1 a1 a1 a1 a1 a1 a1 a1 a2 a2],13*N+1);
a1_10a2 = cauchy_ext([a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a2],11*N+1);
a2_3 = cauchy_ext([a2 a2 a2],11*N+1);

phi1 = (1/4)*a1_2 - (1/4)*a1_6 - (1/4)*[a2;zeros(10*N,1)] + (1/2)*a1_2a2 + (1/12)*a1_4a2 - a1_6a2 + (1/2)*a1_10a2 + (3/4)*a1_2a2_2 - (3/4)*a2_3 - a1_2a2_3 + a1_6a2_3 + (1/2)*a1_2a2_5;
phi2 = -(1/4)*[a1;zeros(12*N,1)] + (1/12)*a1_3 - (1/4)*a1_5 + (1/4)*a1_7 + (5/4)*a1_9 - (3/4)*a1_13 -  2*a1_5a2 + (3/4)*a1a2_2 + (23/12)*a1_3a2_2 + (1/2)*a1_5a2_2 - (5/4)*a1_9a2_2 - (3/4)*a1a2_4 - (1/4)*a1_5a2_4 + (1/4)*a1a2_6;

f1 = [lambda*n.*a1(3:end);zeros(10*N,1)] - phi1(3:end); 
f2 = [lambda*n.*a2(3:end);zeros(12*N,1)] - phi2(3:end);

end

